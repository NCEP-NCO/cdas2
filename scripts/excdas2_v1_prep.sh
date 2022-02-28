#!/bin/sh

#############################################################
# Script: excdas2_v1_prep.sh.sms
# Purpose: to obtain the observation and model precipitation
#          analysis from CDAS and DCOM
# Auther: Wesley Ebisuzaki
# Log:   2004-06-09 Initial script from Wesley Ebisuzaki
#        2004-06-10 Modified for production--Julia Zhu
#############################################################
set -x
export OMP_NUM_THREADS=1
cd $DATA

###################################################################

# get needed files from CDAS and GFS: prepbufr, snow, sst, ice
#   for day+1  (R2 needs +/- 24hr, +/- 12hr and 00hr prepbufr files
###################################################################

ierr=0

export snow=${RUN}.t${cyc}z.snogrb
$USHcdas2/cdas2_v1_get_afwa_snow.sh
if [ $? -ne 0 ] ; then
   echo "SNOW not created"
   msg="snow file not created"
   postmsg "$msg"
   exit 9
fi

export ice=${RUN}.t${cyc}z.engicegrb
$USHcdas2/cdas2_v1_get_seaice.sh
if [ $? -ne 0 ] ; then
   echo "Sea Ice not created"
   msg="Sea Ice file not created"
   postmsg "$msg"
   exit 9
fi

export sst=${RUN}.t${cyc}z.sstgrb
$USHcdas2/cdas2_v1_get_sst.sh
if [ $? -ne 0 ] ; then
   echo "SST not created"
   msg="SST file not created"
   postmsg "$msg"
   exit 9
fi

cpreq $COMobsproc/cdas.${PDY}/cdas.t${cyc}z.prepbufr_pre-qc $COMcdas2/${RUN}.t${cyc}z.prepbufr_pre-qc

if [ "$CHGRP_RSTPROD" = 'YES' ]; then
   chgrp rstprod $COMcdas2/${RUN}.t${cyc}z.prepbufr_pre-qc
   errch=$?
   if [ $errch -eq 0 ]; then
      chmod 640 $COMcdas2/${RUN}.t${cyc}z.prepbufr_pre-qc
   else
      cp /dev/null $COMcdas2/${RUN}.t${cyc}z.prepbufr_pre-qc
      warning=yes
   fi
fi


#####################################
# get observerd (pentad) precip file
#####################################
# find number of the current pentad
mmdd=`echo $PDY | cut -c5-8`
year=`echo $PDY | cut -c1-4`
nnumm=`grep -n "$mmdd" $FIXcdas2/cdas2_v1_pentads | cut -d':' -f1`
if [ "$nnumm" = "" ] ; then
   echo "problem in $0 - missing $FIXcdas2/cdas2_v1_pentads?"
   msg="missing $FIXcdas2/cdas2_v1_pentads? nnumm was not calculated"
   postmsg "$msg"
   exit 9
fi

# find previous pentad
if [ $nnumm -eq 1 ] ; then
   nnumm=73
   pyyyy=`expr $year - 1`
else
   nnumm="`expr $nnumm - 1`"
   pyyyy=$year
fi

# find start and finish of the previous pentad
nams=`sed -n "$nnumm p" $FIXcdas2/cdas2_v1_pentads | cut -d' ' -f1`
namf=`sed -n "$nnumm p" $FIXcdas2/cdas2_v1_pentads | sed 's/.* //g'`
ppdate=$pyyyy$nams

# now, ppdate = 1st day of previous pentad
#      nnum = number of previous pentad

#changed line 75 from egrep ":d=${ppdate}..:PRATE:sfc:0-5d ave:" 
cp $COMPCP/pingrainT62.grb.$pyyyy pingraint62

$WGRIB -s -4yr pingraint62 |  \
   egrep ":d=${ppdate}..:PRATE:sfc:0-(5|6)d ave:" \
  | $WGRIB -s  pingraint62 -grib -i -o $DATA/cdas2.t${cyc}z.obs_precip.tmp
$COPYGB -g98 -i3 -x $DATA/cdas2.t${cyc}z.obs_precip.tmp \
       $COMcdas2/cdas2.t${cyc}z.obs_precip
rm $DATA/cdas2.t${cyc}z.obs_precip.tmp

if [ ! -s $COMcdas2/cdas2.t${cyc}z.obs_precip ] ; then

  # Check to see for how many days since the pentad data has been missing
  # The job will be aborted after it reaches the max_nopentad day

  iday=0
  nday=$pyyyy$namf

  while [ $nday -le $PDY ]
  do
     iday=`expr $iday + 1`
     nday=`$NDATE +24 ${nday}00 | cut -c1-8`
  done

  if [ $iday -ge $max_nopentad ]
  then
     msg="The obs precip data has been missing for over the maximum \
          allowed days -- $max_nopentad days, please check with CPC's \ 
          Pingping Xie @571-241-7019 or email him at pingping.xie@noaa.gov"

     postmsg "$msg"
     echo "No Observation Precipitation data"
     exit 9
  else
     echo "No obs precip for $ppdate" 
     echo "CDAS2 will continue without obs precip surface analysis"
  fi
fi

#######################################
# get model precip for previous pentad
#######################################

odate=$pyyyy${nams}00
ldate=$pyyyy${namf}18
cdate=$odate

[ -f $COMcdas2/cdas2.t${cyc}z.mdl_precip ] && rm $COMcdas2/cdas2.t${cyc}z.mdl_precip
while [ $cdate -le $ldate ]
do
   eval filename="$COMarkv/flx.ft06.$cdate.grib"
   if [ ! -s $filename ] ; then
      echo "missing file $filename" 
      msg="missing file $filename" 
      postmsg "$msg"
      exit 9
   fi
   $WGRIB $filename | egrep "(PRATE|RUNOF|LAND)" | \
      $WGRIB -i -grib $filename -append -o $COMcdas2/cdas2.t${cyc}z.mdl_precip
   cdate=`$NDATE +6 $cdate` 
done
echo "done getting the model precipitation"

echo "Complete the CDAS2 PREP job"
exit 0
