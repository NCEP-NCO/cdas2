#!/bin/sh

#################################################################
# Script: excdas2_v1_anl.sh.sms
# Purpose: This is the main driver script for the CDAS-2 Analysis 
#          and Forecast.
# Auther: Wesley Ebisuzaki
# Log:   2004-06-09 Initial script from Wesley Ebisuzaki
#        2004-06-14 Modified for production--Julia Zhu
#        2004-06-28 Added restricted data type for the prepbufr
#                   oiqc and pvents files.
#        2004-06-30 The grib table was modified to accomodate the
#                   the R2 analysis data, consequently the new
#                   wgrib executable was used in this job.
################################################################

#########################################################
echo
echo "------------------------------------------------"
echo "excdas2_v1_fcst.sh.sms - CDAS-2 analysis and forcast"
echo "------------------------------------------------"
echo "History: June 20 2004 Original script."
#########################################################
set -xa
export OMP_NUM_THREADS=1

cd $DATA
msg="Job $job_name HAS BEGUN on `hostname`"
postmsg "$msg"

echo "step ############# break ##############################" > $pgmout

msg="CYCLE TIME FOR CDAS-2 ANALYSIS IS $PDY$cyc"
postmsg "$msg"

# (desired) analysis date
export adate=$PDY$cyc
export odate=`$NDATE -6 $adate`

# Definition of input file names for the prepbufr files
date=$adate
prep_000="$COMCDAS2/${RUN}.`echo $date|cut -c1-8`/${model}.t${cyc}z.prepbufr_pre-qc"
date=`$NDATE 12 $adate`
cyc12=`echo $date | cut -c9-10`
prep_p12="$COMGFS/gfs.`echo $date|cut -c1-8`/$cyc12/atmos/gfs.t${cyc12}z.prepbufr_pre-qc"

# for old disk structure WNE
# [ ! -f $prep_p12 ] && prep_p12="$COMGFS/gfs.`echo $date|cut -c1-8`/$cyc12/gfs.t${cyc12}z.prepbufr_pre-qc"

date=`$NDATE 24 $adate`
prep_p24="$COMGFS/gfs.`echo $date|cut -c1-8`/$cyc/atmos/gfs.t${cyc}z.prepbufr_pre-qc"

date=`$NDATE -12 $adate`
prep_m12="$COMCDAS2/${RUN}.`echo $date|cut -c1-8`/${model}.t${cyc12}z.prepbufr_pre-qc"
date=`$NDATE -24 $adate`
prep_m24="$COMCDAS2/${RUN}.`echo $date|cut -c1-8`/${model}.t${cyc}z.prepbufr_pre-qc"

# Definition of input file names for the sst snow and ice grib files
sst="$COMcdas2/${model}.t${cyc}z.sstgrb"
snow="$COMcdas2/${model}.t${cyc}z.snogrb"
ice="$COMcdas2/${model}.t${cyc}z.engicegrb"
obs_precip="$COMcdas2/${model}.t${cyc}z.obs_precip"
mdl_precip="$COMcdas2/${model}.t${cyc}z.mdl_precip"

# Input sigma and sfc analysis files: 
pdy_m06=`echo $odate|cut -c1-8`
cyc_m06=`echo $odate|cut -c9-10`
oldbges="$COMCDAS2/${RUN}.${pdy_m06}/${model}.t${cyc_m06}z.bf06"
oldsges="$COMCDAS2/${RUN}.${pdy_m06}/${model}.t${cyc_m06}z.sf06"
oldsanl="$COMCDAS2/${RUN}.${pdy_m06}/${model}.t${cyc_m06}z.sanl"

# Definition of output file names

sanl="$COMcdas2/${model}.t${cyc}z.sanl"
sfcanl="$COMcdas2/${model}.t${cyc}z.sfcanl"
sges="$COMcdas2/${model}.t${cyc}z.sf06"
bges="$COMcdas2/${model}.t${cyc}z.bf06"

# Definition of archive file names

prepoiqc="$COMarkv/oiqc.anl.$adate.bufr"
prepqm="$COMarkv/pvents.anl.$adate.bufr"
prep_pre="$COMarkv/prepbufr$adate"

sstout="$COMarkv/sstgrb$adate"
snowout="$COMarkv/snogrb$adate"
iceout="$COMarkv/icegrb$adate"

pgbf00="$COMarkv/pgb.anl.$adate.grib"
pgbf06="$COMarkv/pgb.ft06.$adate.grib"
flxf00="$COMarkv/flx.ft00.$adate.grib"
flxf06="$COMarkv/flx.ft06.$adate.grib"
dg3f00="$COMarkv/dg3.ft00.$adate.grib"
dg3f06="$COMarkv/dg3.ft06.$adate.grib"
sgbf00="$COMarkv/sig.anl.$adate.grib"
sgbf06="$COMarkv/sig.ft06.$adate.grib"
znlf00="$COMarkv/znl.ft00.$adate.native"
znlf06="$COMarkv/znl.ft06.$adate.native"

cqb="$COMarkv/cqb.anl.$adate.ascii"
cqe="$COMarkv/cqe.anl.$adate.ascii"
cqt="$COMarkv/cqt.anl.$adate.ascii"

sfcanl_arc="$COMarkv/sfc.anl.$adate.ieee"
bges_arc="$COMarkv/sfc.ft06.$adate.ieee"
sanl_arc="$COMarkv/sig.anl.$adate.ieee"
sges_arc="$COMarkv/sig.ft06.$adate.ieee"

txtout="$COMarkv/stdout1.anl.$adate.ascii"

# misc definitions

restart_step=${restart_step:-1}
PRECIP_SOIL_ADJ=${PRECIP_SOIL_ADJ:-yes}

#
# in-line routines
#

# find date of prepbufr file

prepdate() {
    export FORT11=$1
    $EXECcdas2/cdas2_v1_prepdate
}

# eliminate any FORTnn environment variables

undef_FORTunit() {
   touch fort.0 ; rm fort.*
   export FORT0=abc
   for f in `set | sed 's/=.*//' | grep FORT`
   do
      echo "unset $f"
      unset $f
   done
}

#     check for prepbufr-pre_qc files .. copy to $DATA

hr=-24
for p in $prep_m24 $prep_m12 $prep_000 $prep_p12 $prep_p24
do
   date=`$NDATE $hr $adate`
   if [ ! -f "$p" ] ; then
      
      msg="file: $p is missing"
      postmsg "$msg"
      exit 8
   fi
   ndate=`prepdate $p`
   if [ $date -ne "$ndate" ] ; then
      echo "$p file wrong date wanted $date found $ndate"
      msg="$p file wrong date wanted $date found $ndate"
      postmsg "$msg"
      exit 8
   fi
   hr=$(( $hr + 12 ))
done

cp $prep_m24 $DATA/prepbufr.t-24
cp $prep_m12 $DATA/prepbufr.t-12
cp $prep_000 $DATA/prepbufr.t00
cp $prep_p12 $DATA/prepbufr.t+12
cp $prep_p24 $DATA/prepbufr.t+24

# convert oldsges, oldbges and oldsanl to native format (R8)
set -x
echo "$oldsges $oldbges $oldsanl"

for f in "$oldsges" "$oldbges" "$oldsanl"
do
   if [ ! -s "$f" ] ; then
       echo "missing file $f"
       msg="missing file $f"
       postmsg "$msg"
       exit 8
   fi
done
echo "finished old sanl test"

$EXECcdas2/cdas2_v1_sig2dbl  $oldsges oldsges
$EXECcdas2/cdas2_v1_sig2dbl  $oldbges oldbges
$EXECcdas2/cdas2_v1_sig2dbl  $oldsanl oldsanl

ls -l oldsges

echo "date  c2`echo $adate | cut -c3-10`    washington  " >nmcdate

if [ $restart_step -le 1 ];then
   pgm=prevents
   msg=" `date`  -- $pgm for $adate started "
   postmsg "$msg"

   rm -f bufr_adpupa bufr_aircft bufr_sfcsat
   undef_FORTunit
   export FORT11=prepbufr.t00
   export FORT12=oldsges
   export FORT13=$FIXcdas2/cdas2_v1_ssierr
   export FORT14=nmcdate
   export FORT50=bufr_adpupa
   export FORT51=bufr_aircft
   export FORT52=bufr_sfcsat

   startmsg
   $EXECcdas2/cdas2_v1_$pgm >$pgm.out 2>&1
   err=$?;export err; err_chk
   cat $pgm.out >>$pgmout
fi

#
# 2. cqc
#
if [ $restart_step -le 2 ] ; then

   pgm=cqc
   msg=" `date` -- $pgm started "
   postmsg "$msg"

   rm -f cqe.anl.ascii cqb.anl.ascii cqt.anl.ascii
   undef_FORTunit
   export FORT4=fort.4.cqc
   export FORT14=bufr_adpupa
   export FORT17='prepbufr.t-24'
   export FORT18='prepbufr.t-12'
   export FORT19='prepbufr.t+12'
   export FORT20='prepbufr.t+24'
   export FORT51=bufr_adpupa_cqc_output

   export FORT11=cqc11
   export FORT12=cqe.anl.ascii
   export FORT13=cqc13
   export FORT15=cqb.anl.ascii
   export FORT16=cqc16
   export FORT60=cqt.anl.ascii
   export FORT61=cqc61
   export FORT62=cqc62
   export FORT63=cqc64
   export FORT64=cqc65
   export FORT83=fort.83.cqc
   export FORT84=fort.84.cqc
   export FORT85=fort.85.cqc

   touch cqc11
   touch cqe.anl.ascii
   touch cqc13
   touch cqb.anl.ascii
   touch cqc16
   touch cqt.anl.ascii
   touch cqc61
   touch cqc62
   touch cqc64
   touch cqc65
   touch fort.4.cqc

   startmsg
   $EXECcdas2/cdas2_v1_$pgm  >$pgm.out 2>errfile
   err=$?;export err;err_chk
   cat $pgm.out >>$pgmout

   # not fatal error if cqe.anl.ascii is missing
   [ ! -s cqe.anl.ascii ] && echo >cqe.anl.ascii
   rm FORT4.cqc
fi
#
# 3. acqc
#
if [ $restart_step -le 3 ] ; then

   pgm=acqc
   msg=" `date` -- $pgm started "
   postmsg "$msg"

   rm -f bufr_aircft_acqc_output
   undef_FORTunit
   export FORT14=bufr_aircft
   export FORT15=$FIXcdas2/cdas2_v1_landsea
   export FORT23=$FIXcdas2/cdas2_v1_waypts
   export FORT52=sdmacqc
   export FORT53=sdmstac
   export FORT61=bufr_aircft_acqc_output
   export FORT88=debugout
   export FORT83=fort.83.acqc
   export FORT84=fort.84.acqc
   export FORT85=fort.85.acqc

# note: code from kana has IFLGUS=1
# nersc output has IFLGUS=0
# 1 -> geographical test

   startmsg
   $EXECcdas2/cdas2_v1_$pgm >$pgm.out 2>errfile <<EOF
&INPUT
      DOSPOB = .FALSE., DOACRS=.FALSE., WINDOW=3.00, TIMINC=1.00,
      STCLIM = 41.9,    WAYPIN=.TRUE.,  INIDST= 2,   IFLGUS=0,
      JAMASS = 6*0,     JAWIND= 6*0,
      FWRITE = .TRUE.,  SWRITE=.FALSE., IWRITE=.FALSE., EWRITE=.FALSE.,
      RCPTST=.FALSE.,
/
EOF
   err=$?;export err;err_chk
   cat $pgm.out >>$pgmout
fi

#
#  4.  combbufr
#
#  combines the sfc/sat, adpupa and aircft bufr files
#
if [ $restart_step -le 4 ] ; then

   pgm=combbufr
   msg=" `date` -- $pgm started "
   postmsg "$msg"

   rm -f combbufr_bufr_output
   undef_FORTunit
   export FORT20=bufr_sfcsat
   export FORT21=bufr_adpupa_cqc_output
   export FORT22=bufr_aircft_acqc_output
   export FORT50=combbufr_bufr_output

   startmsg
   echo 3|$EXECcdas2/cdas2_v1_$pgm >$pgm.out 2>errfile
   err=$?;export err;err_chk
   cat $pgm.out >>$pgmout
fi
#
#  5.  oiqc
#
if [ $restart_step -le 5 ] ; then

# changed 11/2012
# for wcoss, bufr library outputs unblocked
#  but oiqc uses the obsolete routine ufbrew
#  which gets the message count of combufr by counting
#  f77 records.  Use cwordsh to add f77 blocking

   $USHcdas2/cwordsh block combbufr_bufr_output combbufr.tmp
   mv combbufr.tmp combbufr_bufr_output

   pgm=oiqc
   msg=" `date` -- $pgm started "
   postmsg "$msg"

   rm -f oiqc.anl.bufr
   undef_FORTunit
   export FORT11=nmcdate
   export FORT14=combbufr_bufr_output
   export FORT17=$FIXcdas2/cdas2_v1_oiqcerr
   export FORT18=obprt.wrk
   export FORT20=tolls.wrk
   export FORT60=obcbt.out
   export FORT61=toss.sfz
   export FORT62=toss.upa
   export FORT63=toss.sat
   export FORT64=toss.smi
   export FORT65=tosslist
   export FORT70=oiqc.anl.bufr
   export FORT81=obogram.out
   export FORT82=obogram.bin
   export FORT83=fort.83.oiqc
   export FORT84=fort.84.oiqc
   export FORT85=fort.85.oiqc

   touch obprt.wrk
   touch tolls.wrk

   startmsg
#  set OpenMP environment
   export OMP_NUM_THREADS=$OMP_NUM_NTHREADS_OIQC
   export OMP_STACKSIZE=$OMP_STACKSIZE_OIQC
   $EXECcdas2/cdas2_v1_$pgm >$pgm.out 2>errfile
   err=$?;export err
#  turn off OpenMP threading
   export OMP_NUM_THREADS=1
   err_chk
   cat $pgm.out >>$pgmout
fi
#
#  6.  ssi
#
if [ $restart_step -le 6 ] ; then

   pgm=ssi
   msg=" `date` -- $pgm started "
   postmsg "$msg"

   # T62 L28 values
   JCAP=62
   LATG=94
   LONF=192
   LEVS=28
   NITER=100

   NLATH=`expr $LATG \/ 2 + 1`
   cat <<EOF >$pgm.parm
&NAMANAL
      JCAP=$JCAP,NLATH=$NLATH,NLON=$LONF,NSIG=$LEVS,
      niter=$NITER,miter=1,
      a=.25,.33,.42,.45,ampdivt=.7,dampdivt=.8,grosst=10.,grossst=10.,
      grossw=10.,grossp=10.,grossq=5.,grosspw=10.,
/
EOF

   rm -f sanl
   undef_FORTunit

   export FORT30=oiqc.anl.bufr
   export FORT35=oldsges
   export FORT36=oldsanl
   export FORT37=oldbges
   export FORT47=$FIXcdas2/cdas2_v1_divterrs2812645
   export FORT48=$FIXcdas2/cdas2_v1_v28newx
   export FORT49=$FIXcdas2/cdas2_v1_eofs28126
   export FORT51=sanl
   export FORT61=fort.61
   export FORT98=scratch3
   export FORT81=fort.81.ssi
   export FORT82=fort.82.ssi
   export FORT83=fort.83.ssi
   export FORT84=fort.84.ssi
   export FORT85=fort.85.ssi

   startmsg
   $EXECcdas2/cdas2_v1_$pgm <$pgm.parm >$pgm.out 2>errfile
   err=$?;
   echo "pgb=$pgm err=$err"
   export err; err_chk
#   err=$?;export err; err_chk

   rm -f fort.61 scratch3

   cat $pgm.out >>$pgmout
fi
ls -l
pwd
#
# 7.  sfc .. creates sfc.anl
#

if [ $restart_step -le 7 ] ; then

   pgm=sfc
   msg=" `date` -- $pgm started "
   postmsg "$msg"

   # DATE CHECK of sanl produced by ssi
   
   export FORT11=sanl
   
   startmsg
   $EXECcdas2/cdas2_v1_sigdate >date.out

   read ndate fhour <date.out
   if [ "$fhour" -ne 0 ] ; then
      msg="date check failed: sanl fhour not zero" 
      postmsg "$msg"
      err=1 ; err_chk
      exit 8
   fi

   if [ "$ndate" != $adate ] ; then
      msg="date check failed: sanl-$ndate vs $adate" 
      postmsg "$msg"
      err=1 ; err_chk
      exit 8
   fi

   #
   #  The following surface file name specification is site dependent
   #

   FNTSFA=$sst
   FNSCVA=$snow
   FNACNA=$ice
   year=`echo $adate | cut -c1-4`
   month=`echo $adate | cut -c5-6`
   day=`echo $adate | cut -c7-8`
   hour=`echo $adate | cut -c9-10`
   undef_FORTunit

   echo "&NAMMAIN"						> $pgm.parm
   echo " IY=$year,IM=$month,ID=$day,IH=$hour,FH=0.,"		>>$pgm.parm
   echo "/"							>>$pgm.parm
   echo "&NAMSFC"						>>$pgm.parm
   echo " FNBGSI='oldbges'"					>>$pgm.parm
   echo " FNBGSO='sfc.anl'"					>>$pgm.parm
   echo " FNOROG='$FIXcdas2/cdas2_v1_orogrd.smth',"		>>$pgm.parm
   echo " FNMASK='$FIXcdas2/cdas2_v1_slmsk',"			>>$pgm.parm
   echo " FNGLAC='$FIXcdas2/cdas2_v1_clim.glacier',"		>>$pgm.parm
   echo " FNMXIC='$FIXcdas2/cdas2_v1_clim.maxice',"		>>$pgm.parm
   echo " FNTSFC='$FIXcdas2/cdas2_v1_clim.sst',"			>>$pgm.parm
   echo " FNWETC='                                  ',"		>>$pgm.parm
   echo " FNSNOC='$FIXcdas2/cdas2_v1_clim.snow',"			>>$pgm.parm
   echo " FNZORC='$FIXcdas2/cdas2_v1_clim.sibrough',"		>>$pgm.parm
   echo " FNALBC='$FIXcdas2/cdas2_v1_clim.matthewalb',"		>>$pgm.parm
   echo " FNAISC='$FIXcdas2/cdas2_v1_clim.ice',"			>>$pgm.parm
   echo " FNPLRC='$FIXcdas2/cdas2_v1_clim.sibresis',"		>>$pgm.parm
   echo " FNTG3C='$FIXcdas2/cdas2_v1_clim.tg3',"			>>$pgm.parm
   echo " FNSCVC='                                   ',"	>>$pgm.parm
   echo " FNSMCC='$FIXcdas2/cdas2_v1_clim.deepsoil',"		>>$pgm.parm
   echo " FNSTCC='                                   ',"	>>$pgm.parm
   echo " FNACNC='                                   ',"	>>$pgm.parm
   echo " FNTSFA='$FNTSFA',"					>>$pgm.parm
   echo " FNAISA='                                   ',"	>>$pgm.parm
   echo " FNSNOA='                                   ',"	>>$pgm.parm
   echo " FNSCVA='$FNSCVA',"					>>$pgm.parm
   echo " FNACNA='$FNACNA',"					>>$pgm.parm
   echo "/"      >>$pgm.parm

   [ -f sfc.anl ] && rm sfc.anl

   startmsg
   $EXECcdas2/cdas2_v1_$pgm <$pgm.parm >$pgm.out 2>errfile
   err=$?; export err; err_chk
   cat $pgm.out >>$pgmout

#  7.5 Soil wetness adjustment
#
# note: converted code to read precip_adj to read single precip wgrib output 5/2004
#

 # find number of the current pentad
   nnumm=`grep -n "$month$day" $FIXcdas2/cdas2_v1_pentads | cut -d':' -f1`

 # find previous pentad
   if [ $nnumm -eq 1 ] ; then
      nnumm=73
      pyyyy=$(( $year - 1 ))
   else
      nnumm=$(( $nnumm - 1 ))
      pyyyy=$year
   fi
 # now, pyyyy, nnumm = year and pentad number of previous pentad

 # find start and finish of the previous pentad
   nams=`sed -n "$nnumm p" $FIXcdas2/cdas2_v1_pentads | cut -d' ' -f1`
   namf=`sed -n "$nnumm p" $FIXcdas2/cdas2_v1_pentads | sed 's/.* //g'`
   pdate=$pyyyy$nams

 # extract obs precip

   $WGRIB -s -4yr $obs_precip | grep ":PRATE:" | grep ":d=$pdate" | \
   $WGRIB -i -s -4yr $obs_precip -ieee -o obs_precip

   if [ ! -s obs_precip ] ; then
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
          msg="The obs precip data has been missing for over the \
               maximum allowed days--$max_nopentad days! Please \
               check with CPC for the proper data"
          postmsg "$msg"
	  echo "No Observation Precipitation data"
          export err=1 ; err_chk
	  exit 9
       else
          echo "No obs precip for $ppdate"
          echo "CDAS2 will continue without obs precip surface analysis"
	  export PRECIP_SOIL_ADJ=no
       fi
   fi

 # extract mdl precip

   $WGRIB -4yr $mdl_precip | sort -t: -k3,3 -k5,5 | \
   $WGRIB -i $mdl_precip -grib -o mdl_precip.grb

   $WGRIB -4yr $mdl_precip | sort -t: -k3,3 -k5,5 | \
   $WGRIB -i $mdl_precip -ieee -o mdl_precip

   if [ ! -s mdl_precip ] ; then
      echo "Missing model precip, should only occur on restart days 1..11"
      msg="Missing model precip, should only occur on restart days 1..11"
      echo "CDAS2 will continue without obs precip surface analysis"
      postmsg "$msg"
      export PRECIP_SOIL_ADJ=no
   fi

   if [ $PRECIP_SOIL_ADJ = yes ] ; then
      pgm=precipadj
      msg="`date` --$pgm started"
      postmsg "$msg"

      undef_FORTunit
      export FORT10=sfc.anl
      export FORT11=obs_precip
      export FORT12=mdl_precip
      export FORT51=sfc.anl.precip_adj
      export FORT52=precipadj.grib

      startmsg
      $EXECcdas2/cdas2_v1_$pgm >$pgm.out 2>errfile
      err=$?; export err; err_chk

      if [ $err -eq 0 ]; then
          mv sfc.anl sfc.anl.noprecip_adj
          mv sfc.anl.precip_adj sfc.anl
      fi
  
      cat $pgm.out >>$pgmout
   fi
fi
pwd
ls -l

#
#  8.  fcst
#
if [ $restart_step -le 8 ] ; then

   pgm=fcst
   msg="`date` --$pgm started"
   postmsg "$msg"


   cp sanl sigit
   cp sfc.anl sfci
   cp sigit sigitdt
   
   #  Forecast parameters

   swhr_gbl=${SHORT_WAVE_INTVL:-1}
   lwhr_gbl=${LONG_WAVE_INTVL:-1}
   
   INCHOUR=6
   PRTHOUR=6
   ENDHOUR=$INCHOUR
   DELTAT=${DELTAT:-1800}
   INTHOUR=${INTHOUR:-0}
   HDIFF_Q_T_RATIO=${HDIFF_Q_T_RATIO:-1.0}
   
   echo "&NAMSMF"                                >$pgm.parm
   echo " NUM(5)=0, "                            >>$pgm.parm
   echo " CON(1)=$DELTAT.,"                      >>$pgm.parm
   echo " CON(3)=0., "                           >>$pgm.parm
   echo " CON(4)=$swhr_gbl.,"                    >>$pgm.parm
   echo " CON(5)=$lwhr_gbl., "                   >>$pgm.parm
   echo " CON(7)=$INCHOUR., "                    >>$pgm.parm
   echo " CON(9)=$PRTHOUR., "                    >>$pgm.parm
   echo " CON(17)=$ENDHOUR.,"                    >>$pgm.parm
   echo " CON(6)=$INTHOUR., "                    >>$pgm.parm
   echo " CON(8)=$HDIFF_Q_T_RATIO,"              >>$pgm.parm
   echo " ICEN2=3, IGEN=195,"                    >>$pgm.parm
   echo "/"                                      >>$pgm.parm

   #
   #  INCHOUR Forecast
   #
   touch sig.ft0 sfc.ft0 flx.ft0 dg3.ft0 znl.ft0
   rm sig.ft* sfc.ft* flx.ft* dg3.ft* znl.ft*
   undef_FORTunit
   
   export FORT11=sigit
   export FORT12=sigitdt
   export FORT14=sfci
   export FORT15=$FIXcdas2/cdas2_v1_co2con
   [ -f heatrate ] && rm heatrate
   export FORT21=heatrate
   export FORT24=$FIXcdas2/cdas2_v1_mtnvar
   export FORT43=$FIXcdas2/cdas2_v1_tune_nmax_1979.ewmrge.sngl.dbl
   export FORT48=$FIXcdas2/cdas2_v1_gcmo3.asc
   export FORT49=$FIXcdas2/cdas2_v1_albaer.snl
   export FORT51=sig.ft06
   export FORT52=sigdt
   export FORT53=sfc.ft06
   export FORT61=znl.ft00.native
   export FORT62=flx.ft00.grib
   export FORT63=flx.ft06.grib
   export FORT64=znl.ft06.native
   export FORT65=dg3.ft00.grib
   export FORT66=dg3.ft06.grib
   export FORT67=ken.ft06.native
   export FORT81=diabanl
   export FORT82=adiages
   export FORT83=fullges
   export FORT92=diagscr
   export FORT98=radscr
   export FORT99=w3out
   touch heatrate

   echo "`date` - Forecast" >>$pgmout

   startmsg
   $EXECcdas2/cdas2_v1_$pgm <$pgm.parm >$pgm.out 2>errfile
   err=$?;export err; err_chk

   cat $pgm.out >>$pgmout
   #
   # append precipadj grib file to the flx file
   #
   if [ $PRECIP_SOIL_ADJ = yes ] ; then
      cat precipadj.grib >>flx.ft06.grib
   fi
fi
#
#  9. Post processing of analysis (PGB)
#
if [ $restart_step -le 9 ] ; then

   pgm=pgb
   msg="`date` --$pgm started"
   postmsg "$msg"

   echo "&NAMPGB" >$pgm.parm
   echo "   ICEN2=3,">>$pgm.parm
   echo "   IGEN=195,">>$pgm.parm
   echo "/" >>$pgm.parm

   [ -f pgb.anl.grib ] && rm pgb.anl.grib 
   [ -f pgb.ft06.grib ] && rm pgb.ft06.grib

   undef_FORTunit
   export FORT11=sanl
   export FORT51=pgb.anl.grib

   startmsg
   $EXECcdas2/cdas2_v1_$pgm <$pgm.parm >$pgm.out 2>errfile
   err=$?;export err; err_chk 
   cat $pgm.out >>$pgmout

   undef_FORTunit
   export FORT11=sig.ft06
   export FORT51=pgb.ft06.grib

   startmsg
   $EXECcdas2/cdas2_v1_$pgm <$pgm.parm >$pgm.out 2>errfile
   err=$?;export err; err_chk 
   cat $pgm.out >>$pgmout

fi

#
# 10. Post processing of (SGB)
#
set -x

if [ $restart_step -le 10 ] ; then

   pgm=sgb
   msg="`date` --$pgm started"
   postmsg "$msg"

   echo "&NAMSGB" >$pgm.parm
   echo "   ICEN2=3,">>$pgm.parm
   echo "   IGEN=195,">>$pgm.parm
   echo "/" >>$pgm.parm

   [ -f sgb.anl.grib  ] && rm sgb.anl.grib 
   [ -f sgb.ft06.grib  ] && rm sgb.ft06.grib 

   undef_FORTunit
   export FORT11=sanl
   export FORT51=sig.anl.grib

   startmsg
   $EXECcdas2/cdas2_v1_$pgm <$pgm.parm >$pgm.out 2>errfile
   err=$?;export err; err_chk
   cat $pgm.out >>$pgmout

   undef_FORTunit
   export FORT11=sig.ft06
   export FORT51=sig.ft06.grib

   startmsg
   $EXECcdas2/cdas2_v1_$pgm <$pgm.parm >$pgm.out 2>errfile
   err=$?;export err; err_chk
   cat $pgm.out >>$pgmout
fi

#
# 11. Postvents
#
# #if [ $restart_step -le 11 ] ; then
if [ "$POSTVENT" = "YES" ] ; then

   pgm=postvents
   msg="`date` --$pgm started"
   postmsg "$msg"

   [ -f pvents.anl.bufr ] && rm pvents.anl.bufr 

   undef_FORTunit
   export FORT20=oiqc.anl.bufr
   export FORT21=sanl
   export FORT51=GFIT.anl.jwpk
   export FORT52=AFIT.anl.jwpk
   export FORT53=dummy
   export FORT50=pvents.anl.bufr
   export FORT54=dummy2
  
   startmsg
   $EXECcdas2/cdas2_v1_$pgm >$pgm.out 2>errfile
   err=$?;export err
   err_chk
   cat $pgm.out >>$pgmout
fi

#
# 12. copy results to $COMcdas2
#

echo "SENDCOM=$SENDCOM COMcdas2=$COMcdas2 COMarkv=$COMarkv"

if [ $restart_step -le 12 -a $SENDCOM = YES ] ; then

   # files to copy to $COMcdas2

   $EXECcdas2/cdas2_v1_sig2sngl sfc.anl $sfcanl
   $EXECcdas2/cdas2_v1_sig2sngl sanl $sanl
   $EXECcdas2/cdas2_v1_sig2sngl sig.ft06 $sges
   $EXECcdas2/cdas2_v1_sig2sngl sfc.ft06 $bges

   # Alert $sanl
   if [ $SENDDBN = YES ]; then
      $DBNROOT/bin/dbn_alert MODEL CDAS2_SA $job $sanl
   fi

   # files to copy to $COMARC

   cp $prep_000 $prep_pre
   [ -f  $prep_pre.gz ] && rm $prep_pre.gz
   gzip $prep_pre

   cp oiqc.anl.bufr $prepoiqc
   [ -f $prepoiqc.gz ] && rm $prepoiqc.gz
   gzip $prepoiqc

   if [ "$POSTVENT" = 'YES' ] ; then
      cp pvents.anl.bufr $prepqm
      [ -f $prepqm.gz ] && rm $prepqm.gz
      gzip $prepqm
   fi
   cp $sst $sstout
   cp $snow $snowout
   cp $ice $iceout

   cp pgb.anl.grib  $pgbf00
   cp pgb.ft06.grib $pgbf06
   cp flx.ft00.grib $flxf00
   cp flx.ft06.grib $flxf06
   cp dg3.ft00.grib $dg3f00
   cp dg3.ft06.grib $dg3f06
   cp sig.anl.grib $sgbf00
   cp sig.ft06.grib $sgbf06
   cp znl.ft00.native $znlf00
   cp znl.ft06.native $znlf06

   cp cqb.anl.ascii $cqb
   [ -f $cqb.gz ] && rm $cqb.gz
   gzip $cqb

   cp cqe.anl.ascii $cqe
   [ -f $cqe.gz ] && rm $cqe.gz
   gzip $cqe 

   cp cqt.anl.ascii $cqt
   [ -f $cqt.gz ] && rm $cqt.gz
   gzip $cqt

   cp $sfcanl $sfcanl_arc
   cp $sanl $sanl_arc
   cp $sges $sges_arc
   cp $bges $bges_arc
   cat $pgmout >stdout1.anl.ascii
   cp stdout1.anl.ascii $txtout
   [ -f $txtout.gz ] && rm $txtout.gz
   gzip $txtout
   if [ "$CHGRP_RSTPROD" = 'YES' ]; then
       for file in "$prepoiqc.gz" "$prepqm.gz" "$prep_pre.gz"
       do
          if [ -s "$file" ] ; then
             chgrp rstprod "$file"
             errch=$?
             if [ $errch -eq 0 ]; then
                chmod 640 $file
             else
                cp /dev/null $file
                warning=yes
             fi 
          fi
       done
   fi
fi

