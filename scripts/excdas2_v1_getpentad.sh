#!/bin/sh
############################################################
# Script -- cdas2_v1_getpentad.sh.sms
# Purpose: to retrieve the pentad data from CPC workstation
#          and convert it to GRIB format
# History: 2005-02-01 First implementation
############################################################

set -x
export OMP_NUM_THREADS=1
cd $DATA

#########################################################
# Start retrieving the pentad file (in grads format)from 
# the /dcom/us00700003/pentad.
#########################################################
echo "PDY=$PDY"
if [ "$PDY" = '' ] ; then
  echo '$PDY is undefined'
  exit 8
fi

date=`echo $PDY | cut -c5-8`
# the pentad file generated by CPC on Jan 01 contains the previous year's data
# so the filename is appended with the previous year. On Jan 06, the script 
# needs to bring over that file. A new file will be generated on the 6th with
# the new year appended.
if [ "$date" = "0106" ]; then
   yyyy=`echo $PDYm6 |cut -c1-4` 
else
   yyyy=`echo $PDYm5 |cut -c1-4`
fi
ls -l $DCOM_DIR

remote_file=${DCOM_DIR}/cmap_pen_rt_v0011_out.lnx.$yyyy
local_file=pingrain.$yyyy
cp $remote_file $local_file

if [ $? -ne 0 ] || [ ! -s $local_file ]
then
  echo "Copy failed $remote_file to `pwd`/$local_file"
  exit 8
else
  echo "The pentad data is updated"
fi

###################################
# Convert the file into GRIB format
###################################
export pgm=cdas2_v1_gribpentad

$EXECcdas2/cdas2_v1_gribpentad $local_file pingrain.grb $yyyy


#####################################
# convert data to T62 gaussian grid
#####################################
[ -f pingrainT62.grb ] && rm pingrainT62.grb
$COPYGB -g98 -i3 -x pingrain.grb pingrainT62.grb

[ "$?" -ne 0 ] && exit 9

if [ $SENDCOM = "YES" ]
then
   cp pingrainT62.grb $COMOUT/.
   cp pingrainT62.grb $COMOUT/pingrainT62.grb.$yyyy
fi

msg="$job HAS COMPLETED NORMALLY."
echo $msg

postmsg "$msg"

############## END OF SCRIPT #######################

