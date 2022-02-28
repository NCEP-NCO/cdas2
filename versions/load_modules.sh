#!/bin/sh
#
# load build.ver or run.ver
#
# . load_modules.sh (name of file(
#
set +x
in=$1
while read line
do
   $line
   n=`echo $line | egrep -c ' (gfs|cdas|cdas2|rcdas|nam|pcpanl|omb)_ver='`
   if [ "$n" -eq 0 ] ; then
      line=`echo $line | sed -e 's/PrgEnv_intel/PrgEnv-intel/' -e 's/cray_mpich/cray-mpich/' -e 's/cray_pals/cray-pals/'`
      cmd=`echo $line | sed -e 's/export/module load/' -e 's,_ver=v,/,'`
      # echo "cmd=$cmd"
   $cmd
   fi
done  <$in
module list
