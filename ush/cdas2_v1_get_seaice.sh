#!/bin/sh
#                                            11/2019 Wesley Ebisuzaki
#
# get sea ice concentration
#
# variable: COM_SICE             
#           COMcdas2             output of sea ice
#           PDY                  YYYYMMDD
#           cyc                  00/06/12/18
#
# INPUT:    $COM_SICE/seaice.t00z.grb
#
# OUTPUT:   $COMcdas2/$ice             (grib1)
#           $COMcdas2/$ice.grib2       (grib2)
#
# routines used
#          $WGRIB, $NDATE
#
# modules used
#          grib_util, prod_util
#
# returns 0 if ok
# returns 1 if failure
#
# v1.0 11/25/2019               Initial version Wesley Ebisuzaki
# 

set -x

export in_ice=seaice.t00z.grb

if [ "$NDATE" = '' -o "$WGRIB" = '' ] ; then
   echo "modules not loaded"
   exit 1
fi

# if file already exists, return
if [ -f $COMcdas2/$ice ] ; then
   n=`$WGRIB $COMcdas2/$ice -s | grep -c ":ICEC:"`
   [ "$n" -eq 1 ] && exit 0
fi

for i in 0 1 2 3 4 5 6 7 8 9 10 11
do
   hour=`expr $i \* 24`
   date=`$NDATE -$hour ${PDY}${cyc} | cut -c1-8`
   in=$COM_SICE.$date/$in_ice
   [ ! -f $in ] && continue

   # check the date code
   ice_date=`$WGRIB -4yr $in | head -n 1 | cut -f3 -d: | cut -c3-`
   [ "$ice_date" = "" ] && continue
   if [ "$ice_date" -gt "$PDY$cyc" ] ; then
      echo "ice_date ($ice_date) is not good"
      continue
   fi

   n=`$WGRIB $in | grep -c ':ICEC:'`
   [ "$n" -eq 0 ] && continue

   cp $in $COMcdas2/$ice
   $CNVGRIB -g12 -p40 $in $COMcdas2/${ice}.grib2
   exit 0

done

echo "FAILED: ICEC not found"
# failed
exit 1
