
set -x
source ../versions/build.ver
module reset
module use `pwd`
module load cdas2.module.lua
module list

make clean
make
