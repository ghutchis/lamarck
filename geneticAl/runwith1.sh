#!/bin/sh
#PBS -N OneLastJob
#PBS -A ucche002c
#PBS -r n
#PBS -j oe
##PBS -m a
#PBS -M baoilleach@gmail.com
#PBS -l nodes=8:ppn=12,walltime=0:6:00
##PBS -q  DevQ

cd $PBS_O_WORKDIR

here=/ichec/work/ucche002c
ob=$here/Tools/ob-gccinstall
repo=$here/molwire/fromWork
export LD_LIBRARY_PATH=$ob/lib
export PYTHONPATH=$ob/lib64/python2.6/site-packages:$repo/geneticAl:$here/Tools/cclibtrunk/src

module load gaussian/09b01
ulimit -Ss 1048576
ulimit -Sl 524288
ulimit -c 0

module load taskfarm2
taskfarm tasks



