#!/bin/bash

#SBATCH -D /tmp
#SBATCH --time=0:15:00

INFILE=$1
OUTFILE=$2

LOG=${OUTFILE%.*}.log

cd /lustre/nyx/hades/user/parfenov/mpd_new/get-dca/

. /lustre/nyx/hades/user/parfenov/Soft/MPDROOT/mpdroot_dev_040418/build/config.sh

./get_dca -i $INFILE -o $OUTFILE
