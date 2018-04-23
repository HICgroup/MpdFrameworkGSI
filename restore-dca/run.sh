#!/bin/bash

#SBATCH -D /tmp
#SBATCH --time=0:25:00

INFILE=$1
OUTFILE=$2

LOG=${OUTFILE%.*}.log

cd /lustre/nyx/hades/user/parfenov/mpd_new/restore-dca/

source /lustre/nyx/hades/user/parfenov/Soft/MPDROOT/mpdroot_dev_040418/build/config.sh

if [ -f restore_dca ]; then
	./restore_dca -i $INFILE -o $OUTFILE > $LOG
else
	make
	./restore_dca -i $INFILE -o $OUTFILE > $LOG
fi
