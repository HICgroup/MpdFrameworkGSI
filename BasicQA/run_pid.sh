#!/bin/bash

#SBATCH -D /tmp
#SBATCH --time=0:15:00

INFILE=$1
OUTFILE=$2

LOG=${OUTFILE%.*}.log

MAINDIR=/lustre/nyx/hades/user/parfenov/mpd_new/BasicQA
SOFTDIR=/lustre/nyx/hades/user/parfenov/Soft/MPDROOT/mpdroot_dev_040418

cd $MAINDIR

. $SOFTDIR/SetEnv.sh
. $SOFTDIR/build/config.sh

./basicqa -i $INFILE -o $OUTFILE --pid &> $LOG
