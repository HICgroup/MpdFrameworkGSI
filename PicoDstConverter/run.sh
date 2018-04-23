#!/bin/bash

#SBATCH -D /tmp
#SBATCH --time=0:15:00

INFILEHIST=$1
INFILE=$2
OUTFILE=$3
DCAFILE=$4

LOG=${OUTFILE%.*}.log

cd /lustre/nyx/hades/user/parfenov/mpd_new/PicoDstConverter

. /lustre/nyx/hades/user/parfenov/Soft/MPDROOT/mpdroot_dev_040418/SetEnv.sh
. /lustre/nyx/hades/user/parfenov/Soft/MPDROOT/mpdroot_dev_040418/build/config.sh

./reducedTreeCreator -i $INFILE -o $OUTFILE -dca $DCAFILE -centrality $INFILEHIST &> $LOG
