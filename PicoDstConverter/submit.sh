#!/bin/bash

FILE_LIST=$1
INFILEHIST=$2
DCAFILE=$3
MOD=1

cd /lustre/nyx/hades/user/parfenov/mpd_new/PicoDstConverter/

i=0
while read INFILE; do
  #~ i=$((i+1))
  #~ MOD=$((i%100))
  #~ if [ $MOD == 0 ]
  #~ then
    #~ sleep 10m
  #~ fi
  OUTFILE=/lustre/nyx/hades/user/parfenov/mpd_new/PicoDstConverter/TMP/${INFILE##*/}
  OUTFILE=${OUTFILE%.*}_reduced.root
  sbatch run.sh $INFILEHIST $INFILE $OUTFILE $DCAFILE
  #sleep 1
done <${FILE_LIST}
