#!/bin/bash

FILE_LIST=$1
MOD=1

MAINDIR=/lustre/nyx/hades/user/parfenov/mpd_new/BasicQA

cd $MAINDIR

i=0
while read INFILE; do
  #~ i=$((i+1))
  #~ MOD=$((i%100))
  #~ if [ $MOD == 0 ]
  #~ then
    #~ sleep 10m
  #~ fi
  OUTFILE=$MAINDIR/TMP/${INFILE##*/}
  OUTFILE=${OUTFILE%.*}_qa_pid.root
  sbatch run_pid.sh $INFILE $OUTFILE
  #sleep 1
done <${FILE_LIST}
