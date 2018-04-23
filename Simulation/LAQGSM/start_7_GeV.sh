#!/bin/bash

#SBATCH -D /lustre/nyx/hades/user/parfenov/TMP/
#SBATCH -J reco_run
#SBATCH -p long
#SBATCH --time=7-00:00:00
#SBATCH -a 1-3

#Energy in cms
ecm=7
#Number of events
nev=1000

LAQGSM_FILE=$1

export OUT=/lustre/nyx/hades/user/$USER/mpd_new/DST/MpdDst/${ecm}gev
export OUT_FAM=/lustre/nyx/hades/user/$USER/mpd_new/DST/PicoDst/${ecm}gev
export OUTLOG=$OUT/log
export OUTFAMLOG=$OUT_FAM/log
export ORIGDIR=$PWD

mkdir -p $OUTLOG
mkdir -p $OUTFAMLOG


_log() {

local format='+%Y/%m/%d-%H:%M:%S'
echo [`date $format`] "$@"

}

module use /cvmfs/it.gsi.de/modulefiles/
module load compiler/gcc/6.3.0

. /lustre/nyx/hades/user/$USER/Soft/MPDROOT/mpdroot_dev_040418/SetEnv.sh

export Root=${ROOTSYS}/bin/root

export TMPDIR=/lustre/nyx/hades/user/$USER/TMP/TMP_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

export TMPALL=/lustre/nyx/hades/user/$USER/TMP/

export MACRO_DIR=/lustre/nyx/hades/user/$USER/Soft/MPDROOT/mpdroot_dev_040418/macro/mpd

export BUILD_DIR=/lustre/nyx/hades/user/$USER/Soft/MPDROOT/mpdroot_dev_040418/build

export LAQGSM_DATA_DIR=/lustre/nyx/hades/user/$USER/mpd_new/Simulation/LAQGSM/Data

. ${BUILD_DIR}/config.sh

mkdir -p $TMPDIR

cd $TMPDIR

_log ${TMPDIR}

_log ${INPUTFILE}

_log `ls ${TMPDIR}`

_log ${ROOTSYS_cvmfs}

_log ${SIMPATH}

_log ${ROOTSYS}

_log ${VMCWORKDIR}

_log ${LD_LIBRARY_PATH}

_log ${PATH}

#_log `head -n20 runMC.C`



_log -------


_log `ls $TMPDIR`

cd ${MACRO_DIR}

START_EVENT=`echo "( ${SLURM_ARRAY_TASK_ID} - 1 ) * $nev" | bc`

_log `echo "runMC.C(${LAQGSM_FILE},${TMPDIR}/evetest.root,$START_EVENT,$nev) | tee >> $OUTLOG/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log "`

#$Root -b -q 'runMC.C("'${LAQGSM_FILE}'","'${TMPDIR}'/evetest.root",'`echo "${SLURM_ARRAY_TASK_ID}*$nev" | bc`','$nev')' | tee >> $OUTLOG/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log

_log `echo "reco.C(${TMPDIR}/evetest.root,${TMPDIR}/mpddst.root,$START_EVENT,$nev) | tee >> $OUTLOG/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"`

#$Root -b -q 'reco.C("'${TMPDIR}'/evetest.root","'${TMPDIR}'/mpddst.root",'`echo "${SLURM_ARRAY_TASK_ID}*$nev" | bc`','$nev')' | tee >> $OUTLOG/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log

cd ${TMPDIR}

_log Moving to $OUT/mpddst_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root
#_log `ls .`

mv ${TMPDIR}/mpddst.root $OUT/mpddst_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root

#_log Moving to $OUT_FAM/picodst_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root
#mv ${TMPDIR}/picodst.root $OUT_FAM/picodst_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root



#rm -rf $TMPDIR

#rm -f ${TMPALL}/slurm*out

cd $ORIGDIR
