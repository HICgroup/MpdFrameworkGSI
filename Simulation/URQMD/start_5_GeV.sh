#!/bin/bash

#SBATCH -D /lustre/nyx/hades/user/parfenov/TMP/
#SBATCH -J reco_run
#SBATCH -p long
#SBATCH --time=7-00:00:00
#SBATCH -a 1-1000

#Energy in cms
ecm=5
#Number of events
nev=1000

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

# ROOTSYS_cvmfs=/cvmfs/fairroot.gsi.de/fairsoft/may16_root5

# export ROOTSYS=${ROOTSYS_cvmfs}

# export PATH=${PATH}:${ROOTSYS}/bin

# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ROOTSYS}/lib/root:${ROOTSYS}/lib

export Root=${ROOTSYS}/bin/root

export TMPDIR=/lustre/nyx/hades/user/$USER/TMP/TMP_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

export TMPALL=/lustre/nyx/hades/user/$USER/TMP/



export MACRO_DIR=/lustre/nyx/hades/user/$USER/Soft/MPDROOT/mpdroot_dev_040418/macro/mpd

export BUILD_DIR=/lustre/nyx/hades/user/$USER/Soft/MPDROOT/mpdroot_dev_040418/build

export URQMD_DIR=/lustre/nyx/hades/user/$USER/Soft/urqmd-3.4

#export URQMD_DIR_KIREEV=/hera/fopi/kireyeu/bin/urqmd33p2

. ${BUILD_DIR}/config.sh

export INPUTFILE=inputfile

mkdir -p $TMPDIR

cd $TMPDIR

# cp ${MACRO_DIR}/runMC.C ${TMPDIR}/

# cp ${MACRO_DIR}/reco.C ${TMPDIR}/

# cp ${MACRO_DIR}/mpdloadlibs.C ${TMPDIR}

if [[ -f /lustre/nyx/hades/user/$USER/mpd_new/FAMTree/FAMTree.so && -f /lustre/nyx/hades/user/$USER/mpd_new/FAMTree/FAMTree ]]; then
	cp /lustre/nyx/hades/user/$USER/mpd_new/FAMTree/FAMTree.so ${TMPDIR}/
	cp /lustre/nyx/hades/user/$USER/mpd_new/FAMTree/FAMTree ${TMPDIR}/
else
	
	cp /lustre/nyx/hades/user/$USER/mpd_new/FAMTree/FAMTree.LinkDef.h ${TMPDIR}/
	cp /lustre/nyx/hades/user/$USER/mpd_new/FAMTree/FAMTree.h ${TMPDIR}/
	cp /lustre/nyx/hades/user/$USER/mpd_new/FAMTree/FAMTree.cxx ${TMPDIR}/
	cp /lustre/nyx/hades/user/$USER/mpd_new/FAMTree/Makefile ${TMPDIR}/

	make
	rm Makefile FAMTree.h FAMTree.cxx FAMTree.LinkDef.h

fi

cp ${URQMD_DIR}/runqmd.bash ${TMPDIR}/

cp ${URQMD_DIR}/urqmd.x86_64 ${TMPDIR}/

cp ${URQMD_DIR}/$INPUTFILE ${TMPDIR}/inputfile

sed -e "s|energyincms|$ecm|" -i inputfile
sed -e "s|numberofevents|$nev|" -i inputfile
sed -e "s|randomrandom|`shuf -i 1-1000000 -n 1`|" -i inputfile

_log ${TMPDIR}

_log ${INPUTFILE}

_log `ls ${TMPDIR}`

_log ${ROOTSYS_cvmfs}

_log ${SIMPATH}

_log ${ROOTSYS}

_log ${VMCWORKDIR}

_log ${LD_LIBRARY_PATH}

_log ${PATH}

_log `head -n20 runMC.C`



_log -------

#cat inputfile >> $OUT/JOB_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log
cat inputfile >> $OUTLOG/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log

_log -------



. $TMPDIR/runqmd.bash 2>>$OUTLOG/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log

#child=$!

#_log URQMD PID $child

#wait $child

#state=$?

#_log Finishing with $state

_log -------

_log `ls $TMPDIR`

cd ${MACRO_DIR}

$Root -b -q 'runMC.C("'${TMPDIR}'/test.f14","'${TMPDIR}'/evetest.root",0,'$nev')' | tee >> $OUTLOG/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log

$Root -b -q 'reco.C("'${TMPDIR}'/evetest.root","'${TMPDIR}'/mpddst.root",0,'$nev')' | tee >> $OUTLOG/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log

cd ${TMPDIR}

./FAMTree picodst.root mpddst.root | tee >> $OUTFAMLOG/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log

#$Root -b -q 'create_reduced_tree.C("./mpddst.root","./mpddst_reduced.root")' | tee >> $OUT/JOB_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log

_log Moving to $OUT/mpddst_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root
#_log `ls .`

mv ${TMPDIR}/mpddst.root $OUT/mpddst_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root

_log Moving to $OUT_FAM/picodst_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root
mv ${TMPDIR}/picodst.root $OUT_FAM/picodst_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root



rm -rf $TMPDIR

rm -f ${TMPALL}/slurm*out

cd $ORIGDIR
