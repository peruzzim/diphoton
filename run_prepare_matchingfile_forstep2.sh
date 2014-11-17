#!/bin/bash
REMDIR=gg_minitree_data_141106_step1_isogeomcorr
T3HOME=srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/peruzzi
WORKDIR=/scratch/`whoami`/sgejob-${JOB_ID}-$1
source /swshare/ROOT/thisroot.sh
mkdir -p ${WORKDIR}
cd ${WORKDIR} || exit 1
root -q -b -l /shome/peruzzi/CMSSW_5_3_10_patch2/src/DiLeptonAnalysis/NTupleProducer/test/ASAnalysis/diphoton/prepare_matchingfile_forstep2.C+O\(\"/scratch/peruzzi/matchingfile_data_141106_step1_isogeomcorr.root\",\"/scratch/peruzzi/data_step1_141106_isogeomcorr_EXTRA.root\",$1\)
cd ${WORKDIR} || exit 1
for i in matchingtree_*root; do lcg-cp -b -D srmv2 file:///`pwd`/${i} ${T3HOME}/${REMDIR}/matchingtrees/${i}; done
rm -rf ${WORKDIR}
