#!/bin/bash

REMDIR=gg_minitree_data_141106_step1_isogeomcorr
FILE1=matchingfile_data_141106_step1_isogeomcorr.root
FILE2=data_step1_141106_isogeomcorr_EXTRA.root

T3DCAPHOME=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/peruzzi/

mkdir -p /scratch/peruzzi
cd /scratch/peruzzi || exit 1
rm ${FILE1} ${FILE2}
dccp ${T3DCAPHOME}/${REMDIR}/${FILE1} .
dccp ${T3DCAPHOME}/${REMDIR}/${FILE2} .

