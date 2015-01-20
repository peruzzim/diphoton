#!/bin/bash

source /swshare/ROOT/thisroot.sh

T3DCAPHOME=dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat//store/user/peruzzi
T3HOME=srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/peruzzi

for sample in mu1 mu1o2 mu2
do
for file in virt int LO real
do

#qsub -q all.q ~/shape_studies/theory_marco/run_gosam_converter.sh ${T3DCAPHOME}/GoSam_diphoton_ntuples/pp_aajj/${sample}/${file}.root ${T3HOME}/GoSam_diphoton_ntuples/pp_aajj/${sample}/${file}_converted.root
#qsub -q all.q ~/shape_studies/theory_marco/run_genlevelacceptance.sh ${T3DCAPHOME}/GoSam_diphoton_ntuples/pp_aajj/${sample}/${file}_converted.root ${T3HOME}/GoSam_diphoton_ntuples/pp_aajj/${sample}/${file}_genlevelacc.root 1 1 GoSam
qsub -q all.q ~/shape_studies/theory_marco/run_outphotonproducer.sh ${T3DCAPHOME}/GoSam_diphoton_ntuples/pp_aajj/${sample}/${file}_genlevelacc.root ${PWD}/outphoton_effunf_${sample}_${file}.root

done
done