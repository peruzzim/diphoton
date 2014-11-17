#Set MAX to listauuid.size()/5 approx. by defect
MAX=521
for i in $(seq 0 ${MAX}); do qsub -q all.q run_prepare_matchingfile_forstep2.sh $i; done
