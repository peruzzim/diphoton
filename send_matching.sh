for i in 41 43 44 50; do qsub -q all.q -l hostname=t3wn$i fetch_matching.sh; done
