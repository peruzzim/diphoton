#!/usr/bin/env python
import os

lista_cats=['EBEB','EBEE','EEEE']

lista_vars=[
"invmass",
"diphotonpt",
"costhetastar",
"dphi",
"dR",
"njets",
"1jet_jpt",
"1jet_dR_lead_j",
"1jet_dR_trail_j",
"1jet_dR_close_j",
"1jet_dR_far_j",
"2jet_j1pt",
"2jet_j2pt",
"2jet_deta_jj",
"2jet_dphi_jj",
"2jet_dR_jj",
"2jet_mjj",
"2jet_zeppen",
"2jet_dphi_gg_jj"
]


lista_nbins=[16,21,8,14,22,5,9,13,14,14,13,5,3,5,5,5,5,5,5]

option='' # default
#option='newtemplates_1event'

counter=0

for cat in lista_cats:
    for i in range(len(lista_vars)):
        for j in range(lista_nbins[i]-1):
            if (option=='templateshapeMCpromptdrivenEB' or option=='templateshapeMCfakedrivenEB'):
                if (cat=='EEEE'):
                    continue
            if (option=='templateshapeMCpromptdrivenEE' or option=='templateshapeMCfakedrivenEE'):
                if (cat=='EBEB'):
                    continue
            os.system('qsub -q short.q -o /dev/null -e /dev/null /shome/peruzzi/shape_studies/run_fits.sh ' + lista_vars[i] + ' ' + cat + ' ' + str(j) + option)
            counter=counter+1

print str(counter)+' jobs submitted'


