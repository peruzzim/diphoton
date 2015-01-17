#!/usr/bin/env python
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', help='Force recreation of allbins files', action='store_true')
args = parser.parse_args()

lista_cats=['EBEB','EBEE','EEEE']

lista_vars=[
"invmass",
"diphotonpt",
"costhetastar",
"dphi",
"dR",
"njets",
"dR_gg_closestjet",
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


forcestring=''
if (args.f):
	forcestring=' -f '

os.system('mkdir histos_purity_eachbin')
for var in lista_vars:
	for cat in lista_cats:
		mycomm='mv histo_purity_'+var+'_'+cat+'_b*.root histos_purity_eachbin'
		os.system('hadd '+forcestring+' histo_purity_'+var+'_'+cat+'_allbins.root'+' '+'histo_purity_'+var+'_'+cat+'_b*.root && '+mycomm)

