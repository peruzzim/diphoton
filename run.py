#!/usr/bin/python

from subprocess import Popen
from time import sleep
import sys
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('what', help='data/mc/data_step2')
args = parser.parse_args()

##################################################################

# limit statistics in templates
limitstat = True

if (args.what=='data'):
    file='./gg_minitree_data_141106_step1_isogeomcorr/data_step1_141106_isogeomcorr.root'
    modes=['sigsig','sigbkg','bkgbkg','standard','randomcone','sieiesideband','zee']
    activate12events=False
    isdata=1
    number=-1
elif (args.what=='mc'):
    file='./gg_minitree_mc_141106_step1_isogeomcorr/allmc.root'
    modes=['signal','randomcone','background','sieiesideband','sigsig','sigbkg','bkgbkg','standard','2pgen','2fgen','1p1fbothgen','1pgen1fside','1prcone1fgen','fragmentation','nofragmentation']
    activate12events=False
    isdata=0
    number=-1
elif (args.what=='data_step2'):
    file='./gg_minitree_data_141106_step2_isogeomcorr/data_step2_141106_isogeomcorr.root'
    modes=['standard_newtemplates_sigsig','standard_newtemplates_sigbkg','standard_newtemplates_bkgbkg']
    activate12events=True
    isdata=1
    number=-1
else:
    print 'Wrong usage, exiting.'
    sys.exit()

##################################################################

spawn=8
lista_processi=[]


def wait_processes():
    sleep(1)
    for i in lista_processi:
        if (i.poll()!=None):
            lista_processi.remove(i)
    if (len(lista_processi)>=spawn):
        wait_processes()
    else:
        return

if (isdata==0):
    strdata='allmc'
else:
    strdata='data'

for mode in modes:
    wait_processes()
    thisnumber = number
    if (mode.rfind('standard')>=0):
        thisnumber = -1
    if limitstat:
        if (mode.rfind('randomcone')>=0 or mode.rfind('signal')>=0 or mode.rfind('background')>=0 or mode.rfind('nofragmentation')>=0 or mode.rfind('sieiesideband')>=0 or mode.rfind('sigsig')>=0):
            thisnumber=1e5
    if (activate12events):
        args = ['root','-q','-b','-l','template_production.C+O("'+file+'","'+mode+'",'+str(isdata)+',"outphoton/outphoton_'+strdata+'_'+mode+'_1event.root","photoniso",'+str(thisnumber)+',false);']
        lista_processi.append(Popen(args))
        if (mode=='standard_newtemplates_bkgbkg'):
            wait_processes()
            args = ['root','-q','-b','-l','template_production.C+O("'+file+'","'+mode+'",'+str(isdata)+',"outphoton/outphoton_'+strdata+'_'+mode+'_2events.root","photoniso",'+str(thisnumber)+',true);']
            lista_processi.append(Popen(args))
    else:
        args = ['root','-q','-b','-l','template_production.C+O("'+file+'","'+mode+'",'+str(isdata)+',"outphoton/outphoton_'+strdata+'_'+mode+'.root","photoniso",'+str(thisnumber)+');']
        lista_processi.append(Popen(args))

for i in lista_processi:
    i.wait()




