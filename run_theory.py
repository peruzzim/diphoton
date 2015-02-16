#!/usr/bin/python

from subprocess import Popen
from time import sleep
import sys
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('what', help='sherpa/amcatnlo')
args = parser.parse_args()

##################################################################

if (args.what=='sherpa'):
    dir='./theory_marco'
    files=['sherpa_central.root','sherpa_scaledown.root','sherpa_scaleup.root']
    scales = 0
    pdfs = 0
elif (args.what=='amcatnlo'):
    dir='./theory_marco'
    files=['amcatnlo_012j_tuneCUETP8M1_looseacc_pythia8_205.root','amcatnlo_012j_tune4C_looseacc_pythia8_205.root','amcatnlo_012j_noMPI_looseacc_pythia8_205.root']
    scales = 8
    pdfs = 100
elif (args.what=='box'):
    dir='./theory_marco'
    files=['pythia8box.root']
    scales = 0
    pdfs = 0
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

for file in files:
    wait_processes()
    args = ['root','-q','-b','-l','template_production.C+O("'+dir+'/'+file+'","effunf",0,"theory_marco/outphoton_theory_'+file+'","photoniso",-1,false,"","Default",'+str(scales)+','+str(pdfs)+');']
    lista_processi.append(Popen(args))
        
for i in lista_processi:
    i.wait()




