#!/usr/bin/python

##################################################################

modes=[]
isdata=[]

# DATA (do not swap with mc)
modes.append(['ESCALEup','ESCALEdown','JECup','JECdown'])
isdata.append(1)

# MC (do not swap with data)
modes.append(['ESMEARup','ESMEARdown','JERup','JERdown','PUup','PUdown'])
isdata.append(0)


#################################################################

number=-1
spawn=4
lista_processi=[]

from subprocess import Popen
from time import sleep

def wait_processes():
    sleep(1)
    for i in lista_processi:
        if (i.poll()!=None):
            lista_processi.remove(i)
    if (len(lista_processi)>=spawn):
        wait_processes()
    else:
        return


for i in xrange(2):
    if (isdata[i]==0):
        strdata='sig'
    else:
        strdata='data'
    for mode in modes[i]:
        wait_processes()
        file1='outphoton/outphoton_effunf_'+strdata+'_'+mode+'.root'
        file2='outphoton/outphoton_effunf_'+strdata+'_'+'Default'+'.root'
        args = ['root','-q','-b','-l','make_histos_syst_rawyield.C("'+file1+'","'+file2+'","'+mode+'");"']
        lista_processi.append(Popen(args))
            
for i in lista_processi:
    i.wait()


