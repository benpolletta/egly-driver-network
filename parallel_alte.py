#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 13:56:56 2020

@author: amelie
"""

print('Starting')

import os

os.system('cd ~/models/model_FEF_LIP')

import matplotlib
import matplotlib.pyplot as plt
import time
#matplotlib.use('agg')
plt.switch_backend('agg')
import sys
sys.path.insert(0, '~/models/model_FEF_LIP')

from brian2 import *

from joblib import Parallel, delayed
import multiprocessing
import os

from FEF_and_LIP_parallel_alte import *

os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_DYNAMIC'] = 'FALSE'

import time
import ntpath
from itertools import *


path=""
if os.name == 'nt':
    path=os.path.join(ntpath.dirname(os.path.abspath(__file__)),"results_"+str(datetime.datetime.now()).replace(':','-'))
else :
#    path="/projectnb/crc-nak/aaussel/results_alte_"+str(datetime.datetime.now())
    #path="/project/crc-nak/aaussel/results_alte_"+str(datetime.datetime.now())
    path="/results_alte_"+str(datetime.datetime.now())

os.mkdir(path)

#Network parameters
#N=50
#liste_target_time=[300*msecond,400*msecond,500*msecond,600*msecond,700*msecond,800*msecond,900*msecond,1000*msecond,1100*msecond,1200*msecond,1300*msecond,1400*msecond,1500*msecond,1600*msecond,1700*msecond]
#liste_target_time+=[350*msecond,450*msecond,550*msecond,650*msecond,750*msecond,850*msecond,950*msecond,1050*msecond,1150*msecond,1250*msecond,1350*msecond,1450*msecond,1550*msecond,1650*msecond]
#
#liste_target_time+=[325*msecond,425*msecond,525*msecond,625*msecond,725*msecond,825*msecond,925*msecond,1025*msecond,1125*msecond,1225*msecond,1325*msecond,1425*msecond,1525*msecond,1625*msecond,1725*msecond]
#liste_target_time+=[375*msecond,475*msecond,575*msecond,675*msecond,775*msecond,875*msecond,975*msecond,1075*msecond,1175*msecond,1275*msecond,1375*msecond,1475*msecond,1575*msecond,1675*msecond]

N=1
liste_target_time=[500*msecond]
                   
Atheta=1
Asqrtheta=0
Aalpha=5/2* msiemens * cm **-2

liste_simus=[]
for t in liste_target_time:
    liste_simus+=[t]*N

liste_simus=[[liste_simus[i],Atheta,Asqrtheta,Aalpha,i] for i in range(len(liste_simus))]

simus_pas_faites=list(range(len(liste_simus)))

order=[50*i for i in range(len(liste_target_time))]
for ind in range(1,50):
    order+=[50*i+ind for i in range(len(liste_target_time))]

liste_simus=[liste_simus[i] for i in order if i in simus_pas_faites]

print('Number of simulations: '+str(len(liste_simus)))

#setting the number of cores to used (all cpus by default)
#num_cores = multiprocessing.cpu_count()
#if os.name == 'nt':
#    num_cores=-3 #using all cpus on a windows does not work for an unknown reason

for simu in liste_simus:
    FEF_and_LIP(simu,path)

clear_cache('cython')



