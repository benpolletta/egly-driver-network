

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 13:56:56 2020

@author: amelie
"""

#This is used for simulations that change the time constants of all interneurons of one type at the same time.

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

from FEF_and_LIP_no_SI_rec_parallel_alternate_v2 import *

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
    path="/projectnb2/crc-nak/aaussel/results_tSIv2_"+str(datetime.datetime.now())

os.mkdir(path)

N=20
liste_target_time=[300*msecond,400*msecond,500*msecond,600*msecond,700*msecond,800*msecond,900*msecond,1000*msecond,1100*msecond,1200*msecond,1300*msecond,1400*msecond,1500*msecond,1600*msecond,1700*msecond]

liste_target_time+=[350*msecond,450*msecond,550*msecond,650*msecond,750*msecond,850*msecond,950*msecond,1050*msecond,1150*msecond,1250*msecond,1350*msecond,1450*msecond,1550*msecond,1650*msecond]

liste_simus=[]
for t in liste_target_time:
    liste_simus+=[t]*N

#liste_t_SI=[5*ms,10*ms,15*ms,20*ms,25*ms]
#liste_t_FS=[5*ms]
liste_t_SI=[30*ms]
liste_t_FS=[5*ms]
liste_simus=[[i,j,k] for i in liste_simus for j in liste_t_SI for k in liste_t_FS]
liste_simus=[[liste_simus[i][0],i,liste_simus[i][1],liste_simus[i][2]] for i in range(len(liste_simus))]

simus_pas_faites=list(range(len(liste_simus)))
#simus_pas_faites=[]

#order=[20*i for i in range(len(liste_target_time))]
#for ind in range(1,20):
#    order+=[20*i+ind for i in range(len(liste_target_time))]

#liste_simus=[liste_simus[i] for i in order if i in simus_pas_faites]

#setting the number of cores to used (all cpus by default)
num_cores = multiprocessing.cpu_count()
if os.name == 'nt':
    num_cores=-3 #using all cpus on a windows does not work for an unknown reason

Parallel(n_jobs=num_cores)(delayed(FEF_and_LIP)(simu,path) for simu in liste_simus)





path=""
if os.name == 'nt':
    path=os.path.join(ntpath.dirname(os.path.abspath(__file__)),"results_"+str(datetime.datetime.now()).replace(':','-'))
else :
    path="/projectnb2/crc-nak/aaussel/results_tFSv2_"+str(datetime.datetime.now())

os.mkdir(path)

N=20
liste_target_time=[300*msecond,400*msecond,500*msecond,600*msecond,700*msecond,800*msecond,900*msecond,1000*msecond,1100*msecond,1200*msecond,1300*msecond,1400*msecond,1500*msecond,1600*msecond,1700*msecond]

liste_target_time+=[350*msecond,450*msecond,550*msecond,650*msecond,750*msecond,850*msecond,950*msecond,1050*msecond,1150*msecond,1250*msecond,1350*msecond,1450*msecond,1550*msecond,1650*msecond]

liste_simus=[]
for t in liste_target_time:
    liste_simus+=[t]*N

#liste_t_SI=[20*ms]
#liste_t_FS=[3*ms,5*ms,10*ms,15*ms,20*ms]
liste_t_SI=[20*ms]
liste_t_FS=[25*ms,30*ms]
liste_simus=[[i,j,k] for i in liste_simus for j in liste_t_SI for k in liste_t_FS]
liste_simus=[[liste_simus[i][0],i,liste_simus[i][1],liste_simus[i][2]] for i in range(len(liste_simus))]

simus_pas_faites=list(range(len(liste_simus)))

#order=[20*i for i in range(len(liste_target_time))]
#for ind in range(1,20):
#    order+=[20*i+ind for i in range(len(liste_target_time))]

#liste_simus=[liste_simus[i] for i in order if i in simus_pas_faites]
liste_simus=[liste_simus[i] for i in simus_pas_faites]


#setting the number of cores to used (all cpus by default)
num_cores = multiprocessing.cpu_count()
if os.name == 'nt':
    num_cores=-3 #using all cpus on a windows does not work for an unknown reason

Parallel(n_jobs=num_cores)(delayed(FEF_and_LIP)(simu,path) for simu in liste_simus)

clear_cache('cython')