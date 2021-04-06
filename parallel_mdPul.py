#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 09:06:37 2021

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

from FEF_and_LIP_and_mdpul_2col_parallel import *

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
    path="/project/crc-nak/aaussel/results_"+str(datetime.datetime.now())

os.mkdir(path)

#Basic network parameters
N=10
liste_target_time=[300*msecond,400*msecond,500*msecond,600*msecond,700*msecond,800*msecond,900*msecond,1000*msecond,1100*msecond,1200*msecond,1300*msecond,1400*msecond,1500*msecond,1600*msecond,1700*msecond]
liste_J_RS_inter=['100 * uA * cmeter ** -2','80 * uA * cmeter ** -2','60 * uA * cmeter ** -2','40 * uA * cmeter ** -2','20 * uA * cmeter ** -2','0 * uA * cmeter ** -2']

liste_simus_time=[]
for t in liste_target_time:
    liste_simus_time+=[t]*N
    
liste_simus=[[liste_simus_time[i], liste_J_RS_inter[0]] for i in range(len(liste_simus_time))]
for J_RS_inter in liste_J_RS_inter[1:]:
    liste_simus=liste_simus+[[liste_simus_time[i], J_RS_inter] for i in range(len(liste_simus_time))]

liste_simus=[[liste_simus[i][0],liste_simus[i][1],i] for i in range(len(liste_simus))]

print('Number of simulations: '+str(len(liste_simus)))

#setting the number of cores to used (all cpus by default)
num_cores = multiprocessing.cpu_count()
if os.name == 'nt':
    num_cores=-3 #using all cpus on a windows does not work for an unknown reason

Parallel(n_jobs=num_cores)(delayed(run_simu)(simu,path) for simu in liste_simus)

clear_cache('cython')