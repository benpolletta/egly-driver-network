#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 13:56:56 2020

@author: amelie
"""

print('Starting')

import os

os.system('cd /projectnb/crc-nak/brpp/Full_model_04_06')

import matplotlib
import matplotlib.pyplot as plt
import time
#matplotlib.use('agg')
plt.switch_backend('agg')
import sys
sys.path.insert(0, '/projectnb/crc-nak/brpp/Full_model_04_06')

from brian2 import *

from joblib import Parallel, delayed
import multiprocessing
import os
import sys

from FEF_and_LIP_parallel import *

if sys.platform=='linux':
    cache_dir=os.environ['TMPDIR']
    prefs.codegen.runtime.cython.cache_dir = cache_dir
    prefs.codegen.runtime.cython.multiprocess_safe = False

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
    path="/project/crc-nak/brpp/results_"+str(datetime.datetime.now()).replace(':','-')

os.mkdir(path)

#Basic network parameters
N=50
liste_target_time=[300*msecond,400*msecond,500*msecond,600*msecond,700*msecond,800*msecond,900*msecond,1000*msecond,1100*msecond,1200*msecond,1300*msecond,1400*msecond,1500*msecond,1600*msecond,1700*msecond]

liste_simus=[]
for t in liste_target_time:
    liste_simus+=[t]*N

liste_simus=[[liste_simus[i],i] for i in range(len(liste_simus))]

print('Number of simulations: '+str(len(liste_simus)))

#setting the number of cores to used (all cpus by default)
num_cores = multiprocessing.cpu_count()
if os.name == 'nt':
    num_cores=-3 #using all cpus on a windows does not work for an unknown reason

Parallel(n_jobs=num_cores)(delayed(FEF_and_LIP)(simu,path) for simu in liste_simus)

clear_cache('cython')