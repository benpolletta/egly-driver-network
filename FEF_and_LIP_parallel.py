#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 13:56:08 2020

@author: amelie
"""

from brian2 import *

from scipy import signal

from FEF_full3 import *
from LIP_full import *

from itertools import *
from joblib import Parallel, delayed
import multiprocessing


def save_raster(name,raster_i,raster_t,path):
    raster_file=open(path+'/raster_'+name+'_i.txt','w')
    for elem in raster_i:
        raster_file.write(str(elem)+',')
    raster_file.close()
    raster_file=open(path+'/raster_'+name+'_t.txt','w')
    for elem in raster_t:
        raster_file.write(str(elem)+',')
    raster_file.close()
    return

def generate_syn(source,target,syntype,connection_pattern,g_i,taur_i,taud_i,V_i):
    eq_syn='''_post=s_i*g_i*(V_post-V_i) : amp * meter ** -2 (summed)
        ds_i/dt=-s_i/taud_i+(1-s_i)/taur_i*0.5*(1+tanh(V_pre/10/mV)) : 1
        g_i : siemens * meter**-2
        V_i : volt
        taud_i : second
        taur_i : second
    '''

    S=Synapses(source,target,model=syntype+eq_syn,method='exact')
    if connection_pattern=='':
        S.connect()
    else :
        S.connect(j=connection_pattern, skip_if_invalid=True)
    S.g_i=g_i
    S.taur_i=taur_i
    S.taud_i=taud_i
    S.V_i=V_i  
    return S

def FEF_and_LIP(simu,path):
    prefs.codegen.target = 'numpy' 
    target_time,N_simu=simu[0],simu[1]
    
    new_path=path+"/results_"+str(N_simu)
    os.mkdir(new_path)
    
    theta_phase='mixed'
    target_on=True  
#    target_time = 500*msecond
    
    start_scope()
    close('all')

    runtime=2*second
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    
    NN=1 #multiplicative factor on the number of neurons
    N_RS,N_FS,N_SI,N_IB= NN*80,NN*20,NN*20,NN*20 #Number of neurons of RE, TC, and HTC type
    N_SI,N_RS_gran,N_SI_gran=20,20,20
    N_RS_vis,N_FS_vis,N_RS_mot,N_SI_mot,N_dSI_vm,N_RS_vm,N_gSI_vm=[20]*7

    
    all_SIdFSg=[2*msiemens * cm **-2] #1
    all_FSgRSg=[1* msiemens * cm **-2]
    all_RSgFSg=[1*msiemens * cm **-2]
    all_RSgRSg=[0.5*msiemens * cm **-2]
    all_FSgFSg=[0.3* msiemens * cm **-2]
    all_RSgRSs=[2*msiemens * cm **-2]
    all_RSgFSs=[0.1*msiemens * cm **-2]
    all_FSgRSs=[0.1* msiemens * cm **-2]
    all_J_RSg=['15 * uA * cmeter ** -2']
    all_J_FSg=['-5 * uA * cmeter ** -2']
    all_thal=[10* msiemens * cm **-2]
    thal=all_thal[0]
    
    all_syn_cond=list(product(all_SIdFSg,all_FSgRSg,all_RSgFSg,all_RSgRSg,all_FSgFSg,all_RSgRSs,all_RSgFSs,all_FSgRSs))
    all_J=list(product(all_J_RSg,all_J_FSg))
    syn_cond=all_syn_cond[0]
    J=all_J[0]
    
    slot_duration = 100*ms
    timeslots=zeros((int(around(runtime/slot_duration)),))
    target_index = int(around(target_time/slot_duration))
    timeslots[target_index]=1
    sinp_SI=TimedArray(array(timeslots), dt=slot_duration)
    
    if theta_phase=='bad':
        input_beta2_IB=False
        input_beta2_RS=False
        input_beta2_FS_SI=True
        input_thalamus_gran=True
        gFS=0* msiemens * cm **-2
        ginp_SI=0* msiemens * cm **-2
        ginpSIdeep=0* msiemens * cm **-2
        thal_cond=2* msiemens * cm **-2
        kainate='low'
        
    if theta_phase=='good':
#        input_beta2_IB=True
        input_beta2_IB=False
        ginp_IB=500* msiemens * cm **-2
        ginpSIdeep=500* msiemens * cm **-2
        input_beta2_RS=False
        input_beta2_FS_SI=False
        input_thalamus_gran=True
        thal_cond=thal
        kainate='low'
        
    if theta_phase=='mixed':
        input_mixed=True
        ginp_IB=500* msiemens * cm **-2
        ginpSIdeep=500* msiemens * cm **-2
        input_beta2_IB=False
        input_beta2_RS=False
        input_beta2_RS=False
        input_beta2_FS_SI=False
        input_thalamus_gran=False
        kainate='low'
    
    print('Network setup')
    
    net=Network()
    
    all_neurons_FEF,all_synapses_FEF,all_monitors_FEF=create_FEF_full2(N_RS_vis,N_FS_vis,N_RS_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime,target_time)
    R8,R9,R10,V_RS,V_FS,V_SI,R11,R12,R13,R14,mon_RS=all_monitors_FEF
    RSvm_FEF,SIvm_FEF,RSv_FEF,SIv_FEF,VIPv_FEF=all_neurons_FEF[1],all_neurons_FEF[2],all_neurons_FEF[6],all_neurons_FEF[9],all_neurons_FEF[8]
    
    all_neurons_LIP, all_synapses_LIP, all_gap_junctions_LIP, all_monitors_LIP=make_full_network(syn_cond,J,thal,theta_phase)
    V1,V2,V3,R1,R2,R3,I1,I2,I3,V4,R4,I4s,I4a,I4ad,I4bd,R5,R6,R7,V5,V6,V7,inpmon,inpIBmon=all_monitors_LIP
    RS_sup_LIP,SI_sup_LIP,IB_LIP,SI_deep_LIP=all_neurons_LIP[0],all_neurons_LIP[2],all_neurons_LIP[5],all_neurons_LIP[9]
    RS_gran_LIP,FS_gran_LIP=all_neurons_LIP[7],all_neurons_LIP[8]
    
    IB_LIP.ginp_IB=0* msiemens * cm **-2 #the input to RS_sup_LIP is provided with synapses from FEF 
    SI_sup_LIP.ginp_SI=10* msiemens * cm **-2
    SI_deep_LIP.ginp_SI=0* msiemens * cm **-2
#    RSvm_FEF.ginp_RS=0* msiemens * cm **-2
    SIvm_FEF.ginp_SI=0* msiemens * cm **-2
    RSv_FEF.ginp_RS=0* msiemens * cm **-2
    SIv_FEF.ginp_SI=0* msiemens * cm **-2
    VIPv_FEF.ginp_VIP_good=0* msiemens * cm **-2
    VIPv_FEF.ginp_VIP_bad=0* msiemens * cm **-2
    
    if theta_phase=='good':
        VIP_FEF=all_neurons_FEF[0]
        VIP_FEF.ginp_VIP_good=10* msiemens * cm **-2
        RS_gran_LIP.ginp_RS_good=15* msiemens * cm **-2
        FS_gran_LIP.ginp_FS_good=15* msiemens * cm **-2
        VIP_FEF.ginp_VIP_bad=10* msiemens * cm **-2
        RS_gran_LIP.ginp_RS_bad=15* msiemens * cm **-2
        FS_gran_LIP.ginp_FS_bad=15* msiemens * cm **-2
    if theta_phase=='mixed':
        VIP_FEF=all_neurons_FEF[0]
        VIP_FEF.ginp_VIP_good=10* msiemens * cm **-2
        RS_gran_LIP.ginp_RS_good=5* msiemens * cm **-2
        FS_gran_LIP.ginp_FS_good=5* msiemens * cm **-2
        RS_gran_LIP.ginp_RS_bad=5* msiemens * cm **-2
        FS_gran_LIP.ginp_FS_bad=5* msiemens * cm **-2
        VIP_FEF.ginp_VIP_bad=10* msiemens * cm **-2

    
    net.add(all_neurons_FEF)
    net.add(all_synapses_FEF)
    net.add(all_monitors_FEF)    
    
    net.add(all_neurons_LIP)
    net.add(all_synapses_LIP)
    net.add(all_gap_junctions_LIP)
    net.add(all_monitors_LIP)
    
    S_FEF_IB_LIP=generate_syn(RSvm_FEF,IB_LIP,'Isyn_FEF','',0*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_FEF_SIdeep_LIP=generate_syn(RSvm_FEF,SI_deep_LIP,'Isyn_FEF','',0.05*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_LIP_RS_FEF=generate_syn(RS_sup_LIP,RSvm_FEF,'Isyn_LIP','',0.009*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    S_LIP_FS_FEF=generate_syn(RS_sup_LIP,SIvm_FEF,'Isyn_LIP','',0.009*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   

    S_LIP_RSv_FEF=generate_syn(RS_sup_LIP,RSv_FEF,'Isyn_LIP','',0.025*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    S_LIP_SIv_FEF=generate_syn(RS_sup_LIP,SIv_FEF,'Isyn_LIP','',0.025*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    S_LIP_VIPv_FEF=generate_syn(RS_sup_LIP,VIPv_FEF,'Isyn_LIP','',0.005*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
   
    RSv_FEF.ginp_RS2=2.5* msiemens * cm **-2
    SIv_FEF.ginp_SI2=2.5* msiemens * cm **-2
    VIPv_FEF.ginp_VIP2=2.5* msiemens * cm **-2
        
    net.add(S_FEF_IB_LIP)
    net.add(S_FEF_SIdeep_LIP)
    net.add(S_LIP_RS_FEF)
    net.add(S_LIP_FS_FEF)
    net.add(S_LIP_RSv_FEF)
    net.add(S_LIP_SIv_FEF)
    net.add(S_LIP_VIPv_FEF)
    
    print('Compiling with cython')
    
    prefs.codegen.target = 'cython' #cython=faster, numpy = default python
    
#    defaultclock.dt = 0.01*ms
    net.run(runtime,report='text',report_period=300*second)
    
    save_raster('LIP RS',R1.i,R1.t,new_path)
    save_raster('LIP FS',R2.i,R2.t,new_path)
    save_raster('LIP SI',R3.i,R3.t,new_path)
    save_raster('LIP IB',R4.i,R4.t,new_path)
    save_raster('LIP RS gran',R5.i,R5.t,new_path)
    save_raster('LIP FS gran',R6.i,R6.t,new_path)
    save_raster('LIP SI deep',R7.i,R7.t,new_path)
    save_raster('FEF RS vm',R8.i,R8.t,new_path)
    save_raster('FEF SI2 vm',R9.i,R9.t,new_path)
    save_raster('FEF SI1 vm',R10.i,R10.t,new_path)
    save_raster('FEF RS v',R11.i,R11.t,new_path)
    save_raster('FEF FS v',R12.i,R12.t,new_path)
    save_raster('FEF VIP v',R13.i,R13.t,new_path)
    save_raster('FEF SI v',R14.i,R14.t,new_path)
    save_raster('FEF RS m',mon_RS.i,mon_RS.t,new_path)

    
 #   clear_cache('cython')


if __name__=='__main__':
    path=""
    if os.name == 'nt':
        path=os.path.join(ntpath.dirname(os.path.abspath(__file__)),"results_"+str(datetime.datetime.now()).replace(':','-'))
    else :
        path="./results_"+str(datetime.datetime.now())
    
    os.mkdir(path)
    
    N=1#50
    liste_target_time=[350*msecond]#[350*msecond,450*msecond,550*msecond,650*msecond,750*msecond,850*msecond,950*msecond,1050*msecond,1150*msecond,1250*msecond,1350*msecond,1450*msecond,1550*msecond,1650*msecond]
 
    liste_simus=[]
    for t in liste_target_time:
        liste_simus+=[t]*N
    
    liste_simus=[[liste_simus[i],i+750] for i in range(len(liste_simus))]
    
    simus_pas_faites=list(range(700))
    
    order=[N*i for i in range(len(liste_target_time))]
    for ind in range(1,N):
        order+=[N*i+ind for i in range(len(liste_target_time))]
    
    liste_simus=[liste_simus[i] for i in order if i in simus_pas_faites]

    liste_simus.reverse()
    
    print(liste_simus)
    
    print('Number of simulations: '+str(len(liste_simus)))
    
    for simu in liste_simus:
        FEF_and_LIP(simu,path)
        clear_cache('cython')
    
#    clear_cache('cython')