#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 11:24:49 2020

@author: amelie
"""

from brian2 import *

from scipy import signal
from cells.RS_FEF import *
from cells.FS_FEF import *
from cells.SI_FEF import *
from cells.VIP_FEF import *

from FEF_visuomotor_tests import *
from FEF_contrast_detection_3 import *

runtime=1*second
    
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
        S.connect(condition=connection_pattern, skip_if_invalid=True)
    S.g_i=g_i
    S.taur_i=taur_i
    S.taud_i=taud_i
    S.V_i=V_i  
    return S


def generate_spike_timing(N,f,start_time,end_time=runtime):
    list_time_and_i=[]
    for i in range(N):
        list_time=[(start_time,i)]
        next_spike=list_time[-1][0]+(1+0.01*rand())/f
        while next_spike<end_time:
            list_time.append((next_spike,i))
            next_spike=list_time[-1][0]+(1+0.01*rand())/f
        list_time_and_i+=list_time
    return array(list_time_and_i)



def create_FEF_full2(N_RS_vis,N_FS_vis,N_RS_mot,N_VIP_vm,N_RS_vm,N_FS_vm,N_gSI_vm,theta_phase,target_on,runtime,target_time):
    
    #create each functional group of neurons individually
    all_neurons_vm,all_synapses_vm,all_monitors_vm=generate_deepSI_and_gran_layers(theta_phase,N_VIP_vm,N_RS_vm,N_FS_vm,N_gSI_vm,runtime)
    RS_vm=all_neurons_vm[1]
    
    all_neurons_v,all_synapses_v,all_monitors_v=generate_visual_neurons(theta_phase,N_FS_vis,N_RS_vis,runtime,target_on,target_time)
    RS_vis=all_neurons_v[0]
#    print(all_neurons_v)
    
    RS_mot=NeuronGroup(N_RS_mot,eq_RS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    RS_mot.V = '-70*mvolt+10*rand()*mvolt'
    RS_mot.h = '0+0.05*rand()'
    RS_mot.m = '0+0.05*rand()'
    RS_mot.mAR = '0.035+0.025*rand()'
#    RS_mot.J='50 * uA * cmeter ** -2'  
    RS_mot.J='50 * uA * cmeter ** -2'  
    
    #From visual to visual-motor
##    S_RSvRSm=generate_syn(RS_vis,RS_mot,'IsynRS_FEF_V','',0.15*msiemens * cm **-2,12.5*ms,125*ms,0*mV)
##    S_RSvRSm_AMPA=generate_syn(RS_vis,RS_mot,'IsynRS_FEF_V','',0.04*msiemens * cm **-2,0.125*ms,1*ms,0*mV) #AMPA 0.125*ms,1*ms NMDA 12.5*ms,125*ms
#    S_RSvRSm_AMPA=generate_syn(RS_vis,RS_mot,'IsynRS_FEF_V','i<10',0.08*msiemens * cm **-2,0.125*ms,1*ms,0*mV) #AMPA 0.125*ms,1*ms NMDA 12.5*ms,125*ms
##    S_RSvRSm_AMPA=generate_syn(RS_vis,RS_vm,'IsynRS_FEF_V','i<10',0.1*msiemens * cm **-2,0.125*ms,1*ms,0*mV) #AMPA 0.125*ms,1*ms NMDA 12.5*ms,125*ms
##0.06
##    S_RSvmRSm_AMPA=generate_syn(RS_vm,RS_mot,'IsynFS_FEF_V','',0*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
#    S_RSvmRSm_AMPA=generate_syn(RS_vm,RS_mot,'IsynFS_FEF_V','',0.01*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
#    S_RSvRSm_NMDA=generate_syn(RS_vis,RS_mot,'IsynSI_FEF_VM','',0.01*msiemens * cm **-2,12.5*ms,125*ms,0*mV) #AMPA 0.125*ms,1*ms NMDA 12.5*ms,125*ms
##    S_RSvRSm_NMDA=generate_syn(RS_vis,RS_vm,'IsynFS_FEF_V','',0.05*msiemens * cm **-2,12.5*ms,125*ms,0*mV) #AMPA 0.125*ms,1*ms NMDA 12.5*ms,125*ms
##    S_RSvmRSm_NMDA=generate_syn(RS_vm,RS_mot,'IsynSI2_FEF_VM','',0*msiemens * cm **-2,12.5*ms,125*ms,0*mV)
#    S_RSvmRSm_NMDA=generate_syn(RS_vm,RS_mot,'IsynSI2_FEF_VM','',0.1*msiemens * cm **-2,12.5*ms,125*ms,0*mV)


    S_RSvRSm_AMPA=generate_syn(RS_vis,RS_mot,'IsynRS_FEF_V','i<10',0.08*msiemens * cm **-2,0.125*ms,1*ms,0*mV) #AMPA 0.125*ms,1*ms NMDA 12.5*ms,125*ms
    S_RSvmRSm_AMPA=generate_syn(RS_vm,RS_mot,'IsynFS_FEF_V','',0*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_RSvRSm_NMDA=generate_syn(RS_vis,RS_mot,'IsynSI_FEF_VM','',0*msiemens * cm **-2,12.5*ms,125*ms,0*mV) #AMPA 0.125*ms,1*ms NMDA 12.5*ms,125*ms
    S_RSvmRSm_NMDA=generate_syn(RS_vm,RS_mot,'IsynSI2_FEF_VM','',0.08*msiemens * cm **-2,12.5*ms,125*ms,0*mV)
  
    eq_gap='''_post=g_i*(V_post-V_pre) : amp * meter ** -2 (summed)
        g_i : siemens * meter**-2
    '''
    
    gapRSmot=Synapses(RS_mot,RS_mot,model='Igap'+eq_gap,method='exact')
    gapRSmot.connect(j='i')
    gapRSmot.g_i=0* msiemens * cm **-2

       
    mon_RS=SpikeMonitor(RS_mot,record=True)
        
    all_monitors=all_monitors_vm+all_monitors_v+(mon_RS,)
    all_neurons=all_neurons_vm+all_neurons_v+(RS_mot,)
    all_synapses=all_synapses_vm+all_synapses_v+(S_RSvRSm_AMPA,S_RSvmRSm_AMPA,S_RSvRSm_NMDA,S_RSvmRSm_NMDA,gapRSmot)
    
    
    return all_neurons,all_synapses,all_monitors


if __name__=='__main__':
    close('all')
    prefs.codegen.target = 'numpy'
    defaultclock.dt = 0.01*ms
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp=0* msiemens * cm **-2
    ginp_SI=0* msiemens * cm **-2
    
    print('Creating the network')
    N_RS_vis,N_FS_vis,N_RS_mot,N_VIP_vm,N_RS_vm,N_FS_vm,N_gSI_vm=[20]*7
    
    theta_phase='mixed'
    target_on=True
    runtime=2*second
    target_time=650*msecond
#    target_time=1500*msecond
    
    net=Network()
    
#    net,all_monitors=create_network(N_RS_vis,N_FS_vis,N_RS_mot,N_SI_mot,N_VIP_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime)
    all_neurons,all_synapses,all_monitors=create_FEF_full2(N_RS_vis,N_FS_vis,N_RS_mot,N_VIP_vm,N_RS_vm,N_FS_vm,N_gSI_vm,theta_phase,target_on,runtime,target_time)
    
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_monitors)
    
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp       
    taurinp2=2*ms
    taudinp2=10*ms
    tauinp2=taudinp2   
    
    print('Compiling with cython')
    prefs.codegen.target = 'cython'
    net.run(runtime,report='text',report_period=300*second)
    
#    R1,R2,R3,V1,V2,V3,R4,R5,V4,V5,R6,R7,V6,V7,R8,mon_RS=all_monitors
    RSvm,SIvm,VIPvm,FSvm,V_RS,V_SI,V_VIP,BUinp,TDinp,RSv,FSv,VIPv,SIv,RSm=all_monitors
    
    figure(figsize=(9,4))
    subplot(131)
    title('Visual Neurons')
    plot(RSv.t,RSv.i+60,'r.',label='RS')
    plot(FSv.t,FSv.i+40,'b.',label='FS')
    plot(SIv.t,SIv.i+20,'g.',label='SI')
    plot(VIPv.t,VIPv.i+0,'k.',label='VIP')
    xlim(0,runtime/second)
    legend(loc='upper left')   
    
    subplot(132)
    title('Visual-Motor Neurons')
    plot(RSvm.t,RSvm.i+60,'r.',label='RS')
    plot(FSvm.t,FSvm.i+40,'b.',label='FS')
    plot(SIvm.t,SIvm.i+20,'g.',label='SI')
    plot(VIPvm.t,VIPvm.i+0,'k.',label='VIP')
    xlim(0,runtime/second)
    legend(loc='upper left') 
    
    subplot(133)
    title('Motor Neurons')
    plot(RSm.t,RSm.i+0,'r.',label='RS')
    xlim(0,runtime/second)
    legend(loc='upper left')
    
    tight_layout()
    
    savefig('sims/FEF_full3')
    
    clear_cache('cython')
    