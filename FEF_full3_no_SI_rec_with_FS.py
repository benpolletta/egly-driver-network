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
from cells.SI_FEF_more_h import *
from cells.VIP_FEF import *

from FEF_visuomotor_no_SI_rec_with_FS import *
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



def create_FEF_full2(N_RS_vis,N_FS_vis,N_RS_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime,target_time):
    
    #create each functional group of neurons individually
    all_neurons_vm,all_synapses_vm,all_monitors_vm=generate_deepSI_and_gran_layers(theta_phase,N_dSI_vm,N_RS_vm,N_gSI_vm,runtime)
    RS_vm=all_neurons_vm[1]
    
    all_neurons_v,all_synapses_v,all_monitors_v=generate_visual_neurons(theta_phase,N_FS_vis,N_RS_vis,runtime,target_on,target_time)
    RS_vis=all_neurons_v[0]
    VIP,SI,Poisson_target=all_neurons_v[2],all_neurons_v[3],all_neurons_v[5]
#    print(all_neurons_v)
#    print(VIP.ginp_VIP2,SI.ginp_SI2)
#    print(SI.ginp_SI,VIP.ginp_VIP_good,VIP.ginp_VIP_bad)
#    print(Poisson_target.spike_time)
#    print(VIP.Iapp,SI.J)
#    print()
#    print(all_synapses_v)
    
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
    taurinp2=2*ms
    taudinp2=10*ms
    tauinp2=taudinp2  
#    ginp_SI=0* msiemens * cm **-2
    
    print('Creating the network')
    N_RS_vis,N_FS_vis,N_RS_mot,N_dSI_vm,N_RS_vm,N_gSI_vm=[20]*6
    
#    theta_phase='mixed'
#    theta_phase='good'
    theta_phase='bad'
    target_on=True
    runtime=1*second
    target_time=500*msecond
#    target_time=1500*msecond
    
    net=Network()
    
#    net,all_monitors=create_network(N_RS_vis,N_FS_vis,N_RS_mot,N_SI_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime)
    all_neurons,all_synapses,all_monitors=create_FEF_full2(N_RS_vis,N_FS_vis,N_RS_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime,target_time)
    
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
    R1,R2,R3,R_FS,V1,V2,V3,R4,R5,R6,R7,mon_RS=all_monitors
    
    figure(figsize=(10,4))
    subplot(131)
    title('Visual Neurons')
    plot(R4.t,R4.i+20,'r.',label='RS')
    plot(R5.t,R5.i+0,'k.',label='FS')
    plot(R6.t,R6.i+40,'b.',label='VIP')
    plot(R7.t,R7.i+60,'g.',label='SI')
    xlim(0,runtime/second)
    legend(loc='upper left')   
    
    subplot(132)
    title('Visual-Motor Neurons')
    plot(R3.t,R3.i+0,'c.',label='VIP')
    plot(R1.t,R1.i+60,'r.',label='RS')
    plot(R2.t,R2.i+40,'b.',label='SOM')
    plot(R_FS.t,R_FS.i+80,'g.',label='FS')
    xlim(0,runtime/second)
    legend(loc='upper left') 
    
    subplot(133)
    title('Motor Neurons')
    plot(mon_RS.t,mon_RS.i+0,'r.',label='RS')
    xlim(0,runtime/second)
    legend(loc='upper left') 
    
    
    figure(figsize=(9,4))
    subplot(131)
    title('Visual Neurons')
    plot(R4.t,R4.i+0,'r.',label='RS')
    plot(R5.t,R5.i+20,'b.',label='FS')
    plot(R6.t,R6.i+40,'.',label='VIP',color='lime')
    plot(R7.t,R7.i+60,'g.',label='SOM')
    xlim(0.2,runtime/second)
    legend(loc='upper left') 
    xlabel('Time (s)')
    ylabel('Neuron index')
    
    subplot(132)
    title('Visual-Motor Neurons')
    plot(R3.t,R3.i+0,'g.',label='VIP')
    plot(R1.t,R1.i+60,'r.',label='RS')
    plot(R2.t,R2.i+40,'.',label='SOM',color='lime')
    plot(R_FS.t,R_FS.i+80,'g.',label='FS')
    xlim(0.2,runtime/second)
    legend(loc='upper left') 
    xlabel('Time (s)')
    ylabel('Neuron index')
    
    subplot(133)
    title('Motor Neurons')
    plot(mon_RS.t,mon_RS.i+0,'r.',label='RS')
    xlim(0.2,runtime/second)
    legend(loc='upper left') 
    xlabel('Time (s)')
    ylabel('Neuron index')
    
    tight_layout()
    
    
#    min_t=int(50*ms*100000*Hz)
#    LFP_V1=1/20*sum(V1.V,axis=0)[min_t:]
#    LFP_V2=1/20*sum(V2.V,axis=0)[min_t:]
#    LFP_V3=1/20*sum(V3.V,axis=0)[min_t:]
#    LFP_V4=1/20*sum(V4.V,axis=0)[min_t:]
#    LFP_V5=1/20*sum(V5.V,axis=0)[min_t:]
##    LFP_V6=1/20*sum(V6.V,axis=0)[min_t:]
##    LFP_V7=1/20*sum(V7.V,axis=0)[min_t:]
#    
#    f,Spectrum_LFP_V1=signal.periodogram(LFP_V1, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V2=signal.periodogram(LFP_V2, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V3=signal.periodogram(LFP_V3, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V4=signal.periodogram(LFP_V4, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V5=signal.periodogram(LFP_V5, 100000,'flattop', scaling='spectrum')
##    f,Spectrum_LFP_V6=signal.periodogram(LFP_V6, 100000,'flattop', scaling='spectrum')
##    f,Spectrum_LFP_V7=signal.periodogram(LFP_V7, 100000,'flattop', scaling='spectrum')
#
#    figure(figsize=(10,4))
#    subplot(331)
#    plot(f,Spectrum_LFP_V4)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('visual RS')
#    subplot(334)
#    plot(f,Spectrum_LFP_V5)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('visual FS')  
#    
#    subplot(332)
#    plot(f,Spectrum_LFP_V1)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('visual-motor gran RS')
#    subplot(335)
#    plot(f,Spectrum_LFP_V2)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('visual-motor gran SI')  
#    subplot(338)
#    plot(f,Spectrum_LFP_V3)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('visual-motor deep SI')  
    
#    subplot(333)
#    plot(f,Spectrum_LFP_V6)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('motor RS')
#    subplot(336)
#    plot(f,Spectrum_LFP_V7)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('motor SI')
#    tight_layout()
    
#    figure()
#    plot(mon_RS.t,mon_RS.Isyn[0])
#    plot(mon_RS.t,mon_RS.Isyn[10])
    
    clear_cache('cython')
    