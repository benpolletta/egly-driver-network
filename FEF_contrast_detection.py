#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 08:58:41 2020

@author: amelie
"""

from brian2 import *

from scipy import signal
from cells.RS_FEF import *
from cells.FS_FEF import *
from cells.SI_FEF import *
from cells.VIP_FEF_vis import *

def generate_visual_neurons(theta_phase,N_FS,N_RS,runtime,target_on):
    
    if theta_phase=='bad':
        ginp_IB=0* msiemens * cm **-2
        input_beta2_RS=False
        input_beta2_FS_SI=True
        input_thalamus_gran=True
        gFS=0* msiemens * cm **-2
#        thal_cond=10* msiemens * cm **-2
        thal_cond=10* msiemens * cm **-2
        kainate='low'
        
    if theta_phase=='good' or theta_phase=='mixed':
        ginp_IB=10* msiemens * cm **-2
        input_beta2_RS=False
        input_beta2_FS_SI=False
        input_thalamus_gran=True
#        thal_cond=10* msiemens * cm **-2
        thal_cond=7.5* msiemens * cm **-2
        kainate='low'
        

    
    prefs.codegen.target = 'numpy'
    
    defaultclock.dt = 0.01*ms
    
    #Single column network
    
    ##Define neuron groups
    RS=NeuronGroup(N_RS,eq_RS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    RS.V = '-70*mvolt+10*rand()*mvolt'
    RS.h = '0+0.05*rand()'
    RS.m = '0+0.05*rand()'
    RS.mAR = '0.035+0.025*rand()'
    RS.J='20 * uA * cmeter ** -2'  #article SI=25, code=1
    #0 
#    RS.J='125 * uA * cmeter ** -2'  #article SI=25, code=1
    
#    FS=NeuronGroup(N_FS,eq_VIP_vis,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    FS=NeuronGroup(N_FS,eq_FS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    FS.V = '-70*mvolt+10*rand()*mvolt'
    FS.h = '0+0.05*rand()'
    FS.m = '0+0.05*rand()'
#    FS.J='6 * uA * cmeter ** -2' #article=code=35
    FS.J='4 * uA * cmeter ** -2' #article=code=35
#     #-30
#    FS.Iapp='0 * uA * cmeter ** -2' #article=code=35
    
    ##Synapses
    eq_syn='''_post=s_i*g_i*(V_post-V_i) : amp * meter ** -2 (summed)
        ds_i/dt=-s_i/taud_i+(1-s_i)/taur_i*0.5*(1+tanh(V_pre/10/mV)) : 1
        g_i : siemens * meter**-2
        V_i : volt
        taud_i : second
        taur_i : second
    '''
    
    def generate_syn(source,target,syntype,connection_pattern,g_i,taur_i,taud_i,V_i):
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
    
    
    #From E (granular layer) cells
    S_RSFS=generate_syn(RS,FS,'IsynRS_FEF_VM','',0.4* msiemens * cm **-2,0.125*ms,1*ms,0*mV) #0.35
    S_FSRS=generate_syn(FS,RS,'IsynSI_FEF_VM','',0.4* msiemens * cm **-2,0.25*ms,5*ms,-80*mV) #0.35

    
    def generate_spike_timing(N,f,start_time,end_time=runtime):
        list_time_and_i=[]
        for i in range(N):
            list_time=[(start_time,i)]
            next_spike=list_time[-1][0]+(1+0*rand())/f  #0.01
            while next_spike<end_time:
                list_time.append((next_spike,i))
                next_spike=list_time[-1][0]+(1+0*rand())/f
            list_time_and_i+=list_time
        return array(list_time_and_i)
    
        
    if input_thalamus_gran:    
#        FS.ginp_VIP_good=thal_cond*1
#        FS.ginp_VIP_bad=thal_cond*1
        FS.ginp_FS=thal_cond*1
        RS.ginp_RS=thal_cond*1
        RS.ginp_RS2=thal_cond*1.5
#        RS.ginp_RS2=thal_cond*0.5
        if theta_phase=='good':
            fLIP=50*Hz
        else :
#            fLIP=50*Hz
            fLIP=13*Hz
#        print(fLIP)
        gamma_background=generate_spike_timing(N_FS,fLIP,0*ms,end_time=2100*ms)
        gamma_target=generate_spike_timing(20,50*Hz,600*ms,end_time=700*ms)
        
        Poisson_background = SpikeGeneratorGroup(N_FS, gamma_background[:,1], gamma_background[:,0]*second)
        if target_on :
            Poisson_target = SpikeGeneratorGroup(20, gamma_target[:,1], gamma_target[:,0]*second)
        else :
            Poisson_target = SpikeGeneratorGroup(20, [], []*second)
        S_in_bg_FS=Synapses(Poisson_background,FS,on_pre='Vinp=Vhigh')
        S_in_bg_FS.connect(j='i')
        S_in_bg_RS=Synapses(Poisson_background,RS,on_pre='Vinp=Vhigh')
        S_in_bg_RS.connect(j='i')
        
        S_in_target_FS=Synapses(Poisson_target,FS,on_pre='Vinp2=Vhigh')
        S_in_target_FS.connect(j='i')
        S_in_target_RS=Synapses(Poisson_target,RS,on_pre='Vinp2=Vhigh')
        S_in_target_RS.connect(j='i')
    
    #Define monitors and run network :
    R5=SpikeMonitor(RS,record=True)
    R6=SpikeMonitor(FS,record=True)
    
    V_RS=StateMonitor(RS,'V',record=True)
    V_FS=StateMonitor(FS,'V',record=True)
    
    all_neurons=RS,FS,Poisson_background,Poisson_target
    all_synapses=S_FSRS,S_RSFS,S_in_bg_FS,S_in_bg_RS,S_in_target_FS,S_in_target_RS
#    all_synapses=S_FSRS,S_RSFS,S_in_bg_RS,S_in_target_RS
    all_monitors=R5,R6,V_RS,V_FS
    
    return all_neurons,all_synapses,all_monitors



if __name__=='__main__':
    close('all')
    start_scope()    
    
    prefs.codegen.target = 'numpy'
    defaultclock.dt = 0.01*ms
    
    FLee=(0.05*mS/cm**2)/(0.4*uS/cm**2)*0.5
    theta_phase='bad' #'good' or 'bad' or 'mixed'
    runtime=1*second
    target_on=True
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp=0* msiemens * cm **-2
    
    taurinp2=0.1*ms
    taudinp2=0.5*ms
    tauinp2=taudinp2
    taurinp3=0.1*ms
    taudinp3=0.5*ms
    tauinp3=taudinp3
    
    N_FS,N_RS=20,20
    all_neurons,all_synapses,all_monitors=generate_visual_neurons(theta_phase,N_FS,N_RS,runtime,target_on)    
    
    net=Network()
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_monitors)
    
    
#    RS=all_neurons[0]
#    M=StateMonitor(RS,['Iinp1','Iinp2'],record=[0])
#    net.add(M)
    
    prefs.codegen.target = 'cython' #cython=faster, numpy = default python
    
    net.run(runtime,report='text',report_period=300*second)

    R5,R6,V_RS,V_FS=all_monitors
    
    figure()
    plot(R5.t,R5.i+0,'r.',label='RS')
    plot(R6.t,R6.i+20,'b.',label='FS')
    xlim(0,runtime/second)
    legend(loc='upper left',fontsize=12)
    xlabel('Time (s)',fontsize=12)
    ylabel('Neuron index',fontsize=12)
    xticks(fontsize=12)
    yticks(fontsize=12)
    
#    figure()
#    plot(M.t,M.Iinp1[0])
#    plot(M.t,M.Iinp2[0])
    
    clear_cache('cython')