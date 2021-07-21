#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 13:56:08 2020

@author: amelie
"""

from brian2 import *

#from brian2genn import *

from scipy import signal

from FEF_full3_wFS import *
from LIP_full import *

from itertools import *

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

if __name__=='__main__':
    theta_phase='mixed'
    target_on=True  
#    target_on=False 
    target_time = 500*msecond
#    target_time = 625*msecond
#    target_time = 575*msecond #good to bad
#    target_time = 450*msecond #bad to good
#    target_time = 1500*msecond
    
    start_scope()
#    set_device('genn')
    close('all')

    runtime=2*second
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
#    taurinp=1*ms
#    taudinp=5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    
    NN=1 #multiplicative factor on the number of neurons
    N_RS,N_FS,N_SI,N_IB= NN*80,NN*20,NN*20,NN*20 #Neurons for LIP
    N_SI,N_RS_gran,N_SI_gran=20,20,20 #Neurons for LIP
    N_RS_vis,N_FS_vis,N_RS_mot,N_SI_mot,N_dSI_vm,N_RS_vm,N_FS_vm,N_gSI_vm=[20]*8 #Neurons for FEF

    
    all_SIdFSg=[2*msiemens * cm **-2] #1
    all_FSgRSg=[1* msiemens * cm **-2]
    all_RSgFSg=[1*msiemens * cm **-2]
    all_RSgRSg=[0.5*msiemens * cm **-2]
    all_FSgFSg=[0.3* msiemens * cm **-2]
    all_RSgRSs=[2*msiemens * cm **-2]
    all_RSgFSs=[0.1*msiemens * cm **-2]
    all_FSgRSs=[0.1* msiemens * cm **-2]
    all_J_RSg=['10 * uA * cmeter ** -2']
    all_J_FSg=['-5 * uA * cmeter ** -2']
    all_thal=[10* msiemens * cm **-2]
    thal=all_thal[0]
    
    all_syn_cond=list(product(all_SIdFSg,all_FSgRSg,all_RSgFSg,all_RSgRSg,all_FSgFSg,all_RSgRSs,all_RSgFSs,all_FSgRSs))
    all_J=list(product(all_J_RSg,all_J_FSg))
    syn_cond=all_syn_cond[0]
    J=all_J[0]
    
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
    
    all_neurons_FEF,all_synapses_FEF,all_monitors_FEF=create_FEF_full2(N_RS_vis,N_FS_vis,N_RS_mot,N_dSI_vm,N_RS_vm,N_FS_vm,N_gSI_vm,theta_phase,target_on,runtime,target_time)
    RSvm,SIvm,VIPvm,FSvm,V_RS,V_SI,V_VIP,RSv,FSv,VIPv,SIv,RSm=all_monitors_FEF
    RSvm_FEF,SIvm_FEF,RSv_FEF,SIv_FEF,VIPv_FEF=all_neurons_FEF[1],all_neurons_FEF[3],all_neurons_FEF[7],all_neurons_FEF[10],all_neurons_FEF[9]
    
    all_neurons_LIP, all_synapses_LIP, all_gap_junctions_LIP, all_monitors_LIP=make_full_network(syn_cond,J,thal,theta_phase)
    V1,V2,V3,R1,R2,R3,I1,I2,I3,V4,R4,I4s,I4a,I4ad,I4bd,R5,R6,R7,V5,V6,V7,inpmon,inpIBmon=all_monitors_LIP
    RS_sup_LIP,IB_LIP,SI_deep_LIP=all_neurons_LIP[0],all_neurons_LIP[5],all_neurons_LIP[9]
    RS_gran_LIP,FS_gran_LIP=all_neurons_LIP[7],all_neurons_LIP[8]
    
    IB_LIP.ginp_IB=0* msiemens * cm **-2 #the input to RS_sup_LIP is provided with synapses from FEF 
    SI_deep_LIP.ginp_SI=0* msiemens * cm **-2
    RSvm_FEF.ginp_RS=0* msiemens * cm **-2
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
#    S_LIP_FS_FEF=generate_syn(RS_sup_LIP,SIvm_FEF,'Isyn_LIP','',0.009*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
#    S_LIP_RS_FEF=generate_syn(RS_sup_LIP,RSvm_FEF,'Isyn_LIP','',0.015*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    S_LIP_FS_FEF=generate_syn(RS_sup_LIP,SIvm_FEF,'Isyn_LIP','',0.009*msiemens * cm **-2,0.125*ms,1*ms,0*mV)  
#    S_LIP_FS_FEF=generate_syn(RS_sup_LIP,SIvm_FEF,'Isyn_LIP','',0.015*msiemens * cm **-2,0.125*ms,1*ms,0*mV)  
    
    S_LIP_RSv_FEF=generate_syn(RS_sup_LIP,RSv_FEF,'Isyn_LIP','',0.010*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    S_LIP_SIv_FEF=generate_syn(RS_sup_LIP,SIv_FEF,'Isyn_LIP','',0.020*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    S_LIP_VIPv_FEF=generate_syn(RS_sup_LIP,VIPv_FEF,'Isyn_LIP','',0.005*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
   
    if target_on:
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
#    set_device('genn')
#    defaultclock.dt = 0.01*ms
    
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    taurinp2=2*ms
    taudinp2=10*ms
    tauinp2=taudinp2
    net.run(runtime,report='text',report_period=300*second)
    
    # LIP Plots
    figure()
    plot(R1.t,R1.i+140,'r.',label='RS cells')
    plot(R2.t,R2.i+120,'b.',label='FS cells')
    plot(R3.t,R3.i+100,'g.',label='SI cells')
    plot(R5.t,R5.i+70,'.',label='Granular RS',color='C1')
    plot(R6.t,R6.i+50,'c.',label='Granular FS')
    plot(R4.t,R4.i+20,'m.',label='IB cells')
    plot(R7.t,R7.i,'.',label='Deep SI',color='lime')
    xlim(0.2,runtime/second)
    legend(loc='upper left')
    xlabel('Time (s)')
    ylabel('Neuron index')
    
    min_t=int(100*ms*100000*Hz)
    LFP_LIP=1/80*sum(V1.V,axis=0)[min_t:]
    LFP_FEF=1/20*sum(V_RS.V,axis=0)[min_t:]
    
    record_dt=1/512*second
    t=int(0.3*second/record_dt) #t_debut
    L=int(2*second/record_dt)
    fs = 1/record_dt
    freq = linspace(1/second, fs/2, 100)
    widths = 6*fs/(2*freq*pi)
            
    
    figure()
    subplot(121)
    CWT = signal.cwt(LFP_LIP, signal.morlet2, widths, w=6)
    #f, t, Sxx = signal.spectrogram(LFP_LIP, 100000*Hz,nperseg=30000,noverlap=25000)
    pcolormesh(V1.t[min_t:], freq, CWT, cmap='RdBu')#, shading='gouraud')
    ylabel('Frequency [Hz]')
    xlabel('Time [sec]')
    ylim(0,50)
    
    subplot(122)
    CWT = signal.cwt(LFP_FEF, signal.morlet2, widths, w=6)
    #f, t, Sxx = signal.spectrogram(LFP_FEF, 100000*Hz,nperseg=30000,noverlap=25000)
    pcolormesh(V_RS.t[min_t:], freq, CWT, cmap='RdBu')#, shading='gouraud')
    ylabel('Frequency [Hz]')
    xlabel('Time [sec]')
    ylim(0,50)
    
    savefig(sim_dir+'/spec.eps')
    
#    
#    f,Spectrum_LFP_LIP=signal.periodogram(LFP_LIP, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_FEF=signal.periodogram(LFP_FEF, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_FS=signal.periodogram(LFP_V_FS, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_SI=signal.periodogram(LFP_V_SI, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_IB=signal.periodogram(LFP_V_IB, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_RSg=signal.periodogram(LFP_V_RSg, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_FSg=signal.periodogram(LFP_V_FSg, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_SId=signal.periodogram(LFP_V_SId, 100000,'flattop', scaling='spectrum')
#    
#    figure(figsize=(10,8))    
#    subplot(421)
#    plot(f,Spectrum_LFP_V_RS)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('RS cell')
#    subplot(422)
#    plot(f,Spectrum_LFP_V_FS)
#    yticks([],[])
#    xlim(0,100)
#    title('FS cell')
#    subplot(423)
#    plot(f,Spectrum_LFP_V_SI)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('SI cell')
#    subplot(425)
#    plot(f,Spectrum_LFP_V_RSg)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('gran RS cell')
#    subplot(426)
#    plot(f,Spectrum_LFP_V_FSg)
#    yticks([],[])
#    xlim(0,100)
#    title('gran FS cell')
#    subplot(427)
#    plot(f,Spectrum_LFP_V_IB)
#    xlabel('Frequency (Hz)')
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('IB cell')
#    subplot(428)
#    plot(f,Spectrum_LFP_V_SId)
#    yticks([],[])
#    xlim(0,100)
#    xlabel('Frequency (Hz)')
#    title('deep SI cell')
#    
#    tight_layout()
#    
    
    
    #FEF Plots    
    figure()#figsize=(10,4))
    subplot(311)
    ylabel('Visual Neurons')
    plot(RSv.t,RSv.i+60,'r.',label='RS')
    plot(FSv.t,FSv.i+40,'b.',label='FS')
    plot(VIPv.t,VIPv.i+0,'k.',label='VIP')
    plot(SIv.t,SIv.i+20,'g.',label='SOM')
    xlim(0.2,runtime/second)
    legend(loc='upper left')   
    xlabel('Time (s)')
    #ylabel('Neuron index')
    
    subplot(312)
    ylabel('Visual-Motor Neurons')
    plot(VIPvm.t,VIPvm.i+0,'k.',label='VIP')
    plot(RSvm.t,RSvm.i+60,'r.',label='RS')
    plot(FSvm.t,FSvm.i+40,'b.',label='FS')
    plot(SIvm.t,SIvm.i+20,'g.',label='SOM')
    xlim(0.2,runtime/second)
    #legend(loc='upper left') 
    xlabel('Time (s)')
    #ylabel('Neuron index')
    
    subplot(313)
    ylabel('Decision cells')
    plot(RSm.t,RSm.i+0,'r.',label='RS')
    xlim(0.2,runtime/second)
    #legend(loc='upper left') 
    xlabel('Time (s)')
    #ylabel('Neuron index')
    
    tight_layout()
    
    
    # LIP Plots
    figure(figsize=(9,9))
#    subplot(411)
    up=240
    plot(R1.t,R1.i+140+up,'r.',label='RS cells')
    plot(R2.t,R2.i+120+up,'b.',label='FS cells')
    plot(R3.t,R3.i+100+up,'g.',label='SI cells')
    plot([0.2,runtime/second],[95+up,95+up],'k--')
    plot(R5.t,R5.i+70+up,'r.',label='Granular RS')
    plot(R6.t,R6.i+50+up,'b.',label='Granular FS')
    plot([0.2,runtime/second],[45+up,45+up],'k--')
    plot(R4.t,R4.i+20+up,'m.',label='IB cells')
    plot(R7.t,R7.i+up,'g.',label='Deep SI')
    xlim(0.2,runtime/second)
    plot([0.2,runtime/second],[up-10,up-10],'k')
#    xticks([],[])
    #legend(loc='upper left')
    #xlabel('Time (s)')
#    ylabel('Neuron index')
#    title('LIP')
    
#    subplot(413)
    up=40
#    title('FEF Visual Neurons')
    plot(RSv.t,RSv.i+60+up,'r.',label='RS')
    plot(FSv.t,FSv.i+40+up,'b.',label='FS')
    plot(VIPv.t,VIPv.i+0+up,'k.',label='VIP')
    plot(SIv.t,SIv.i+20+up,'g.',label='SOM')
    xlim(0.2,runtime/second)
    plot([0.2,runtime/second],[up-10,up-10],'k')
#    xticks([],[])
#    legend(loc='upper left')   
    #xlabel('Time (s)')
#    ylabel('Neuron index')
    
#    subplot(412)
    up=140
#    title('FEF Visual-Motor Neurons')
    plot(VIPvm.t,VIPvm.i+0+up,'k.',label='VIP')
    plot(RSvm.t,RSvm.i+60+up,'r.',label='RS')
    plot(FSvm.t,FSvm.i+40+up,'b.',label='FS')
    plot(SIvm.t,SIvm.i+20+up,'g.',label='SOM')
    xlim(0.2,runtime/second)
    plot([0.2,runtime/second],[up-10,up-10],'k')
#    xticks([],[])
#    legend(loc='upper left') 
    #xlabel('Time (s)')
#    ylabel('Neuron index')
    
#    subplot(414)
#    title('FEF Decision cells')
    plot(RSm.t,RSm.i+0,'r.',label='RS')
    xlim(0.2,runtime/second)
#    legend(loc='upper left') 
    xlabel('Time (s)')
    ylabel('Neuron index')
    
    tight_layout()    
    subplots_adjust(bottom=0.1)
    
#    subplot(133)
#    title('Motor Neurons')
#    plot(R6.t,R6.i+60,'r.',label='RS')
#    plot(R7.t,R7.i+40,'b.',label='SI')
#    plot(R8.t,R8.i+0,'c.',label='Fix')
#    xlim(0,runtime/second)
#    legend(loc='upper left') 
    
    
#    min_t=int(50*ms*100000*Hz)
#    LFP_V1=1/20*sum(V1FEF.V,axis=0)[min_t:]
#    LFP_V2=1/20*sum(V2FEF.V,axis=0)[min_t:]
#    LFP_V3=1/20*sum(V3FEF.V,axis=0)[min_t:]
#    LFP_V4=1/20*sum(V4FEF.V,axis=0)[min_t:]
#    LFP_V5=1/20*sum(V5FEF.V,axis=0)[min_t:]
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
#    subplot(321)
#    plot(f,Spectrum_LFP_V4)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('visual RS')
#    subplot(323)
#    plot(f,Spectrum_LFP_V5)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('visual FS')  
#    
#    subplot(322)
#    plot(f,Spectrum_LFP_V1)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('visual-motor gran RS')
#    subplot(324)
#    plot(f,Spectrum_LFP_V2)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('visual-motor gran SI')  
#    subplot(326)
#    plot(f,Spectrum_LFP_V3)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('visual-motor deep SI')    
    
    
    # LIP Plots
    figure(figsize=(9,9))
#    subplot(411)
    up=100
    plot(R1.t,R1.i+140+up,'r.',label='RS cells')
    plot(R2.t,R2.i+120+up,'b.',label='FS cells')
    plot(R3.t,R3.i+100+up,'g.',label='SI cells')
    plot([0.2,runtime/second],[95+up,95+up],'k--')
    plot(R5.t,R5.i+70+up,'r.',label='Granular RS')
    plot(R6.t,R6.i+50+up,'b.',label='Granular FS')
    plot([0.2,runtime/second],[45+up,45+up],'k--')
    plot(R4.t,R4.i+20+up,'m.',label='IB cells')
    plot(R7.t,R7.i+up,'g.',label='Deep SI')
    xlim(0.2,runtime/second)
    plot([0.2,runtime/second],[up-10,up-10],'k')

#    subplot(412)
    up=0
#    title('FEF Visual-Motor Neurons')
    plot(VIPvm.t,VIPvm.i+0+up,'k.',label='VIP')
    plot(RSvm.t,RSvm.i+60+up,'r.',label='RS')
    plot(FSvm.t,FSvm.i+40+up,'b.',label='FS')
    plot(SIvm.t,SIvm.i+20+up,'g.',label='SOM')
    xlim(0.2,runtime/second)
#    plot([0.2,runtime/second],[up-10,up-10],'k')
#    xticks([],[])
#    legend(loc='upper left') 
    #xlabel('Time (s)')
#    ylabel('Neuron index')
    xlabel('Time (s)')
    ylabel('Neuron index')
    
    
    figure(figsize=(5,5))
#    subplot(411)
    up=0
    plot(R1.t,R1.i+140+up,'r.',label='RS cells')
    plot(R2.t,R2.i+120+up,'b.',label='FS cells')
    plot(R3.t,R3.i+100+up,'g.',label='SI cells')
    plot([0.2,runtime/second],[95+up,95+up],'k--')
    plot(R5.t,R5.i+70+up,'r.',label='Granular RS')
    plot(R6.t,R6.i+50+up,'b.',label='Granular FS')
    plot([0.2,runtime/second],[45+up,45+up],'k--')
    plot(R4.t,R4.i+20+up,'m.',label='IB cells')
    plot(R7.t,R7.i+up,'g.',label='Deep SI')
    xlim(0.2,runtime/second)

    xlabel('Time (s)')
    ylabel('Neuron index')
    
    #clear_cache('cython')
