# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 14:28:44 2019

@author: aaussel
"""
#This is used for simulations that change the time constants of all SOM and/or FS interneurons at the same time.

from brian2 import *
from scipy import signal
from cells.RS_LIP import *
from cells.FS_LIP import *
from cells.SI_LIP import *
from cells.VIP_LIP import *
from cells.IB_soma_LIP import *
from cells.IB_axon_LIP import *
from cells.IB_apical_dendrite_LIP import *
from cells.IB_basal_dendrite_LIP import *


def create_superficial_layer(t_SI,t_FS,Nf=1):
    #Defines the superficial layer of LIP
    #Mostly adapted from the papers by Mark Kramer and Alex Gelastopoulos
    
    prefs.codegen.target = 'numpy'
    
    defaultclock.dt = 0.01*ms
    
    
    ##Define neuron groups
    N_RS,N_FS,N_SI,N_VIP= Nf*80,Nf*20,Nf*20,Nf*20 #Number of neurons of RE, TC, and HTC type
    
    #RS cells
    RS=NeuronGroup(N_RS,eq_RS_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    RS.V = '-70*mvolt+10*rand()*mvolt'
    RS.h = '0+0.05*rand()'
    RS.m = '0+0.05*rand()'
    RS.mAR = '0.035+0.025*rand()'
    RS.J='1 * uA * cmeter ** -2'
    
    #FS cells
    FS=NeuronGroup(N_FS,eq_FS_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    FS.V = '-110*mvolt+10*rand()*mvolt'
    FS.h = '0+0.05*rand()'
    FS.m = '0+0.05*rand()'
    FS.J='35 * uA * cmeter ** -2'
    
    #SOM cells
    SI=NeuronGroup(N_SI,eq_SI_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    SI.V = '-100*mvolt+10*rand()*mvolt'
    SI.h = '0+0.05*rand()'
    SI.m = '0+0.05*rand()'
    SI.mAR = '0.02+0.04*rand()'
    SI.J='35* uA * cmeter ** -2' 
    
    VIP=NeuronGroup(N_VIP,eq_VIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    VIP.V = '-90*mvolt+10*rand()*mvolt'
    VIP.Iapp='5 * uA * cmeter ** -2' #article=code=35   
    
    ##Define synapses
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
    
    S_RSRS=None
    
    rsfs_g_i=1/40* msiemens * cm **-2
    rssi_g_i=0.225* msiemens * cm **-2
    fsrs_g_i=6.25* msiemens * cm **-2
    fsfs_g_i=2* msiemens * cm **-2
    fssi_g_i=0.4* msiemens * cm **-2
    sirs_g_i=2* msiemens * cm **-2
    sifs_g_i=0.2* msiemens * cm **-2
    sisi_g_i=7* msiemens * cm **-2

    S_RSFS=generate_syn(RS,FS,'IsynRS_LIP_sup','i//40==j//10',2*rsfs_g_i,0.125*ms,1*ms,0*mV)
    S_RSSI=generate_syn(RS,SI,'IsynRS_LIP_sup','i//40==j//10',2*rssi_g_i,1.25*ms,1*ms,0*mV)
    
    S_FSRS=generate_syn(FS,RS,'IsynFS_LIP_sup','i//10==j//40',2*fsrs_g_i,0.25*ms,t_FS,-80*mV)
    S_FSFS=generate_syn(FS,FS,'IsynFS_LIP_sup','j==i',2*fsfs_g_i,0.25*ms,t_FS,-75*mV)
    S_FSSI=generate_syn(FS,SI,'IsynFS_LIP_sup','i//10==j//10',2*fssi_g_i,0.25*ms,t_FS,-80*mV)
    
    S_SIRS=generate_syn(SI,RS,'IsynSI_LIP_sup','i//10==j//40',2*sirs_g_i,0.25*ms,t_SI,-80*mV)
    S_SIFS=generate_syn(SI,FS,'IsynSI_LIP_sup','i//10==j//10',2*sifs_g_i,0.25*ms,t_SI,-80*mV)
    S_SISI=generate_syn(SI,SI,'IsynSI_LIP_sup','j==i',2*sisi_g_i,0.25*ms,t_SI,-80*mV)
    
    #### Synapses (taken from FEF visual module).

#    S_VIPSI=generate_syn(VIP,SI,'IsynSI_FEF_VM','i//10==j//10',0.7* msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
    S_VIPSI=generate_syn(VIP,SI,'IsynVIP_LIP_sup','i//10==j//10',1.5* msiemens * cm **-2,0.25*ms,20*ms,-80*mV) 
    S_VIPFS=generate_syn(VIP,FS,'IsynVIP_LIP_sup','i//10==j//10',0.1*msiemens*cm**-2,0.25*ms,20*ms,-80*mV)
    S_SIVIP=generate_syn(SI,VIP,'IsynSI_LIP_sup','',0.01* msiemens * cm **-2,0.25*ms,20*ms,-80*mV) 
    S_VIPVIP=generate_syn(VIP,VIP,'IsynVIP_LIP_sup','i//10==j//10',0.1*msiemens*cm**-2,0.25*ms,20*ms,-80*mV)
    
    eq_gap='''_post=g_i*(V_post-V_pre) : amp * meter ** -2 (summed)
        g_i : siemens * meter**-2
    '''
    
    gap_SISI=Synapses(SI,SI,model='Igap'+eq_gap,method='exact')
    gap_SISI.connect()
    gap_SISI.g_i=0.2* msiemens * cm **-2
    
    gap_RSRS=Synapses(RS,RS,model='Igap'+eq_gap,method='exact')
    gap_RSRS.connect()
    gap_RSRS.g_i=0.04* msiemens * cm **-2    

    
    #Define monitors:
    V1=StateMonitor(RS,'V',record=True)
    V2=StateMonitor(FS,'V',record=True)
    V3=StateMonitor(SI,'V',record=True)
    V4=StateMonitor(VIP,'V',record=True)
    
    
    R1=SpikeMonitor(RS,record=True)
    R2=SpikeMonitor(FS,record=True)
    R3=SpikeMonitor(SI,record=True)
    R4=SpikeMonitor(VIP,record=True)
    
    I1=StateMonitor(RS,'Isyn',record=True)
    I2=StateMonitor(FS,'Isyn',record=True)
    I3=StateMonitor(SI,'Isyn',record=True)
    I4=StateMonitor(VIP,'Isyn',record=True)
    
    all_neurons=RS, FS, SI, VIP
    all_synapses=S_RSRS, S_RSFS, S_RSSI, S_FSRS, S_FSFS, S_FSSI, S_SIRS, S_SIFS, S_SISI, S_VIPSI, S_SIVIP, S_VIPVIP
    all_synapses=tuple([y for y in all_synapses if y])
    all_gap_junctions=gap_SISI, gap_RSRS
    all_gap_junctions=tuple([y for y in all_gap_junctions if y])
    all_monitors=V1,V2,V3,V4,R1,R2,R3,R4,I1,I2,I3,I4
    #all_monitors=V1,V2,V3,R1,R2,R3,I1,I2,I3

    return all_neurons,all_synapses,all_gap_junctions,all_monitors    

if __name__=='__main__' :
    start_scope()
    
    runtime=1*second
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp_IB=0* msiemens * cm **-2
    ginp_SI=10* msiemens * cm **-2
    ginp=0* msiemens * cm **-2
    
    t_SI,t_FS=20*msecond,5*msecond
    
    NN=1 #multiplicative factor on the number of neurons
    N_RS,N_FS,N_SI,N_VIP= NN*80,NN*20,NN*20,NN*20 #Number of neurons of RE, TC, and HTC type
    
    net = Network(collect())
    all_neurons,all_synapses,all_gap_junctions,all_monitors=create_superficial_layer(t_SI,t_FS,Nf=NN)
    
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_gap_junctions)
    net.add(all_monitors)
    
    V1,V2,V3,V4,R1,R2,R3,R4,I1,I2,I3,I4=all_monitors
    
    prefs.codegen.target = 'cython'  #cython=faster, numpy = default python

    
    net.run(runtime,report='text',report_period=300*second)
    
    
    figure()
    plot(R1.t,R1.i+60,'r.',label='RS cells')
    plot(R2.t,R2.i+40,'b.',label='FS cells')
    plot(R3.t,R3.i+20,'g.',label='SI cells')
    plot(R4.t,R4.i,'k.',label='VIP cells')
    xlim(0,runtime/second)
    legend(loc='upper left')
    
    min_t=int(50*ms*100000*Hz)
    LFP_V_RS=1/N_RS*sum(V1.V,axis=0)[min_t:]
    LFP_V_FS=1/N_FS*sum(V2.V,axis=0)[min_t:]
    LFP_V_SI=1/N_SI*sum(V3.V,axis=0)[min_t:]
    LFP_V_VIP=1/N_VIP*sum(V4.V,axis=0)[min_t:]
    
    f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_FS=signal.periodogram(LFP_V_FS, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_SI=signal.periodogram(LFP_V_SI, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_VIP=signal.periodogram(LFP_V_VIP, 100000,'flattop', scaling='spectrum')
    
    figure()
    subplot(421)
    plot((V1.t/second)[min_t:],LFP_V_RS)
    ylabel('LFP')
    title('RS cell')
    subplot(423)
    plot((V1.t/second)[min_t:],LFP_V_FS)
    ylabel('LFP')
    title('FS cell')
    subplot(425)
    plot((V1.t/second)[min_t:],LFP_V_SI)
    ylabel('LFP')
    title('SI cell')
    subplot(427)
    plot((V1.t/second)[min_t:],LFP_V_VIP)
    ylabel('LFP')
    title('VIP cell')
    
    subplot(422)
    plot(f,Spectrum_LFP_V_RS)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,50)
    title('RS cell')
    subplot(424)
    plot(f,Spectrum_LFP_V_FS)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,50)
    title('FS cell')
    subplot(426)
    plot(f,Spectrum_LFP_V_SI)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,50)
    title('SI cell')
    subplot(428)
    plot(f,Spectrum_LFP_V_VIP)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,50)
    title('VIP cell')
    
