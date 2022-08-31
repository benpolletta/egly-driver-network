# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 14:28:44 2019

@author: aaussel
"""
#This is used for simulations that change the time constants of all SOM and/or FS interneurons at the same time.

from brian2 import *
from scipy import signal
from model_files.cells.RS_LIP import *
from model_files.cells.FS_LIP import *
from model_files.cells.SI_LIP import *
from model_files.cells.IB_soma_LIP import *
from model_files.cells.IB_axon_LIP import *
from model_files.cells.IB_apical_dendrite_LIP import *
from model_files.cells.IB_basal_dendrite_LIP import *


def create_superficial_layer(t_SI,t_FS,Nf=1):
    #Defines the superficial layer of LIP
    #Mostly adapted from the papers by Mark Kramer and Alex Gelastopoulos
    
    prefs.codegen.target = 'numpy'
    
    defaultclock.dt = 0.01*ms
    
    
    ##Define neuron groups
    N_RS,N_FS,N_SI,N_IB= Nf*80,Nf*20,Nf*20,Nf*20 #Number of neurons of RE, TC, and HTC type
    
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
    
    ##Define synapses
    eq_syn='''_post=s_i*g_i*(V_post-V_i) : amp * meter ** -2 (summed)
        ds_i/dt=-s_i/taud_i+(1-s_i)/taur_i*0.5*(1+tanh(V_pre/10/mV)) : 1
        g_i : siemens * meter**-2
        V_i : volt
        taud_i : second
        taur_i : second
    '''
    
    S_RSRS=None
    
    S_RSFS=Synapses(RS,FS,model='IsynRS_LIP_sup'+eq_syn,method='exact')
    S_RSFS.connect()
    S_RSFS.g_i=1/40* msiemens * cm **-2
    S_RSFS.taur_i=0.125*ms
    S_RSFS.taud_i=1*ms
    S_RSFS.V_i=0*mV
    
    S_RSSI=Synapses(RS,SI,model='IsynRS_LIP_sup'+eq_syn,method='exact')
    S_RSSI.connect()
    S_RSSI.g_i=0.225* msiemens * cm **-2
    S_RSSI.taur_i=1.25*ms
    S_RSSI.taud_i=1*ms
    S_RSSI.V_i=0*mV
    
    
    S_FSRS=Synapses(FS,RS,model='IsynFS_LIP_sup'+eq_syn,method='exact')
    S_FSRS.connect()
    S_FSRS.g_i=6.25* msiemens * cm **-2
    S_FSRS.taur_i=0.25*ms
    S_FSRS.taud_i=t_FS
    S_FSRS.V_i=-80*mV
    
    S_FSFS=Synapses(FS,FS,model='IsynFS_LIP_sup'+eq_syn,method='exact')
    S_FSFS.connect(j='i')
    S_FSFS.g_i=2* msiemens * cm **-2
    S_FSFS.taur_i=0.25*ms
    S_FSFS.taud_i=t_FS
    S_FSFS.V_i=-75*mV
    
    S_FSSI=Synapses(FS,SI,model='IsynFS_LIP_sup'+eq_syn,method='exact')
    S_FSSI.connect()
    S_FSSI.g_i=0.4* msiemens * cm **-2
    S_FSSI.taur_i=0.25*ms
    S_FSSI.taud_i=t_FS
    S_FSSI.V_i=-80*mV
    
    S_SIRS=Synapses(SI,RS,model='IsynSI_LIP_sup'+eq_syn,method='exact')
    S_SIRS.connect()
    S_SIRS.g_i=2* msiemens * cm **-2
    S_SIRS.taur_i=0.25*ms
    S_SIRS.taud_i=t_SI
    S_SIRS.V_i=-80*mV
    
    S_SIFS=Synapses(SI,FS,model='IsynSI_LIP_sup'+eq_syn,method='exact')
    S_SIFS.connect()
    S_SIFS.g_i=0.2* msiemens * cm **-2
    S_SIFS.taur_i=0.25*ms
    S_SIFS.taud_i=t_SI
    S_SIFS.V_i=-80*mV
    
    S_SISI=Synapses(SI,SI,model='IsynSI_LIP_sup'+eq_syn,method='exact')
    S_SISI.connect(j='i')
    S_SISI.g_i=7* msiemens * cm **-2
    S_SISI.taur_i=0.25*ms
    S_SISI.taud_i=t_SI
    S_SISI.V_i=-80*mV
    
    ##Define gap junctions
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
    
    R1=SpikeMonitor(RS,record=True)
    R2=SpikeMonitor(FS,record=True)
    R3=SpikeMonitor(SI,record=True)
    
    I1=StateMonitor(RS,'Isyn',record=True)
    I2=StateMonitor(FS,'Isyn',record=True)
    I3=StateMonitor(SI,'Isyn',record=True)
    
    all_neurons=RS, FS, SI
    all_synapses=S_RSRS, S_RSFS, S_RSSI, S_FSRS, S_FSFS, S_FSSI, S_SIRS, S_SIFS, S_SISI
    all_synapses=tuple([y for y in all_synapses if y])
    all_gap_junctions=gap_SISI, gap_RSRS
    all_gap_junctions=tuple([y for y in all_gap_junctions if y])
    all_monitors=V1,V2,V3,R1,R2,R3,I1,I2,I3

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
    N_RS,N_FS,N_SI,N_IB= NN*80,NN*20,NN*20,NN*20 #Number of neurons of each type
    
    net = Network(collect())
    all_neurons,all_synapses,all_gap_junctions,all_monitors=create_superficial_layer(t_SI,t_FS,Nf=NN)
    
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_gap_junctions)
    net.add(all_monitors)
    
    V1,V2,V3,R1,R2,R3,I1,I2,I3=all_monitors
    
    prefs.codegen.target = 'cython'  #cython=faster, numpy = default python

    
    net.run(runtime,report='text',report_period=300*second)
    
    
    figure()
    plot(R1.t,R1.i+40,'r.',label='RS cells')
    plot(R2.t,R2.i+20,'b.',label='FS cells')
    plot(R3.t,R3.i,'g.',label='SI cells')
    xlim(0,runtime/second)
    legend(loc='upper left')
    
    min_t=int(50*ms*100000*Hz)
    LFP_V_RS=1/N_RS*sum(V1.V,axis=0)[min_t:]
    LFP_V_FS=1/N_FS*sum(V2.V,axis=0)[min_t:]
    LFP_V_SI=1/N_SI*sum(V3.V,axis=0)[min_t:]
    
    f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_FS=signal.periodogram(LFP_V_FS, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_SI=signal.periodogram(LFP_V_SI, 100000,'flattop', scaling='spectrum')
    
    figure()
    subplot(321)
    plot((V1.t/second)[min_t:],LFP_V_RS)
    ylabel('LFP')
    title('RS cell')
    subplot(323)
    plot((V1.t/second)[min_t:],LFP_V_FS)
    ylabel('LFP')
    title('FS cell')
    subplot(325)
    plot((V1.t/second)[min_t:],LFP_V_SI)
    ylabel('LFP')
    title('SI cell')
    
    subplot(322)
    plot(f,Spectrum_LFP_V_RS)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,50)
    title('RS cell')
    subplot(324)
    plot(f,Spectrum_LFP_V_FS)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,50)
    title('FS cell')
    subplot(326)
    plot(f,Spectrum_LFP_V_SI)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,50)
    title('SI cell')
    
