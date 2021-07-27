# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 14:28:44 2019

@author: aaussel
"""

from brian2 import *
from scipy import signal
from cells.RS_LIP import *
from cells.FS_LIP import *
from cells.SI_LIP import *
from cells.IB_soma_LIP import *
from cells.IB_axon_LIP import *
from cells.IB_apical_dendrite_LIP import *
from cells.IB_basal_dendrite_LIP import *

def create_superficial_layer(kainate,version,Nf=1):
    
    prefs.codegen.target = 'numpy'
    
    defaultclock.dt = 0.01*ms
    
    #Single column network
    
    ##Define neuron groups
    N_RS,N_FS,N_SI,N_IB= Nf*80,Nf*20,Nf*20,Nf*20 #Number of neurons of RE, TC, and HTC type
    
    RS=NeuronGroup(N_RS,eq_RS_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4',name='RSsupLIP')
    RS.V = '-70*mvolt+10*rand()*mvolt'
    RS.h = '0+0.05*rand()'
    RS.m = '0+0.05*rand()'
    RS.mAR = '0.035+0.025*rand()'
#    RS.V = '-60*mvolt'
#    RS.h = '0.56'
#    RS.m = '0.038'
#    RS.mAR = '0.01'
#    RS.J='1 * uA * cmeter ** -2'
    if kainate=='low':
        RS.J='1 * uA * cmeter ** -2'  #article SI=25, code=1
    elif kainate=='high':
        RS.J='-10 * uA * cmeter ** -2'  #article SI=25, code=1
        
    FS=NeuronGroup(N_FS,eq_FS_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4',name='FSsupLIP')
    FS.h = '0+0.05*rand()'
    FS.m = '0+0.05*rand()'
    if kainate=='low':
        FS.J='35 * uA * cmeter ** -2' #article=code=35
    elif kainate=='high':
        FS.J='16 * uA * cmeter ** -2'
    
    SI=NeuronGroup(N_SI,eq_SI_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4',name='SIsupLIP')
    SI.V = '-100*mvolt+10*rand()*mvolt'
    SI.h = '0+0.05*rand()'
    SI.m = '0+0.05*rand()'
    SI.mAR = '0.02+0.04*rand()'
    if version == 'Alex':
        if kainate=='low':
            SI.J='35* uA * cmeter ** -2' #article SI=50, code=35, Mark = 45
    #        SI.J='45* uA * cmeter ** -2'
        elif kainate=='high':
            SI.J='30* uA * cmeter ** -2' #article SI=50, code=35, Mark = 45
    elif version=='Mark':
        if kainate=='low':
            SI.J='45* uA * cmeter ** -2' #article SI=50, code=35, Mark = 45
        elif kainate=='high':
            SI.J='40* uA * cmeter ** -2' #article SI=50, code=35, Mark = 45

    
    ##Synapses
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
    S_RSFS.g_i=(1/40* msiemens * cm **-2)*int(version=='Alex')+(1*msiemens * cm **-2)*int(version=='Mark')+(0.05*msiemens * cm **-2)*int(version=='Mark/cell')
#    S_RSFS.g_i=0.01*msiemens * cm **-2
    S_RSFS.taur_i=0.125*ms
    S_RSFS.taud_i=1*ms
    S_RSFS.V_i=0*mV
    
    S_RSSI=Synapses(RS,SI,model='IsynRS_LIP_sup'+eq_syn,method='exact')
    S_RSSI.connect()
    S_RSSI.g_i=0.225* msiemens * cm **-2*int(version=='Alex')+(2*msiemens * cm **-2)*int(version=='Mark')+(0.1*msiemens * cm **-2)*int(version=='Mark/cell')
    #S_RSSI.g_i=0* msiemens * cm **-2
    S_RSSI.taur_i=1.25*ms
    S_RSSI.taud_i=1*ms
    S_RSSI.V_i=0*mV
    
    
    S_FSRS=Synapses(FS,RS,model='IsynFS_LIP_sup'+eq_syn,method='exact')
    S_FSRS.connect()
    S_FSRS.g_i=6.25* msiemens * cm **-2*int(version=='Alex')+(25*msiemens * cm **-2)*int(version=='Mark')+(1.25*msiemens * cm **-2)*int(version=='Mark/cell')
#    S_FSRS.g_i=6* msiemens * cm **-2
    S_FSRS.taur_i=0.25*ms
    S_FSRS.taud_i=5*ms
    S_FSRS.V_i=-80*mV
    
    S_FSFS=Synapses(FS,FS,model='IsynFS_LIP_sup'+eq_syn,method='exact')
    S_FSFS.connect(j='i')
    S_FSFS.g_i=2* msiemens * cm **-2*int(version=='Alex')+(20*msiemens * cm **-2)*int(version=='Mark')+(1*msiemens * cm **-2)*int(version=='Mark/cell')
    #S_FSFS.g_i=2* msiemens * cm **-2
    S_FSFS.taur_i=0.25*ms
    S_FSFS.taud_i=5*ms
    S_FSFS.V_i=-75*mV
    
    S_FSSI=Synapses(FS,SI,model='IsynFS_LIP_sup'+eq_syn,method='exact')
    S_FSSI.connect()
    S_FSSI.g_i=0.4* msiemens * cm **-2*int(version=='Alex')+8* msiemens * cm **-2*int(version=='Mark')+0.4* msiemens * cm **-2*int(version=='Mark/cell')
    #S_FSSI.g_i=0.1* msiemens * cm **-2
    S_FSSI.taur_i=0.25*ms
    S_FSSI.taud_i=6*ms
    S_FSSI.V_i=-80*mV
    
    S_SIRS=Synapses(SI,RS,model='IsynSI_LIP_sup'+eq_syn,method='exact')
    S_SIRS.connect()
    S_SIRS.g_i=0.125* msiemens * cm **-2*int(version=='Alex')+2.5* msiemens * cm **-2*int(version=='Mark')+0.125* msiemens * cm **-2*int(version=='Mark/cell')
    S_SIRS.g_i=2* msiemens * cm **-2
    #S_SIRS.g_i=0* msiemens * cm **-2
    S_SIRS.taur_i=0.25*ms
    S_SIRS.taud_i=20*ms
    S_SIRS.V_i=-80*mV
    
    S_SIFS=None
    if version=='Alex' or True:
        S_SIFS=Synapses(SI,FS,model='IsynSI_LIP_sup'+eq_syn,method='exact')
        S_SIFS.connect()
        S_SIFS.g_i=0.2* msiemens * cm **-2
    #    S_SIFS.g_i=0.1* msiemens * cm **-2
    #    S_SIFS.g_i=0* msiemens * cm **-2
        S_SIFS.taur_i=0.25*ms
        S_SIFS.taud_i=20*ms
        S_SIFS.V_i=-80*mV
    
    S_SISI=Synapses(SI,SI,model='IsynSI_LIP_sup'+eq_syn,method='exact')
    S_SISI.connect(j='i')
    S_SISI.g_i=7* msiemens * cm **-2*int(version=='Alex')+(5*msiemens * cm **-2)*int(version=='Mark')+(0.25*msiemens * cm **-2)*int(version=='Mark/cell')
    #S_SISI.g_i=5* msiemens * cm **-2
    S_SISI.taur_i=0.25*ms
    S_SISI.taud_i=20*ms
    S_SISI.V_i=-80*mV
    
    
    eq_gap='''_post=g_i*(V_post-V_pre) : amp * meter ** -2 (summed)
        g_i : siemens * meter**-2
    '''
    
    gap_SISI=None
    if version=='Alex':
        gap_SISI=Synapses(SI,SI,model='Igap'+eq_gap,method='exact')
        gap_SISI.connect()
        gap_SISI.g_i=0.2* msiemens * cm **-2
    
    gap_RSRS=None
    if version=='Mark' or version=='Mark/cell' or True:
        gap_RSRS=Synapses(RS,RS,model='Igap'+eq_gap,method='exact')
        gap_RSRS.connect()
        gap_RSRS.g_i=0.04* msiemens * cm **-2    
    #    gap_RSRS.g_i=0.01* msiemens * cm **-2   
    
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
    version = 'Alex' #'Alex' or 'Mark' or 'Mark/cell'
    #Mark means the exact values from Mark's article SI, Mark/cell means the one where values have been divided by the number of presynaptic cells
    kainate='low' #'low or 'high'
    
    #Other inputs (depends on simulation)
    runtime=1*second
    f=25*Hz #rythmic input frequency
    top_down_binding=False #Figure 2
    bottom_up=False #Figures 5e,f, 6
    distractors_IB=False #Figures 5a,c,e
    distractors_RS=False #Figures 5b,d,f
    top_down_dishinibition=False #Figures 6c,d
    random_init_RS=False #All figures
    random_init_IB=False #All figures except 1e and 2g
    
    input_beta=False
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp_IB=0* msiemens * cm **-2
    ginp_SI=10* msiemens * cm **-2
    ginp=0* msiemens * cm **-2
    
    NN=1 #multiplicative factor on the number of neurons
    N_RS,N_FS,N_SI,N_IB= NN*80,NN*20,NN*20,NN*20 #Number of neurons of RE, TC, and HTC type
    
    net = Network(collect())
    all_neurons,all_synapses,all_gap_junctions,all_monitors=create_superficial_layer(kainate,version,Nf=NN)
    
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_gap_junctions)
    net.add(all_monitors)
    
    V1,V2,V3,R1,R2,R3,I1,I2,I3=all_monitors
    
    prefs.codegen.target = 'cython'  #cython=faster, numpy = default python

    
    net.run(runtime,report='text',report_period=300*second)
    
    #figure()
    #subplot(221)
    #plot(V1.t/second,V1.V[0]/volt)
    #xlabel('Time (s)')
    #ylabel('Membrane potential (V)')
    #title('RS cell')
    #subplot(222)
    #plot(V1.t/second,V2.V[0]/volt)
    #xlabel('Time (s)')
    #ylabel('Membrane potential (V)')
    #title('FS cell')
    #subplot(223)
    #plot(V1.t/second,V3.V[0]/volt)
    #xlabel('Time (s)')
    #ylabel('Membrane potential (V)')
    #title('SI cell')
    
    #figure()
    #subplot(311)
    #plot(R1.t,R1.i,'r.')
    #xlim(0,runtime/second)
    #title('RS cell')
    #subplot(312)
    #plot(R2.t,R2.i,'r.')
    #xlim(0,runtime/second)
    #title('FS cell')
    #subplot(313)
    #plot(R3.t,R3.i,'r.')
    #xlim(0,runtime/second)
    #title('SI cell')
    
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
    
    
    #LFP_I_RS=1/N_RS*sum(I1.Isyn,axis=0)[min_t:]
    #LFP_I_FS=1/N_FS*sum(I2.Isyn,axis=0)[min_t:]
    #LFP_I_SI=1/N_SI*sum(I3.Isyn,axis=0)[min_t:]
    #
    #f,Spectrum_LFP_I_RS=signal.periodogram(LFP_I_RS, 100000,'flattop', scaling='spectrum')
    #f,Spectrum_LFP_I_FS=signal.periodogram(LFP_I_FS, 100000,'flattop', scaling='spectrum')
    #f,Spectrum_LFP_I_SI=signal.periodogram(LFP_I_SI, 100000,'flattop', scaling='spectrum')
    #
    #figure()
    #subplot(321)
    #plot((V1.t/second)[min_t:],LFP_I_RS)
    #ylabel('LFP')
    #title('RS cell')
    #subplot(323)
    #plot((V1.t/second)[min_t:],LFP_I_FS)
    #ylabel('LFP')
    #title('FS cell')
    #subplot(325)
    #plot((V1.t/second)[min_t:],LFP_I_SI)
    #ylabel('LFP')
    #title('SI cell')
    #
    #subplot(322)
    #plot(f,Spectrum_LFP_I_RS)
    #ylabel('Spectrum')
    #yticks([],[])
    #xlim(0,50)
    #title('RS cell')
    #subplot(324)
    #plot(f,Spectrum_LFP_I_FS)
    #ylabel('Spectrum')
    #yticks([],[])
    #xlim(0,50)
    #title('FS cell')
    #subplot(326)
    #plot(f,Spectrum_LFP_I_SI)
    #ylabel('Spectrum')
    #yticks([],[])
    #xlim(0,50)
    #title('SI cell')
