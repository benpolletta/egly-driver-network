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
from cells.VIP_LIP import *
from cells.IB_soma_LIP import *
from cells.IB_axon_LIP import *
from cells.IB_apical_dendrite_LIP import *
from cells.IB_basal_dendrite_LIP import *

def create_superficial_layer(kainate,version,Nf=1):
    
    prefs.codegen.target = 'numpy'
    
    defaultclock.dt = 0.01*ms
    
    #Single column network
    
    ##Define neuron groups
    N_RS,N_FS,N_SI,N_VIP= Nf*80,Nf*20,Nf*20,Nf*20 #Number of neurons of RE, TC, and HTC type
    
    RS=NeuronGroup(N_RS,eq_RS_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
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
        
    FS=NeuronGroup(N_FS,eq_FS_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    FS.V = '-110*mvolt+10*rand()*mvolt'
    FS.h = '0+0.05*rand()'
    FS.m = '0+0.05*rand()'
    if kainate=='low':
        FS.J='35 * uA * cmeter ** -2' #article=code=35
    elif kainate=='high':
        FS.J='16 * uA * cmeter ** -2'
    
    SI=NeuronGroup(N_SI,eq_SI_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
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
    
    VIP=NeuronGroup(N_VIP,eq_VIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    VIP.V = '-90*mvolt+10*rand()*mvolt'
    VIP.Iapp='4 * uA * cmeter ** -2' #article=code=35   

    
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
    
    S_RSRS=None
    
    rsfs_g_i=(1/40* msiemens * cm **-2)*int(version=='Alex')+(1*msiemens * cm **-2)*int(version=='Mark')+(0.05*msiemens * cm **-2)*int(version=='Mark/cell')
    S_RSFS=generate_syn(RS,FS,'IsynRS_LIP_sup','i//40==j//10',2*rsfs_g_i,0.125*ms,1*ms,0*mV)
    
    rssi_g_i=0.225* msiemens * cm **-2*int(version=='Alex')+(2*msiemens * cm **-2)*int(version=='Mark')+(0.1*msiemens * cm **-2)*int(version=='Mark/cell')
    S_RSSI=generate_syn(RS,SI,'IsynRS_LIP_sup','i//40==j//10',2*rssi_g_i,1.25*ms,1*ms,0*mV)
    
    fsrs_g_i=6.25* msiemens * cm **-2*int(version=='Alex')+(25*msiemens * cm **-2)*int(version=='Mark')+(1.25*msiemens * cm **-2)*int(version=='Mark/cell')
    S_FSRS=generate_syn(FS,RS,'IsynFS_LIP_sup','i//10==j//40',2*fsrs_g_i,0.25*ms,5*ms,-80*mV)
    
    fsfs_g_i=2* msiemens * cm **-2*int(version=='Alex')+(20*msiemens * cm **-2)*int(version=='Mark')+(1*msiemens * cm **-2)*int(version=='Mark/cell')
    S_FSFS=generate_syn(FS,FS,'IsynFS_LIP_sup','j==i',fsfs_g_i,0.25*ms,5*ms,-75*mV)
    
    fssi_g_i=0.4* msiemens * cm **-2*int(version=='Alex')+8* msiemens * cm **-2*int(version=='Mark')+0.4* msiemens * cm **-2*int(version=='Mark/cell')
    S_FSSI=generate_syn(FS,SI,'IsynFS_LIP_sup','i//10==j//10',2*fssi_g_i,0.25*ms,6*ms,-80*mV)
    
    sirs_g_i=0.125* msiemens * cm **-2*int(version=='Alex')+2.5* msiemens * cm **-2*int(version=='Mark')+0.125* msiemens * cm **-2*int(version=='Mark/cell')
    S_SIRS=generate_syn(SI,RS,'IsynSI_LIP_sup','i//10==j//40',2*sirs_g_i,0.25*ms,20*ms,-80*mV)
    
    S_SIFS=None
    if version=='Alex' or True:
        S_SIFS=generate_syn(SI,FS,'IsynSI_LIP_sup','i//10==j//10',2*0.2* msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
    
    sisi_g_i=7* msiemens * cm **-2*int(version=='Alex')+(5*msiemens * cm **-2)*int(version=='Mark')+(0.25*msiemens * cm **-2)*int(version=='Mark/cell')
    S_SISI=generate_syn(SI,SI,'IsynSI_LIP_sup','j==i',sisi_g_i,0.25*ms,20*ms,-80*mV)
    
    #### Synapses (taken from FEF visual module).

    S_VIPSI=generate_syn(VIP,SI,'IsynVIP_LIP_sup','i//10==j//10',0.7* msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
#    S_VIPSI=generate_syn(VIP,SI,'IsynVIP_LIP_sup','i//10==j//10',0.01* msiemens * cm **-2,0.25*ms,20*ms,-80*mV) 
    S_VIPFS=generate_syn(VIP,FS,'IsynVIP_LIP_sup','i//10==j//10',0*msiemens*cm**-2,0.25*ms,20*ms,-80*mV)
    S_SIVIP=generate_syn(SI,VIP,'IsynSI_LIP_sup','',0.01* msiemens * cm **-2,0.25*ms,20*ms,-80*mV) 
    S_VIPVIP=generate_syn(VIP,VIP,'IsynVIP_LIP_sup','i//10==j//10',0*msiemens*cm**-2,0.25*ms,20*ms,-80*mV)
    
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
            
#     #### Inputs (taken from FEF visual module). 
            
#     if input_FEF:    
#         FS.ginp_FS=thal_cond*0
#         RS.ginp_RS=7.5* msiemens * cm **-2
#         SI.ginp_SI=7.5* msiemens * cm **-2
#         VIP.ginp_VIP_good=3* msiemens * cm **-2
#         VIP.ginp_VIP_bad=3* msiemens * cm **-2
#         if theta_phase=='good':
#             fFEF=25*Hz
#         else :
#             fFEF=0*Hz
            
#         gamma_background=generate_spike_timing(N_FS,fFEF,0*ms,end_time=3000*ms)
        
#         if theta_phase=='mixed':
#             t0=0*ms
#             t1=125*ms
#             gamma_background=generate_spike_timing(N_FS,fLIP,t0,end_time=t1)
#             while t0+125*ms<runtime:
#                 fLIP=50*Hz*int(fLIP==13*Hz)+13*Hz*int(fLIP==50*Hz)
#                 t0,t1=t0+125*ms,t1+125*ms
#                 gamma_background=vstack((gamma_background,generate_spike_timing(N_FS,fLIP,t0,end_time=t1)))
                
#         gamma_target=generate_spike_timing(10,50*Hz,target_time,end_time=target_time+100*ms)
        
#         Poisson_background = SpikeGeneratorGroup(N_FS, gamma_background[:,1], gamma_background[:,0]*second)
#         if target_on :
# #            Poisson_target = SpikeGeneratorGroup(20, gamma_target[:,1], gamma_target[:,0]*second)
#             Poisson_target = SpikeGeneratorGroup(10, gamma_target[:,1], gamma_target[:,0]*second)
#         else :
# #            Poisson_target = SpikeGeneratorGroup(20, [], []*second)
#             Poisson_target = SpikeGeneratorGroup(10, [], []*second)
# #        S_in_bg_FS=Synapses(Poisson_background,FS,on_pre='Vinp=Vhigh')
# #        S_in_bg_FS.connect(j='i')
#         S_in_bg_RS=Synapses(Poisson_background,RS,on_pre='Vinp=Vhigh')
#         S_in_bg_RS.connect(j='i')
#         S_in_bg_SI=Synapses(Poisson_background,SI,on_pre='Vinp=Vhigh')
#         S_in_bg_SI.connect(j='i')
#         S_in_bg_VIP=Synapses(Poisson_background,VIP,on_pre='Vinp=Vhigh')
#         S_in_bg_VIP.connect(j='i')
        
# #        S_in_target_FS=Synapses(Poisson_target,FS,on_pre='Vinp2=Vhigh')
# #        S_in_target_FS.connect(j='i')
# #        S_in_target_RS=Synapses(Poisson_target,RS,on_pre='Vinp2=Vhigh')
# #        S_in_target_RS.connect(j='i')
#         S_in_target_VIP=Synapses(Poisson_target,VIP,on_pre='Vinp2=Vhigh')
#         S_in_target_VIP.connect(j='i')
#         S_in_target_SI=Synapses(Poisson_target,SI,on_pre='Vinp2=Vhigh')
#         S_in_target_SI.connect(j='i')
#         SI.ginp_SI2=2.5* msiemens * cm **-2
#         VIP.ginp_VIP2=2.5* msiemens * cm **-2
# ##        SI.ginp_SI2=3* msiemens * cm **-2
# ##        VIP.ginp_VIP2=3* msiemens * cm **-2
#         RS.ginp_RS2=2.5* msiemens * cm **-2
    
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
    version = 'Alex' #'Alex' or 'Mark' or 'Mark/cell'
    #Mark means the exact values from Mark's article SI, Mark/cell means the one where values have been divided by the number of presynaptic cells
    kainate='low' #'low or 'high'
    
    #Other inputs (depends on simulation)
    runtime=1*second
    f=25*Hz #rythmic input frequency
    
    input_FEF=False
    theta_phase='good'
    
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
    N_RS,N_FS,N_SI,N_VIP= NN*80,NN*20,NN*20,NN*20 #Number of neurons of RE, TC, and HTC type
    
    net = Network(collect())
    all_neurons,all_synapses,all_gap_junctions,all_monitors=create_superficial_layer(kainate,version,Nf=NN)
    
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_gap_junctions)
    net.add(all_monitors)
    
    V1,V2,V3,V4,R1,R2,R3,R4,I1,I2,I3,I4=all_monitors
    
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
