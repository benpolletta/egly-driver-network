# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 14:28:44 2019

@author: aaussel
"""

from brian2 import *
from scipy import signal
from cells.RS_LIP_altd import *
from cells.FS_LIP_altd import *
from cells.SI_LIP import *
from cells.IB_soma_LIP import *
from cells.IB_axon_LIP import *
from cells.IB_apical_dendrite_LIP import *
from cells.IB_basal_dendrite_LIP import *

from LIP_superficial_layer import *

def create_Mark_Alex_network(kainate,version,Nf=1):

    prefs.codegen.target = 'numpy'
    
    defaultclock.dt = 0.01*ms

    N_RS,N_FS,N_SI,N_IB= Nf*80,Nf*20,Nf*20,Nf*20 #Number of neurons of RE, TC, and HTC type
    
    
    all_neurons,all_synapses,all_gap_junctions,all_monitors=create_superficial_layer(kainate,version,Nf)
    RS, FS, SI=all_neurons
    #Single column network
    
    ##Define neuron groups 
    IB_soma=NeuronGroup(N_IB,eq_IB_soma,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    IB_soma.V = '-100*mvolt+10*rand()*mvolt'
    IB_soma.h = '0+0.05*rand()'
    IB_soma.m = '0+0.05*rand()'
    IB_soma.J='-4.5 * uA * cmeter ** -2' #article SI=-3.5, code=-4.5
    
    IB_axon=NeuronGroup(N_IB,eq_IB_axon,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    IB_axon.V = '-100*mvolt+10*rand()*mvolt'
    IB_axon.h = '0+0.05*rand()'
    IB_axon.m = '0+0.05*rand()'
    IB_axon.mKM = '0+0.05*rand()'
    if kainate=='low':
        IB_axon.J='-0.4 * uA * cmeter ** -2' #article SI=+0.1, code=-0.4
    elif kainate=='high':
        IB_axon.J='-6 * uA * cmeter ** -2' #article SI=+0.1, code=-0.4
    
    IB_ad=NeuronGroup(N_IB,eq_IB_ad,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    IB_ad.V = '-100*mvolt+10*rand()*mvolt'
    IB_ad.h = '0+0.05*rand()'
    IB_ad.m = '0+0.05*rand()'
    IB_ad.mAR = '0+0.001*rand()'
    IB_ad.mKM = '0+0.05*rand()'
    IB_ad.mCaH = '0+0.01*rand()'
    if kainate=='low':
        IB_ad.J='25.5 * uA * cmeter ** -2'  #article SI=27.5, code=25.5
    elif kainate=='high':
        IB_ad.J='23.5 * uA * cmeter ** -2'  #article SI=27.5, code=25.5
    
    IB_bd=NeuronGroup(N_IB,eq_IB_bd,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    IB_bd.V = '-100*mvolt+10*rand()*mvolt'
    IB_bd.h = '0+0.05*rand()'
    IB_bd.m = '0+0.05*rand()'
    IB_bd.mAR = '0+0.001*rand()'
    IB_bd.mKM = '0+0.05*rand()'
    IB_bd.mCaH = '0+0.01*rand()'
    if kainate=='low':
        IB_bd.J='42.5 * uA * cmeter ** -2' #article SI=44.5, code=42.5
    elif kainate=='high':
        IB_bd.J='23.5 * uA * cmeter ** -2' #article SI=44.5, code=42.5
    
    
    ##Synapses
    eq_syn='''_post=s_i*g_i*(V_post-V_i) : amp * meter ** -2 (summed)
        ds_i/dt=-s_i/taud_i+(1-s_i)/taur_i*0.5*(1+tanh(V_pre/10/mV)) : 1
        g_i : siemens * meter**-2
        V_i : volt
        taud_i : second
        taur_i : second
    '''
    
    S_RSIB_AMPA=None
    S_RSIB_NMDA=None
    if version=='Alex' or True: #RS -> IB connections in Mark version not explicit
        S_RSIB_AMPA=Synapses(RS,IB_ad,model='IsynRS_LIP_sup_AMPA'+eq_syn,method='exact')
        S_RSIB_AMPA.connect(j='i%20+k for k in range(-1,3)', skip_if_invalid=True)
        S_RSIB_AMPA.g_i=1/60* msiemens * cm **-2
    #    S_RSIB_AMPA.g_i=0.1* msiemens * cm **-2
        S_RSIB_AMPA.taur_i=0.125*ms
        S_RSIB_AMPA.taud_i=1*ms
        S_RSIB_AMPA.V_i=0*mV
        
        S_RSIB_NMDA=Synapses(RS,IB_ad,model='IsynRS_LIP_sup_NMDA'+eq_syn,method='exact')
        S_RSIB_NMDA.connect(j='i%20+k for k in range(-1,3)', skip_if_invalid=True)
        S_RSIB_NMDA.g_i=1/240* msiemens * cm **-2
    #    S_RSIB_NMDA.g_i=0.05* msiemens * cm **-2
        S_RSIB_NMDA.taur_i=12.5*ms
        S_RSIB_NMDA.taud_i=125*ms
        S_RSIB_NMDA.V_i=0*mV
    
    
    S_SIIB=Synapses(SI,IB_ad,model='IsynSI_LIP_sup'+eq_syn,method='exact')
    S_SIIB.connect()
    S_SIIB.g_i=0.4* msiemens * cm **-2*int(version=='Alex')+(4*msiemens * cm **-2)*int(version=='Mark')+(0.2*msiemens * cm **-2)*int(version=='Mark/cell')
    #S_SIIB.g_i=4*msiemens * cm **-2
    S_SIIB.taur_i=0.25*ms
    S_SIIB.taud_i=20*ms
    S_SIIB.V_i=-80*mV
    
    S_IBFS=Synapses(IB_axon,FS,model='IsynIB_LIP'+eq_syn,method='exact')
    S_IBFS.connect()
    S_IBFS.g_i=0.2* msiemens * cm **-2*int(version=='Alex')+(0.045*msiemens * cm **-2)*int(version=='Mark')+(0.00225*msiemens * cm **-2)*int(version=='Mark/cell')
    S_IBFS.g_i=0.08* msiemens * cm **-2
#    S_IBFS.g_i=0.1* msiemens * cm **-2
    S_IBFS.taur_i=0.125*ms
    S_IBFS.taud_i=1*ms
    S_IBFS.V_i=0*mV
    
    S_IBSI=Synapses(IB_axon,SI,model='IsynIB_LIP'+eq_syn,method='exact')
    S_IBSI.connect()
    S_IBSI.g_i=0.045* msiemens * cm **-2*int(version=='Alex')+(0.04*msiemens * cm **-2)*int(version=='Mark')+(0.002*msiemens * cm **-2)*int(version=='Mark/cell')
#    S_IBSI.g_i=0.08* msiemens * cm **-2
    S_IBSI.taur_i=1.25*ms
    S_IBSI.taud_i=50*ms
    S_IBSI.V_i=0*mV
    
    S_IBIB=None
    if version=='Alex': #NMDA synapses from IB to IB also exist in Mark's model of beta1
        S_IBIB=Synapses(IB_axon,IB_bd,model='IsynIB_LIP'+eq_syn,method='exact')
        S_IBIB.connect()
        S_IBIB.g_i=1/500* msiemens * cm **-2
        S_IBIB.taur_i=0.25*ms  #0.5 in Mark
        S_IBIB.taud_i=100*ms
        S_IBIB.V_i=0*mV
    
    
    eq_gap='''_post=g_i*(V_post-V_pre) : amp * meter ** -2 (summed)
        g_i : siemens * meter**-2
    '''
    
    gapIB_SomaAd=Synapses(IB_soma,IB_ad,model='Igap_soma'+eq_gap,method='exact')
    gapIB_SomaAd.connect(j='i')
    gapIB_SomaAd.g_i=0.2* msiemens * cm **-2
    
    gapIB_SomaBd=Synapses(IB_soma,IB_bd,model='Igap_soma'+eq_gap,method='exact')
    gapIB_SomaBd.connect(j='i')
    gapIB_SomaBd.g_i=0.2* msiemens * cm **-2
    
    gapIB_SomaAxon=Synapses(IB_soma,IB_axon,model='Igap_soma'+eq_gap,method='exact')
    gapIB_SomaAxon.connect(j='i')
    gapIB_SomaAxon.g_i=0.3* msiemens * cm **-2
    
    gapIB_AdSoma=Synapses(IB_ad,IB_soma,model='Igap_ad'+eq_gap,method='exact')
    gapIB_AdSoma.connect(j='i')
    gapIB_AdSoma.g_i=0.4* msiemens * cm **-2
    
    gapIB_BdSoma=Synapses(IB_bd,IB_soma,model='Igap_bd'+eq_gap,method='exact')
    gapIB_BdSoma.connect(j='i')
    gapIB_BdSoma.g_i=0.4* msiemens * cm **-2
    
    gapIB_AxonSoma=Synapses(IB_axon,IB_soma,model='Igap_axon'+eq_gap,method='exact')
    gapIB_AxonSoma.connect(j='i')
    gapIB_AxonSoma.g_i=0.3* msiemens * cm **-2
    
    gap_IBIB=Synapses(IB_axon,IB_axon,model='Igap_axon'+eq_gap,method='exact')
    gap_IBIB.connect()
    gap_IBIB.g_i=0.0025* msiemens * cm **-2*int(version=='Alex')+(0.002*msiemens * cm **-2)*int(version=='Mark')

    #Define monitors :
    V4=StateMonitor(IB_soma,'V',record=True)
    R4=SpikeMonitor(IB_soma,record=True)
    I4s=StateMonitor(IB_soma,'Isyn',record=True)
    I4a=StateMonitor(IB_axon,'Isyn',record=True)
    I4ad=StateMonitor(IB_ad,'Isyn',record=True)
    I4bd=StateMonitor(IB_bd,'Isyn',record=True)

    new_neurons=IB_soma,IB_axon,IB_bd,IB_ad
    new_synapses=S_RSIB_AMPA,S_RSIB_NMDA, S_SIIB, S_IBFS, S_IBSI, S_IBIB
    new_synapses=tuple([y for y in new_synapses if y])
    new_gap_junctions=gapIB_SomaAd,gapIB_SomaBd,gapIB_SomaAxon,gapIB_AdSoma,gapIB_BdSoma,gapIB_AxonSoma,gap_IBIB
    new_monitors=V4,R4,I4s,I4a,I4ad,I4bd
    
    all_neurons=all_neurons+new_neurons
    all_synapses=all_synapses+new_synapses
    all_gap_junctions=all_gap_junctions+new_gap_junctions
    all_monitors=all_monitors+new_monitors
    
    return all_neurons, all_synapses, all_gap_junctions, all_monitors
    

#Other inputs (depends on simulation)
if __name__=='__main__' :
    
    start_scope()
    version = 'Alex'
    kainate='low'
    
    runtime=1*second
    f=25*Hz #rythmic input frequency
    top_down_binding=False #Figure 2
    bottom_up=False #Figures 5e,f, 6
    distractors_IB=False #Figures 5a,c,e
    distractors_RS=False #Figures 5b,d,f
    top_down_dishinibition=False #Figures 6c,d
    random_init_RS=False #All figures
    random_init_IB=False #All figures except 1e and 2g
    
    NN=1 #multiplicative factor on the number of neurons
    N_RS,N_FS,N_SI,N_IB= NN*80,NN*20,NN*20,NN*20 #Number of neurons of RE, TC, and HTC type
    
    net = Network()
    all_neurons, all_synapses, all_gap_junctions, all_monitors=create_Mark_Alex_network(kainate,version,Nf=NN)
    V1,V2,V3,R1,R2,R3,I1,I2,I3,V4,R4,I4s,I4a,I4ad,I4bd=all_monitors

    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_gap_junctions)
    net.add(all_monitors)
    
    input_beta=False
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp_IB=0* msiemens * cm **-2
    ginp_SI=0* msiemens * cm **-2
    ginp=0* msiemens * cm **-2
    
    def generate_spike_timing(N,f,start_time,end_time=runtime):
        list_time_and_i=[]
        for i in range(N):
            list_time=[(start_time,i)]
            next_spike=list_time[-1][0]+(1+0.1*rand())/f
            while next_spike<end_time:
                list_time.append((next_spike,i))
                next_spike=list_time[-1][0]+(1+0.1*rand())/f
            list_time_and_i+=list_time
        return array(list_time_and_i)
    
    if input_beta:
        ginp_IB=2* msiemens * cm **-2
        inputs_topdown=generate_spike_timing(N_IB,f,0*ms,end_time=2100*ms)
        print(inputs_topdown)
        G_topdown = SpikeGeneratorGroup(N_IB, inputs_topdown[:,1], inputs_topdown[:,0]*second)
        topdown_in=Synapses(G_topdown,IB_bd,on_pre='Vinp=Vhigh')
        topdown_in.connect(j='i')
        
    #    RS.ginp_RS=2* msiemens * cm **-2
    #    inputs_topdown2=generate_spike_timing(N_RS,f,0*ms,end_time=2100*ms)
    #    print(inputs_topdown2)
    #    G_topdown2 = SpikeGeneratorGroup(N_RS, inputs_topdown2[:,1], inputs_topdown2[:,0]*second)
    #    topdown_in2=Synapses(G_topdown2,RS,on_pre='Vinp=Vhigh')
    #    topdown_in2.connect(j='i')
    
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
    #subplot(224)
    #plot(V1.t/second,V4.V[0]/volt)
    #xlabel('Time (s)')
    #ylabel('Membrane potential (V)')
    #title('IB cell')
    
    #figure()
    #subplot(411)
    #plot(R1.t,R1.i,'r.')
    #xlim(0,runtime/second)
    #title('RS cell')
    #subplot(412)
    #plot(R2.t,R2.i,'r.')
    #xlim(0,runtime/second)
    #title('FS cell')
    #subplot(413)
    #plot(R3.t,R3.i,'r.')
    #xlim(0,runtime/second)
    #title('SI cell')
    #subplot(414)
    #plot(R4.t,R4.i,'r.')
    #xlim(0,runtime/second)
    #title('IB cell')
    
    figure()
    plot(R1.t,R1.i+60,'r.',label='RS cells')
    plot(R2.t,R2.i+40,'b.',label='FS cells')
    plot(R3.t,R3.i+20,'g.',label='SI cells')
    plot(R4.t,R4.i,'y.',label='IB cells')
    xlim(0,runtime/second)
    legend(loc='upper left')
    #
    min_t=int(50*ms*100000*Hz)
    LFP_V_RS=1/N_RS*sum(V1.V,axis=0)[min_t:]
    LFP_V_FS=1/N_FS*sum(V2.V,axis=0)[min_t:]
    LFP_V_SI=1/N_SI*sum(V3.V,axis=0)[min_t:]
    LFP_V_IB=1/N_IB*sum(V4.V,axis=0)[min_t:]
    
    f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_FS=signal.periodogram(LFP_V_FS, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_SI=signal.periodogram(LFP_V_SI, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_IB=signal.periodogram(LFP_V_IB, 100000,'flattop', scaling='spectrum')
    
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
    plot((V1.t/second)[min_t:],LFP_V_IB)
    xlabel('Time (s)')
    ylabel('LFP')
    title('IB cell')
    
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
    plot(f,Spectrum_LFP_V_IB)
    xlabel('Frequency (Hz)')
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,50)
    title('IB cell')
    #
    #figure()
    #LFP_V_sup=1/3*(LFP_V_RS+LFP_V_FS+LFP_V_SI)
    #f,Spectrum_LFP_V_sup=signal.periodogram(LFP_V_sup, 100000,'flattop', scaling='spectrum')
    #subplot(221)
    #plot((V1.t/second)[min_t:],LFP_V_sup)
    #xlabel('Time (s)')
    #ylabel('LFP')
    #title('Superficial layer')
    #subplot(223)
    #plot((V1.t/second)[min_t:],LFP_V_IB)
    #xlabel('Time (s)')
    #ylabel('LFP')
    #title('Deep layer')
    #subplot(222)
    #plot(f,Spectrum_LFP_V_sup)
    #xlabel('Frequency (Hz)')
    #ylabel('Spectrum')
    #yticks([],[])
    #xlim(0,50)
    #title('Superficial layer')
    #subplot(224)
    #plot(f,Spectrum_LFP_V_IB)
    #xlabel('Frequency (Hz)')
    #ylabel('Spectrum')
    #yticks([],[])
    #xlim(0,50)
    #title('Deep layer')
    #
    #
    #
    #LFP_I_RS=1/N_RS*sum(I1.Isyn,axis=0)[min_t:]
    #LFP_I_FS=1/N_FS*sum(I2.Isyn,axis=0)[min_t:]
    #LFP_I_SI=1/N_SI*sum(I3.Isyn,axis=0)[min_t:]
    #LFP_I_IB=1/N_IB*(sum(I4s.Isyn,axis=0)[min_t:]+sum(I4a.Isyn,axis=0)[min_t:]+sum(I4ad.Isyn,axis=0)[min_t:]+sum(I4bd.Isyn,axis=0)[min_t:])
    #
    #f,Spectrum_LFP_I_RS=signal.periodogram(LFP_I_RS, 100000,'flattop', scaling='spectrum')
    #f,Spectrum_LFP_I_FS=signal.periodogram(LFP_I_FS, 100000,'flattop', scaling='spectrum')
    #f,Spectrum_LFP_I_SI=signal.periodogram(LFP_I_SI, 100000,'flattop', scaling='spectrum')
    #f,Spectrum_LFP_I_IB=signal.periodogram(LFP_I_IB, 100000,'flattop', scaling='spectrum')
    #
    #
    #figure()
    #subplot(421)
    #plot((V1.t/second)[min_t:],LFP_I_RS)
    #ylabel('LFP')
    #title('RS cell')
    #subplot(423)
    #plot((V1.t/second)[min_t:],LFP_I_FS)
    #ylabel('LFP')
    #title('FS cell')
    #subplot(425)
    #plot((V1.t/second)[min_t:],LFP_I_SI)
    #ylabel('LFP')
    #title('SI cell')
    #subplot(427)
    #plot((V1.t/second)[min_t:],LFP_I_IB)
    #xlabel('Time (s)')
    #ylabel('LFP')
    #title('IB cell')
    #
    #subplot(422)
    #plot(f,Spectrum_LFP_I_RS)
    #ylabel('Spectrum')
    #yticks([],[])
    #xlim(0,50)
    #title('RS cell')
    #subplot(424)
    #plot(f,Spectrum_LFP_I_FS)
    #ylabel('Spectrum')
    #yticks([],[])
    #xlim(0,50)
    #title('FS cell')
    #subplot(426)
    #plot(f,Spectrum_LFP_I_SI)
    #ylabel('Spectrum')
    #yticks([],[])
    #xlim(0,50)
    #title('SI cell')
    #subplot(428)
    #plot(f,Spectrum_LFP_I_IB)
    #xlabel('Frequency (Hz)')
    #ylabel('Spectrum')
    #yticks([],[])
    #xlim(0,50)
    #title('IB cell')
    #
    #figure()
    #LFP_I_sup=1/3*(LFP_I_RS+LFP_I_FS+LFP_I_SI)
    #f,Spectrum_LFP_I_sup=signal.periodogram(LFP_I_sup, 100000,'flattop', scaling='spectrum')
    #subplot(221)
    #plot((V1.t/second)[min_t:],LFP_I_sup)
    #xlabel('Time (s)')
    #ylabel('LFP')
    #title('Superficial layer')
    #subplot(223)
    #plot((V1.t/second)[min_t:],LFP_I_IB)
    #xlabel('Time (s)')
    #ylabel('LFP')
    #title('Deep layer')
    #subplot(222)
    #plot(f,Spectrum_LFP_I_sup)
    #xlabel('Frequency (Hz)')
    #ylabel('Spectrum')
    #yticks([],[])
    #xlim(0,50)
    #title('Superficial layer')
    #subplot(224)
    #plot(f,Spectrum_LFP_I_IB)
    #xlabel('Frequency (Hz)')
    #ylabel('Spectrum')
    #yticks([],[])
    #xlim(0,50)
    #title('Deep layer')
    
    
    clear_cache('cython')