# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 14:28:44 2019

@author: aaussel
"""

#This is used for simulations that change the time constants of all interneurons of one type at the same time.


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

from LIP_superficial_layer import *

def create_beta1_network(t_SI,t_FS,Nf=1):
    #Defines the subnetwork of LIP responsible from the beta1 rhythm
    #Mostly adapted from the papers by Mark Kramer and Alex Gelastopoulos
    
    prefs.codegen.target = 'numpy'
    
    defaultclock.dt = 0.01*ms

    N_RS,N_FS,N_SI,N_IB= Nf*80,Nf*20,Nf*20,Nf*20 #Number of neurons of RE, TC, and HTC type
    
    print(t_SI,t_FS)
    
    all_neurons,all_synapses,all_gap_junctions,all_monitors=create_superficial_layer(kainate,version,Nf)
    RS, FS, SI, VIP=all_neurons
    
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
    IB_axon.J='-0.4 * uA * cmeter ** -2' #article SI=+0.1, code=-0.4
    
    IB_ad=NeuronGroup(N_IB,eq_IB_ad,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    IB_ad.V = '-100*mvolt+10*rand()*mvolt'
    IB_ad.h = '0+0.05*rand()'
    IB_ad.m = '0+0.05*rand()'
    IB_ad.mAR = '0+0.001*rand()'
    IB_ad.mKM = '0+0.05*rand()'
    IB_ad.mCaH = '0+0.01*rand()'
    IB_ad.J='25.5 * uA * cmeter ** -2'  #article SI=27.5, code=25.5
    
    IB_bd=NeuronGroup(N_IB,eq_IB_bd,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    IB_bd.V = '-100*mvolt+10*rand()*mvolt'
    IB_bd.h = '0+0.05*rand()'
    IB_bd.m = '0+0.05*rand()'
    IB_bd.mAR = '0+0.001*rand()'
    IB_bd.mKM = '0+0.05*rand()'
    IB_bd.mCaH = '0+0.01*rand()'
    IB_bd.J='42.5 * uA * cmeter ** -2' #article SI=44.5, code=42.5
    
    
    ##Define synapses
    eq_syn='''_post=s_i*g_i*(V_post-V_i) : amp * meter ** -2 (summed)
        ds_i/dt=-s_i/taud_i+(1-s_i)/taur_i*0.5*(1+tanh(V_pre/10/mV)) : 1
        g_i : siemens * meter**-2
        V_i : volt
        taud_i : second
        taur_i : second
    '''
    
    S_RSIB_AMPA=Synapses(RS,IB_ad,model='IsynRS_LIP_sup_AMPA'+eq_syn,method='exact')
    S_RSIB_AMPA.connect(j='i%20+k for k in range(-1,3)', skip_if_invalid=True)
    S_RSIB_AMPA.g_i=1/60* msiemens * cm **-2
    S_RSIB_AMPA.taur_i=0.125*ms
    S_RSIB_AMPA.taud_i=1*ms
    S_RSIB_AMPA.V_i=0*mV
    
    S_RSIB_NMDA=Synapses(RS,IB_ad,model='IsynRS_LIP_sup_NMDA'+eq_syn,method='exact')
    S_RSIB_NMDA.connect(j='i%20+k for k in range(-1,3)', skip_if_invalid=True)
    S_RSIB_NMDA.g_i=1/240* msiemens * cm **-2
    S_RSIB_NMDA.taur_i=12.5*ms
    S_RSIB_NMDA.taud_i=125*ms
    S_RSIB_NMDA.V_i=0*mV
    
    
    S_SIIB=Synapses(SI,IB_ad,model='IsynSI_LIP_sup'+eq_syn,method='exact')
    S_SIIB.connect('i//10==j//10')
    S_SIIB.g_i=2*0.4* msiemens * cm **-2
    S_SIIB.taur_i=0.25*ms
    S_SIIB.taud_i=t_SI
    S_SIIB.V_i=-80*mV
    
    S_IBFS=Synapses(IB_axon,FS,model='IsynIB_LIP'+eq_syn,method='exact')
    S_IBFS.connect('i//10==j//10')
    S_IBFS.g_i=2*0.08* msiemens * cm **-2
    S_IBFS.taur_i=0.125*ms
    S_IBFS.taud_i=1*ms
    S_IBFS.V_i=0*mV
    
    S_IBSI=Synapses(IB_axon,SI,model='IsynIB_LIP'+eq_syn,method='exact')
    S_IBSI.connect('i//10==j//10')
    S_IBSI.g_i=2*0.045* msiemens * cm **-2
    S_IBSI.taur_i=1.25*ms
    S_IBSI.taud_i=50*ms
    S_IBSI.V_i=0*mV
    
    S_IBIB=Synapses(IB_axon,IB_bd,model='IsynIB_LIP'+eq_syn,method='exact')
    S_IBIB.connect('i//10==j//10')
    S_IBIB.g_i=2*1/500* msiemens * cm **-2
    S_IBIB.taur_i=0.25*ms  #0.5 in Mark
    S_IBIB.taud_i=100*ms
    S_IBIB.V_i=0*mV

    S_IBVIP=Synapses(IB_axon,VIP,model='IsynIB_LIP'+eq_syn,method='exact')
    S_IBVIP.connect('i//10==j//10')
#    S_IBVIP.g_i=0.45* msiemens * cm **-2*int(version=='Alex')+(0.4*msiemens * cm **-2)*int(version=='Mark')+(0.02*msiemens * cm **-2)*int(version=='Mark/cell')
    S_IBVIP.g_i=0* msiemens * cm **-2
    S_IBVIP.taur_i=1.25*ms
    S_IBVIP.taud_i=50*ms
    S_IBVIP.V_i=0*mV
    
    
    
    ##Define gap junctions
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
    gap_IBIB.g_i=0.0025* msiemens * cm **-2

    #Define monitors :
    V5=StateMonitor(IB_soma,'V',record=True)
    R5=SpikeMonitor(IB_soma,record=True)
    I5s=StateMonitor(IB_soma,'Isyn',record=True)
    I5a=StateMonitor(IB_axon,'Isyn',record=True)
    I5ad=StateMonitor(IB_ad,'Isyn',record=True)
    I5bd=StateMonitor(IB_bd,'Isyn',record=True)

    new_neurons=IB_soma,IB_axon,IB_bd,IB_ad
    new_synapses=S_RSIB_AMPA,S_RSIB_NMDA, S_SIIB, S_IBFS, S_IBSI, S_IBIB, S_IBVIP
    new_synapses=tuple([y for y in new_synapses if y])
    new_gap_junctions=gapIB_SomaAd,gapIB_SomaBd,gapIB_SomaAxon,gapIB_AdSoma,gapIB_BdSoma,gapIB_AxonSoma,gap_IBIB
    new_monitors=V5,R5,I5s,I5a,I5ad,I5bd
    
    all_neurons=all_neurons+new_neurons
    all_synapses=all_synapses+new_synapses
    all_gap_junctions=all_gap_junctions+new_gap_junctions
    all_monitors=all_monitors+new_monitors
    
    return all_neurons, all_synapses, all_gap_junctions, all_monitors
    

#Other inputs (depends on simulation)
if __name__=='__main__' :
    
    start_scope()
    
    runtime=1*second
    
    NN=1 #multiplicative factor on the number of neurons
    N_RS,N_FS,N_SI,N_VIP,N_IB= NN*80,NN*20,NN*20,NN*20,NN*20 #Number of neurons of RE, TC, and HTC type
    
    net = Network()
    all_neurons, all_synapses, all_gap_junctions, all_monitors=create_Mark_Alex_network(kainate,version,Nf=NN)
    V1,V2,V3,V4,R1,R2,R3,R4,I1,I2,I3,I4,V5,R5,I5s,I5a,I5ad,I5bd=all_monitors

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
        
    
    prefs.codegen.target = 'cython'  #cython=faster, numpy = default python
    
    net.run(runtime,report='text',report_period=300*second)
    
    figure()
    plot(R1.t,R1.i+80,'r.',label='RS cells')
    plot(R2.t,R2.i+60,'b.',label='FS cells')
    plot(R3.t,R3.i+40,'g.',label='SI cells')
    plot(R4.t,R4.i+20,'k.',label='VIP cells')
    plot(R5.t,R5.i,'y.',label='IB cells')
    xlim(0,runtime/second)
    legend(loc='upper left')
    #
    min_t=int(50*ms*100000*Hz)
    LFP_V_RS=1/N_RS*sum(V1.V,axis=0)[min_t:]
    LFP_V_FS=1/N_FS*sum(V2.V,axis=0)[min_t:]
    LFP_V_SI=1/N_SI*sum(V3.V,axis=0)[min_t:]
    LFP_V_VIP=1/N_VIP*sum(V4.V,axis=0)[min_t:]
    LFP_V_IB=1/N_IB*sum(V5.V,axis=0)[min_t:]
    
    f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_FS=signal.periodogram(LFP_V_FS, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_SI=signal.periodogram(LFP_V_SI, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_VIP=signal.periodogram(LFP_V_VIP, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_IB=signal.periodogram(LFP_V_IB, 100000,'flattop', scaling='spectrum')
    
    figure()
    subplot(521)
    plot((V1.t/second)[min_t:],LFP_V_RS)
    ylabel('RS LFP')
    #title('RS cell')
    subplot(523)
    plot((V1.t/second)[min_t:],LFP_V_FS)
    ylabel('FS LFP')
    #title('FS cell')
    subplot(525)
    plot((V1.t/second)[min_t:],LFP_V_SI)
    ylabel('SI LFP')
    #title('SI cell')
    subplot(527)
    plot((V1.t/second)[min_t:],LFP_V_VIP)
    ylabel('VIP LFP')
    #title('VIP cell')
    subplot(529)
    plot((V1.t/second)[min_t:],LFP_V_IB)
    xlabel('Time (s)')
    ylabel('IB LFP')
    #title('IB cell')
    
    subplot(522)
    plot(f,Spectrum_LFP_V_RS)
    ylabel('RS Spectrum')
    yticks([],[])
    xlim(0,50)
    #title('RS cell')
    subplot(524)
    plot(f,Spectrum_LFP_V_FS)
    ylabel('FS Spectrum')
    yticks([],[])
    xlim(0,50)
    #title('FS cell')
    subplot(526)
    plot(f,Spectrum_LFP_V_SI)
    ylabel('SI Spectrum')
    yticks([],[])
    xlim(0,50)
    #title('SI cell')
    subplot(528)
    plot(f,Spectrum_LFP_V_VIP)
    ylabel('VIP Spectrum')
    yticks([],[])
    xlim(0,50)
    #title('VIP cell')
    subplot(5,2,10)
    plot(f,Spectrum_LFP_V_IB)
    xlabel('Frequency (Hz)')
    ylabel('IB Spectrum')
    yticks([],[])
    xlim(0,50)
    title('IB cell')
    
    clear_cache('cython')