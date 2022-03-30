# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:35:36 2019

@author: aaussel
"""

from IB_soma_LIP import *
from IB_axon_LIP import *
from IB_apical_dendrite_LIP import *
from IB_basal_dendrite_LIP import *

from scipy import signal

prefs.codegen.target = 'numpy'

defaultclock.dt = 0.01*ms
runtime=1*second

start_scope()

N_IB=100

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
IB_ad.J='25.5 * uA * cmeter ** -2'  #article SI=27.5, code=25.5 #Changed here to represent external excitation input to the apical dendrite

IB_bd=NeuronGroup(N_IB,eq_IB_bd,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
IB_bd.V = '-100*mvolt+10*rand()*mvolt'
IB_bd.h = '0+0.05*rand()'
IB_bd.m = '0+0.05*rand()'
IB_bd.mAR = '0+0.001*rand()'
IB_bd.mKM = '0+0.05*rand()'
IB_bd.mCaH = '0+0.01*rand()'
#IB_bd.J='42.5 * uA * cmeter ** -2' #article SI=44.5, code=42.5
IB_bd.J='42.5 * uA * cmeter ** -2-i/100 *80 * uA * cmeter ** -2'

##Synapses
eq_syn='''_post=s_i*g_i*(V_post-V_i) : amp * meter ** -2 (summed)
    ds_i/dt=-s_i/taud_i+(1-s_i)/taur_i*0.5*(1+tanh(V_pre/10/mV)) : 1
    g_i : siemens * meter**-2
    V_i : volt
    taud_i : second
    taur_i : second
'''

S_IBIB=Synapses(IB_axon,IB_bd,model='IsynIB_LIP'+eq_syn)
S_IBIB.connect()
S_IBIB.g_i=1/500* msiemens * cm **-2
S_IBIB.taur_i=0.25*ms
S_IBIB.taud_i=100*ms
S_IBIB.V_i=0*mV


##Gap junctions 

eq_gap='''_post=g_i*(V_post-V_pre) : amp * meter ** -2 (summed)
    g_i : siemens * meter**-2
'''

#gapIB_SomaAd=Synapses(IB_soma,IB_ad,model='Igap_soma'+eq_gap)
#gapIB_SomaAd.connect(j='i')
#gapIB_SomaAd.g_i=0.3* msiemens * cm **-2

gapIB_SomaBd=Synapses(IB_soma,IB_bd,model='Igap_soma'+eq_gap)
gapIB_SomaBd.connect(j='i')
gapIB_SomaBd.g_i=0.2* msiemens * cm **-2

gapIB_SomaAxon=Synapses(IB_soma,IB_axon,model='Igap_soma'+eq_gap)
gapIB_SomaAxon.connect(j='i')
gapIB_SomaAxon.g_i=0.3* msiemens * cm **-2

#gapIB_AdSoma=Synapses(IB_ad,IB_soma,model='Igap_ad'+eq_gap)
#gapIB_AdSoma.connect(j='i')
#gapIB_AdSoma.g_i=0.5* msiemens * cm **-2

gapIB_BdSoma=Synapses(IB_bd,IB_soma,model='Igap_bd'+eq_gap)
gapIB_BdSoma.connect(j='i')
gapIB_BdSoma.g_i=0.4* msiemens * cm **-2

gapIB_AxonSoma=Synapses(IB_axon,IB_soma,model='Igap_axon'+eq_gap)
gapIB_AxonSoma.connect(j='i')
gapIB_AxonSoma.g_i=0.3* msiemens * cm **-2

gap_IBIB=Synapses(IB_axon,IB_axon,model='Igap_axon'+eq_gap)
gap_IBIB.connect()
gap_IBIB.g_i=0.0025* msiemens * cm **-2


##Define monitors
#V1=StateMonitor(IB_soma,'V',record=True)
#V2=StateMonitor(IB_axon,'V',record=True)
#V3=StateMonitor(IB_ad,'V',record=True)
#V4=StateMonitor(IB_bd,'V',record=True)

R1=SpikeMonitor(IB_soma,record=True)
#R2=SpikeMonitor(IB_axon,record=True)
#R3=SpikeMonitor(IB_ad,record=True)
#R4=SpikeMonitor(IB_bd,record=True)

if __name__=='__main__':
    runtime=10*second
    f=20*Hz #rythmic input frequency
#    input_on = True
    input_on=False
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    
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
    
    if input_on:
        IB_bd.ginp_IB=-2*msiemens*cmeter**-2
        inputs_topdown=generate_spike_timing(N_IB,f,0*ms,end_time=1200*ms)
#        print(inputs_topdown)
        G_topdown = SpikeGeneratorGroup(N_IB, inputs_topdown[:,1], inputs_topdown[:,0]*second)
        topdown_in=Synapses(G_topdown,IB_bd,on_pre='Vinp=Vhigh')
        topdown_in.connect(j='i')
    
    
#    sig_ranIB_bd=0.005* mamp * cm **-2*0
#    sig_ranIB_ad=0.005* mamp * cm **-2*0
#    sig_ranIB_axon=0.005* mamp * cm **-2*0
#    sig_ranIB_soma=0.005* mamp * cm **-2*0
    run(runtime,report='text',report_period=300*second)
    
#    figure()
#    subplot(221)
#    plot(V1.t/second,V1.V[0]/volt)
#    xlabel('Time (s)')
#    ylabel('Membrane potential (V)')
#    title('Soma')
#    subplot(222)
#    plot(V1.t/second,V2.V[0]/volt)
#    xlabel('Time (s)')
#    ylabel('Membrane potential (V)')
#    title('Axon')
#    subplot(223)
#    plot(V1.t/second,V3.V[0]/volt)
#    xlabel('Time (s)')
#    ylabel('Membrane potential (V)')
#    title('Apical dendrite')
#    subplot(224)
#    plot(V1.t/second,V4.V[0]/volt)
#    xlabel('Time (s)')
#    ylabel('Membrane potential (V)')
#    title('Basal dendrite')
    
#    figure()
#    subplot(411)
#    plot(R1.t,R1.i,'r.')
#    xlim(0,runtime/second)
#    title('Soma')
#    subplot(412)
#    plot(R2.t,R2.i,'r.')
#    xlim(0,runtime/second)
#    title('Axon')
#    subplot(413)
#    plot(R3.t,R3.i,'r.')
#    xlim(0,runtime/second)
#    title('Apical dendrite')
#    subplot(414)
#    plot(R4.t,R4.i,'r.')
#    xlim(0,runtime/second)
#    title('Basal dendrite')
    
#    min_t=int(50*ms*100000*Hz)
#    LFP_V_IB=1/N_IB*sum(V1.V,axis=0)[min_t:]
#    
#    f,Spectrum_LFP_V_IB=signal.periodogram(LFP_V_IB, 100000,'flattop', scaling='spectrum')
#    
#    figure()
#    subplot(121)
#    plot((V1.t/second)[min_t:],LFP_V_IB)
#    title('IB cell')
#    subplot(122)
#    plot(f,Spectrum_LFP_V_IB)
#    xlabel('Frequency (Hz)')
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,50)
#    title('IB cell')
    
    
    figure()
    plot((-(IB_bd.J/ (uA * cmeter ** -2)-42.5)),R1.count/10)
    xlabel(r'I ($\mu A/cm^2$)', fontsize=14)
    ylabel('f (Hz)', fontsize=14)
    tick_params(axis='both', which='major', labelsize=12)