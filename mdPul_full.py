# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:35:36 2019

@author: aaussel
"""

from cells.RE_mdPul import *
from cells.TC_mdPul import *
from scipy import signal


prefs.codegen.target = 'numpy'


defaultclock.dt = 0.01*ms
runtime=3*second

start_scope()

def create_mdPul(N_TC,N_RE):
    HTC_A=NeuronGroup(N_TC,eq_TC_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    HTC_A.V = '-75*mvolt+20*mvolt*rand()'
    #HTC_A.Ca_iTLT = '1e-7 * mole * metre**-3'
    #HTC_A.Ca_iTHT = '1e-7 * mole * metre**-3'
    #HTC_A.J='-2 * uA * cmeter ** -2'
    HTC_A.J='-1 * nA * cmeter ** -2'
    HTC_A.Ca_i = '1e-7 * mole * metre**-3'
    
    HTC_B=NeuronGroup(N_TC,eq_TC_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    HTC_B.V = '-75*mvolt+20*mvolt*rand()'
    #HTC_B.Ca_iTLT = '1e-7 * mole * metre**-3'
    #HTC_B.Ca_iTHT = '1e-7 * mole * metre**-3'
    #HTC_B.J='-2 * uA * cmeter ** -2'
    HTC_B.J='-1 * nA * cmeter ** -2'
    HTC_B.Ca_i = '1e-7 * mole * metre**-3'
    
    HTC_C=NeuronGroup(N_TC,eq_TC_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    HTC_C.V = '-75*mvolt+20*mvolt*rand()'
    #HTC_C.Ca_iTLT = '1e-7 * mole * metre**-3'
    #HTC_C.Ca_iTHT = '1e-7 * mole * metre**-3'
    #HTC_C.J='-2 * uA * cmeter ** -2'
    HTC_C.J='-1 * nA * cmeter ** -2'
    HTC_C.Ca_i = '1e-7 * mole * metre**-3'
        
    
    RE_A=NeuronGroup(N_RE,eq_RE_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    RE_A.V = '-70*mvolt+20*mvolt*rand()'
    RE_A.J = '600 * nA * cmeter ** -2'
    RE_A.Ca_i = '1e-7 * mole * metre**-3'
    
    RE_B=NeuronGroup(N_RE,eq_RE_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    RE_B.V = '-70*mvolt+20*mvolt*rand()'
    RE_B.J = '600 * nA * cmeter ** -2'
    RE_B.Ca_i = '1e-7 * mole * metre**-3'
    
    RE_C=NeuronGroup(N_RE,eq_RE_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    RE_C.V = '-70*mvolt+20*mvolt*rand()'
    RE_C.J = '600 * nA * cmeter ** -2'
    RE_C.Ca_i = '1e-7 * mole * metre**-3'
    
    
    
    ##Synapses
    eq_syn='''_post=s_i*g_i*(V_post-V_i) : amp * meter ** -2 (summed)
        ds_i/dt=-s_i/taud_i+(1-s_i)/taur_i*0.5*(1+tanh(V_pre/10/mV)) : 1
        g_i : siemens * meter**-2
        V_i : volt
        taud_i : second
        taur_i : second
    '''
    #0.006, 0.03
    g_HTCRE=0.6 * msiemens * cm **-2 #0.4
    g_REHTC=0.6 * msiemens * cm **-2
    g_AB=0.9
    g_AC=1.1
    
    
    S_HTC_RE_A=Synapses(HTC_A,RE_A,model='IsynTC'+eq_syn)
    S_HTC_RE_A.connect()
    S_HTC_RE_A.g_i=g_HTCRE
    S_HTC_RE_A.taur_i=0.25*ms
    S_HTC_RE_A.taud_i=5*ms
    S_HTC_RE_A.V_i=0*mV
    
    S_RE_A_HTC_B=Synapses(RE_A,HTC_B,model='IsynREA'+eq_syn)
    S_RE_A_HTC_B.connect()
    S_RE_A_HTC_B.g_i=g_REHTC*g_AB
    S_RE_A_HTC_B.taur_i=0.25*ms
    S_RE_A_HTC_B.taud_i=20*ms
    S_RE_A_HTC_B.V_i=-80*mV
    
    S_RE_A_HTC_C=Synapses(RE_A,HTC_C,model='IsynREA'+eq_syn)
    S_RE_A_HTC_C.connect()
    S_RE_A_HTC_C.g_i=g_REHTC*g_AC
    S_RE_A_HTC_C.taur_i=0.25*ms
    S_RE_A_HTC_C.taud_i=20*ms
    S_RE_A_HTC_C.V_i=-80*mV
    
    
    S_HTC_RE_B=Synapses(HTC_B,RE_B,model='IsynTC'+eq_syn)
    S_HTC_RE_B.connect()
    S_HTC_RE_B.g_i=g_HTCRE
    S_HTC_RE_B.taur_i=0.25*ms
    S_HTC_RE_B.taud_i=5*ms
    S_HTC_RE_B.V_i=0*mV
    
    S_RE_B_HTC_A=Synapses(RE_B,HTC_A,model='IsynREB'+eq_syn)
    S_RE_B_HTC_A.connect()
    S_RE_B_HTC_A.g_i=g_REHTC*g_AB
    S_RE_B_HTC_A.taur_i=0.25*ms
    S_RE_B_HTC_A.taud_i=20*ms
    S_RE_B_HTC_A.V_i=-80*mV
    
    S_RE_B_HTC_C=Synapses(RE_B,HTC_C,model='IsynREB'+eq_syn)
    S_RE_B_HTC_C.connect()
    S_RE_B_HTC_C.g_i=g_REHTC
    S_RE_B_HTC_C.taur_i=0.25*ms
    S_RE_B_HTC_C.taud_i=20*ms
    S_RE_B_HTC_C.V_i=-80*mV
    
    
    S_HTC_RE_C=Synapses(HTC_C,RE_C,model='IsynTC'+eq_syn)
    S_HTC_RE_C.connect()
    S_HTC_RE_C.g_i=g_HTCRE
    S_HTC_RE_C.taur_i=0.25*ms
    S_HTC_RE_C.taud_i=5*ms
    S_HTC_RE_C.V_i=0*mV
    
    S_RE_C_HTC_A=Synapses(RE_C,HTC_A,model='IsynREC'+eq_syn)
    S_RE_C_HTC_A.connect()
    S_RE_C_HTC_A.g_i=g_REHTC*g_AC
    S_RE_C_HTC_A.taur_i=0.25*ms
    S_RE_C_HTC_A.taud_i=20*ms
    S_RE_C_HTC_A.V_i=-80*mV
    
    S_RE_C_HTC_B=Synapses(RE_C,HTC_B,model='IsynREC'+eq_syn)
    S_RE_C_HTC_B.connect()
    S_RE_C_HTC_B.g_i=g_REHTC
    S_RE_C_HTC_B.taur_i=0.25*ms
    S_RE_C_HTC_B.taud_i=20*ms
    S_RE_C_HTC_B.V_i=-80*mV
    
    
    ##Define monitors
    R1A=SpikeMonitor(HTC_A,record=True)
    R2A=SpikeMonitor(RE_A,record=True)
    
    R1B=SpikeMonitor(HTC_B,record=True)
    R2B=SpikeMonitor(RE_B,record=True)
    
    R1C=SpikeMonitor(HTC_C,record=True)
    R2C=SpikeMonitor(RE_C,record=True)
    
    V1A=StateMonitor(HTC_A,'V',record=[0])
    V2A=StateMonitor(RE_A,'V',record=[0])
    
    all_neurons=HTC_A,HTC_B,HTC_C,RE_A,RE_B,RE_C
    all_synapses=S_HTC_RE_A,S_RE_A_HTC_B,S_RE_A_HTC_C,S_HTC_RE_B,S_RE_B_HTC_A,S_RE_B_HTC_C,S_HTC_RE_C,S_RE_C_HTC_A,S_RE_C_HTC_B
    all_monitors=R1A,R2A,R1B,R2B,R1C,R2C,V1A,V2A
    
    return all_neurons,all_synapses,all_monitors


if __name__=='__main__':
    runtime=1*second
    f=13*Hz #rythmic input frequency
    input_on=False
    N_TC,N_RE= 20,20 #Number of neurons of RE, TC, and HTC type

    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp_IB=0* msiemens * cm **-2
    ginp=0* msiemens * cm **-2
    
    net=Network()
    all_neurons,all_synapses,all_monitors=create_mdPul(N_TC,N_RE)
    R1A,R2A,R1B,R2B,R1C,R2C,V1A,V2A=all_monitors
    
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_monitors)
    
#    def generate_spike_timing(N,f,start_time,end_time=runtime):
#        list_time_and_i=[]
#        for i in range(N):
#            list_time=[(start_time,i)]
#            next_spike=list_time[-1][0]+(1+0.1*rand())/f
#            while next_spike<end_time:
#                list_time.append((next_spike,i))
#                next_spike=list_time[-1][0]+(1+0.1*rand())/f
#            list_time_and_i+=list_time
#        return array(list_time_and_i)
#    
#    if input_on:
#        RSA.ginp_RS=1* msiemens * cm **-2 #1
#        inputs_topdown=generate_spike_timing(N_RS,f,0*ms,end_time=3000*ms)
##        print(inputs_topdown)
#        G_topdown = SpikeGeneratorGroup(N_RS, inputs_topdown[:,1], inputs_topdown[:,0]*second)
#        topdown_in=Synapses(G_topdown,RSA,on_pre='Vinp=Vhigh')
#        topdown_in.connect(j='i')
#        
#        RSB.ginp_RS=2* msiemens * cm **-2
#        inputs_topdown2=generate_spike_timing(N_RS,f,0*ms,end_time=3000*ms)
##        print(inputs_topdown2)
#        G_topdown2 = SpikeGeneratorGroup(N_RS, inputs_topdown2[:,1], inputs_topdown2[:,0]*second)
#        topdown_in2=Synapses(G_topdown2,RSB,on_pre='Vinp=Vhigh')
#        topdown_in2.connect(j='i')
#        
#        RSC.ginp_RS=1* msiemens * cm **-2
#        inputs_topdown3=generate_spike_timing(N_RS,f,0*ms,end_time=3000*ms)
##        print(inputs_topdown)
#        G_topdown3 = SpikeGeneratorGroup(N_RS, inputs_topdown3[:,1], inputs_topdown3[:,0]*second)
#        topdown_in3=Synapses(G_topdown3,RSC,on_pre='Vinp=Vhigh')
#        topdown_in3.connect(j='i')
#        
#        for elem in [G_topdown,topdown_in,G_topdown2,topdown_in2,G_topdown3,topdown_in3]:
#            net.add(elem)
            
    prefs.codegen.target = 'cython'
    net.run(runtime,report='text',report_period=300*second)
    
    
    figure()
    plot(R1A.t,R1A.i+0,'r.',label='TC')
    plot(R2A.t,R2A.i+20,'b.',label='RE')

    plot(R1B.t,R1B.i+50,'r.')
    plot(R2B.t,R2B.i+70,'b.')

    plot(R1C.t,R1C.i+100,'r.')
    plot(R2C.t,R2C.i+120,'b.')
    
    xlim(0,runtime/second)  
    yticks([20,70,120],['A','B','C'])
    legend()  
    
    figure()
    plot(V1A.t,V1A.V[0],label='TC V')
    plot(V2A.t,V2A.V[0],label='RE V')
    legend()
    
    f,Spectrum_LFP_V1=signal.periodogram(V1A.V[0], 100000,'flattop', scaling='spectrum')
    figure()
    plot(f,Spectrum_LFP_V1)
    xlim(0,100)
    clear_cache('cython') 
    
