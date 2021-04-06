# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:35:36 2019

@author: aaussel
"""

from cells.RE_mdPul import *
from cells.TC_mdPul import *
from cells.HTC_buffer_mdPul_Destxhe_tests import *
from mdPul_one_column import *

from scipy import signal


prefs.codegen.target = 'numpy'


defaultclock.dt = 0.01*ms
runtime=1*second

start_scope()

def create_mdPul(N_HTC,N_TC,N_RE,condition,in_mode,theta_phase):
    
    all_neuronsA,all_synapsesA,all_gap_junctionsA,all_monitorsA=create_mdPul_column(N_HTC,N_TC,N_RE,condition,in_mode,theta_phase)
    all_neuronsB,all_synapsesB,all_gap_junctionsB,all_monitorsB=create_mdPul_column(N_HTC,N_TC,N_RE,condition,in_mode,theta_phase)
    all_neuronsC,all_synapsesC,all_gap_junctionsC,all_monitorsC=create_mdPul_column(N_HTC,N_TC,N_RE,condition,in_mode,theta_phase)
       
    HTC_A=all_neuronsA[0]
    HTC_B=all_neuronsB[0]
    HTC_C=all_neuronsC[0]
    
    RE_A=NeuronGroup(N_RE,eq_RE_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    RE_A.V = '-70*mvolt+20*mvolt*rand()'
#    RE_A.J = '600 * nA * cmeter ** -2'
    RE_A.J = '0 * nA * cmeter ** -2'
    RE_A.Ca_i = '1e-7 * mole * metre**-3'
    
    RE_B=NeuronGroup(N_RE,eq_RE_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    RE_B.V = '-70*mvolt+20*mvolt*rand()'
#    RE_B.J = '600 * nA * cmeter ** -2'
    RE_B.J = '0 * nA * cmeter ** -2'
    RE_B.Ca_i = '1e-7 * mole * metre**-3'
    
    RE_C=NeuronGroup(N_RE,eq_RE_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    RE_C.V = '-70*mvolt+20*mvolt*rand()'
#    RE_C.J = '600 * nA * cmeter ** -2'
    RE_C.J = '0 * nA * cmeter ** -2'
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
#    g_HTCRE=0.6 * msiemens * cm **-2 #0.4
#    g_REHTC=0.6 * msiemens * cm **-2
#    g_AB=0.9
#    g_AC=1.1
    
    g_HTCRE=0.2 * msiemens * cm **-2 #0.4
    g_REHTC=0.2 * msiemens * cm **-2
    g_RERE=0.05 * msiemens * cm **-2
    g_AB=1
    g_AC=1
    
#    HTC_A.J='1 * nA * cmeter ** -2'
#    HTC_B.J='5 * nA * cmeter ** -2'
#    HTC_C.J='5 * nA * cmeter ** -2'
    
    HTC_A.J='5*(1-0.4*t/second) * nA * cmeter ** -2'
    HTC_B.J='15*(1-0.4*t/second) * nA * cmeter ** -2'
    HTC_C.J='15*(1-0.4*t/second) * nA * cmeter ** -2'

    HTC_A.run_regularly('J = 5*(1-0.4*t/second) * nA * cmeter ** -2', dt=0.1*ms)
    HTC_B.run_regularly('J = 15*(1-0.4*t/second) * nA * cmeter ** -2', dt=0.1*ms)
    HTC_C.run_regularly('J = 15*(1-0.4*t/second) * nA * cmeter ** -2', dt=0.1*ms)
    
    
#    g_HTCRE=0.5 * msiemens * cm **-2 #0.4
#    g_REHTC=0.5 * msiemens * cm **-2
#    g_AB=0.85
#    g_AC=1.15
    
    
    S_HTC_RE_A=Synapses(HTC_A,RE_A,model='IsynHTC'+eq_syn)
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
    
    S_RE_A_RE_A=Synapses(RE_A,RE_A,model='IsynREA'+eq_syn)
    S_RE_A_RE_A.connect()
    S_RE_A_RE_A.g_i=g_RERE
    S_RE_A_RE_A.taur_i=0.25*ms
    S_RE_A_RE_A.taud_i=20*ms
    S_RE_A_RE_A.V_i=-80*mV
    
    
    S_HTC_RE_B=Synapses(HTC_B,RE_B,model='IsynHTC'+eq_syn)
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
    
    S_RE_B_RE_B=Synapses(RE_B,RE_B,model='IsynREB'+eq_syn)
    S_RE_B_RE_B.connect()
    S_RE_B_RE_B.g_i=g_RERE
    S_RE_B_RE_B.taur_i=0.25*ms
    S_RE_B_RE_B.taud_i=20*ms
    S_RE_B_RE_B.V_i=-80*mV
    
    
    S_HTC_RE_C=Synapses(HTC_C,RE_C,model='IsynHTC'+eq_syn)
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
    
    S_RE_C_RE_C=Synapses(RE_C,RE_C,model='IsynREA'+eq_syn)
    S_RE_C_RE_C.connect()
    S_RE_C_RE_C.g_i=g_RERE
    S_RE_C_RE_C.taur_i=0.25*ms
    S_RE_C_RE_C.taud_i=20*ms
    S_RE_C_RE_C.V_i=-80*mV
    
    
    ##Define monitors
    RA=SpikeMonitor(RE_A,record=True)
    
    RB=SpikeMonitor(RE_B,record=True)
    
    RC=SpikeMonitor(RE_C,record=True)
    
    
    all_neurons=all_neuronsA+all_neuronsB+all_neuronsC+(RE_A,RE_B,RE_C,)
    all_synapses=all_synapsesA,all_synapsesB+all_synapsesC+(S_HTC_RE_A,S_RE_A_HTC_B,S_RE_A_HTC_C,S_HTC_RE_B,S_RE_B_HTC_A,S_RE_B_HTC_C,S_HTC_RE_C,S_RE_C_HTC_A,S_RE_C_HTC_B,S_RE_A_RE_A,S_RE_B_RE_B,S_RE_C_RE_C)
    all_monitors=all_monitorsA+all_monitorsB+all_monitorsC+(RA,RB,RC)
    all_gap_junctions=all_gap_junctionsA,all_gap_junctionsB,all_gap_junctionsC
    
    return all_neurons,all_synapses,all_gap_junctions,all_monitors


if __name__=='__main__':
    close('all')
    runtime=1*second
    f=13*Hz #rythmic input frequency
    input_on=False
    N_HTC, N_TC,N_RE= 20,80,100 #Number of neurons of RE, TC, and HTC type

    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp_IB=0* msiemens * cm **-2
    ginp=0* msiemens * cm **-2
    
    
    Vrev_inp2=0*mV
    taurinp2=0.1*ms
    taudinp2=0.5*ms
    tauinp2=taudinp2
    Vhigh2=0*mV
    Vlow2=-80*mV
    
    #condition='mGluR1'
    condition='mAChR'
    in_mode='single_spike'
#    in_mode='burst'
    theta_phase='good'

    if condition=='mGluR1':
        gKL_TC=0.0028e-3 * siemens * cm **-2
        gKL_HTC=0.0069e-3 * siemens * cm **-2
        gKL_RE=0.05e-3 * siemens * cm **-2   
    elif condition=='mAChR':
        gKL_TC=0.0028e-3 * siemens * cm **-2
        gKL_HTC=0.0069e-3 * siemens * cm **-2
        gKL_RE=0.08e-3 * siemens * cm **-2
        

#    gapp=0.1*mamp * cmeter ** -2 # in HTC cells
    gKL_HTC=0.001e-3 * siemens * cm **-2 
    gapp=0.1*mamp * cmeter ** -2 # in HTC cells    
    
    net=Network()
    all_neurons,all_synapses,all_gap_junctions,all_monitors=create_mdPul(N_HTC,N_TC,N_RE,condition,in_mode,theta_phase)
    R1A,R2A,R3A,V1A,V2A,V3A,I1A,I2A,R1B,R2B,R3B,V1B,V2B,V3B,I1B,I2B,R1C,R2C,R3C,V1C,V2C,V3C,I1C,I2C,RA,RB,RC=all_monitors
    

    
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
    
    HTC_A,HTC_B,HTC_C=all_neurons[0],all_neurons[6],all_neurons[12]
    if in_mode=='single_spike':
        HTC_A.delay_steps = [1]  # delay in time steps per neuron
        HTC_B.delay_steps = [1]  # delay in time steps per neuron
        HTC_C.delay_steps = [1]  # delay in time steps per neuron
        buffer_size = 2  # 1+Maximum delay (in time steps)
    else :
        HTC_A.delay_steps = [3999]  # delay in time steps per neuron
        HTC_B.delay_steps = [3999]  # delay in time steps per neuron
        HTC_C.delay_steps = [3999]  # delay in time steps per neuron
        buffer_size = 4000  # 1+Maximum delay (in time steps)
        
    HTC_A.variables.add_array('voltage_buffer', dimensions=volt.dim, size=(buffer_size, len(HTC_A)))
    HTC_B.variables.add_array('voltage_buffer', dimensions=volt.dim, size=(buffer_size, len(HTC_B)))
    HTC_C.variables.add_array('voltage_buffer', dimensions=volt.dim, size=(buffer_size, len(HTC_C)))

    init_array=-70*mV*ones((buffer_size, len(HTC_A)))
    HTC_A.voltage_buffer=init_array
    HTC_B.voltage_buffer=init_array
    HTC_C.voltage_buffer=init_array
    
    update_code = '''buffer_pointer = (buffer_pointer + 1) % buffer_size
                     voltage_delayed = update_voltage_buffer(V, voltage_buffer, buffer_pointer, delay_steps, buffer_size)'''
       
    buffer_updater_A = HTC_A.run_regularly(update_code, codeobj_class=NumpyCodeObject)
    buffer_updater_B = HTC_B.run_regularly(update_code, codeobj_class=NumpyCodeObject)
    buffer_updater_C = HTC_C.run_regularly(update_code, codeobj_class=NumpyCodeObject)
        
    @check_units(V=volt, voltage_buffer=volt, buffer_pointer=1, delay_steps=1, buffer_size=1, result=volt)
    def update_voltage_buffer(V, voltage_buffer, buffer_pointer, delay_steps, buffer_size):
        # Write current rate into the buffer
        voltage_buffer[buffer_pointer, :] = V
        # Get delayed rates 
        rows = (buffer_pointer - delay_steps) % buffer_size    
        return voltage_buffer[rows, arange(len(rows))]

    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_gap_junctions)
    net.add(all_monitors)   

    MHTCA=StateMonitor(HTC_A,['J'],record=[0])  
    net.add(MHTCA)
    MHTCB=StateMonitor(HTC_B,['J'],record=[0])  
    net.add(MHTCB)
            
    prefs.codegen.target = 'cython'
    net.run(runtime,report='text',report_period=300*second)
    
    
    figure()
    plot(R1A.t,R1A.i+0,'r.',label='HTC')
    plot(R2A.t,R2A.i+20,'y.',label='TC')
    plot(R3A.t,R3A.i+100,'g.',label='RE int')
    plot(RA.t,RA.i+200,'b.',label='RE lat')

    plot(R1B.t,R1B.i+350,'r.')
    plot(R2B.t,R2B.i+370,'y.')
    plot(R3B.t,R3B.i+450,'g.')
    plot(RB.t,RB.i+550,'b.')

    plot(R1C.t,R1C.i+700,'r.')
    plot(R2C.t,R2C.i+720,'y.')
    plot(R3C.t,R3C.i+800,'g.')
    plot(RC.t,RC.i+900,'b.')
    
    xlim(0,runtime/second)  
    yticks([150,500,850],['A','B','C'])
    legend()  
    
    figure()
    plot(MHTCA.t,MHTCA.J[0],label='HTC_A.J')
    plot(MHTCB.t,MHTCB.J[0],label='HTC_B.J')
    legend()
    
#    figure()
#    plot(V1A.t,V1A.V[0],label='TC V')
#    plot(V2A.t,V2A.V[0],label='RE V')
#    legend()
    
#    f,Spectrum_LFP_V1=signal.periodogram(V1A.V[0], 100000,'flattop', scaling='spectrum')
#    figure()
#    plot(f,Spectrum_LFP_V1)
#    xlim(0,100)
    clear_cache('cython') 
    
