# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:35:36 2019

@author: aaussel
"""

from cells.RE_mdPul import *
from cells.TC_mdPul import *
#from cells.HTC_mdPul import *
from cells.HTC_buffer_mdPul_Destxhe_tests import *
from scipy import signal


prefs.codegen.target = 'numpy'


defaultclock.dt = 0.01*ms
runtime=2*second


start_scope()

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

def mon_delay(mon,i,t_delay):
    return mon.V[i][-t_delay]

def generate_syn_delay(source,target,syntype,connection_pattern,g_i,taur_i,taud_i,V_i):
    eq_syn_delay='''_post=s_i*g_i*(V_post-V_i) : amp * meter ** -2 (summed)
        ds_i/dt=-s_i/taud_i+(1-s_i)/taur_i*0.5*(1+tanh(voltage_delayed_pre/10/mV)) : 1
        g_i : siemens * meter**-2
        V_i : volt
        taud_i : second
        taur_i : second
    '''

    S=Synapses(source,target,model=syntype+eq_syn_delay,method='exact')
    if connection_pattern=='':
        S.connect()
    else :
        S.connect(j=connection_pattern, skip_if_invalid=True)
    S.g_i=g_i
    S.taur_i=taur_i
    S.taud_i=taud_i
    S.V_i=V_i  
    return S

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

def create_mdPul_column(N_HTC,N_TC,N_RE,condition,in_mode,theta_phase):
    HTC=NeuronGroup(N_HTC,eq_HTC_buffer_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    HTC.V = '-25*mvolt'
    HTC.Ca_iTLT = '1e-7 * mole * metre**-3'
    HTC.Ca_iTHT = '0.01*mmolar'
    HTC.J='0 * uA * cmeter ** -2'
    HTC.mAHP='0.3'
    
    TC=NeuronGroup(N_TC,eq_TC_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    TC.V = '-45*mvolt+20*mvolt*rand()'
    TC.J='0 * nA * cmeter ** -2'
    TC.Ca_i = '1e-7 * mole * metre**-3'

#    HTC.amptheta=3 * uA * cmeter ** -2
#    HTC.ftheta=4*Hz
#    
#    TC.amptheta=1 * uA * cmeter ** -2
#    TC.ftheta=4*Hz

    
    ##Define monitors
    R1=SpikeMonitor(HTC,record=True)
    R2=SpikeMonitor(TC,record=True)
    
    V1=StateMonitor(HTC,'V',record=True)
    V2=StateMonitor(TC,'V',record=True)
    
    I2=StateMonitor(HTC,'ITHT',record=[0])
    I3=StateMonitor(HTC,'ITLT',record=[0])
    I4=StateMonitor(HTC,'Itheta',record=[0])
    
    
    ##Synapses
        

    gGABA_A_HTC_TC = 0.2 * msiemens * cm **-2
    
    synHTCTC=generate_syn_delay(HTC,TC,'IsynHTC','',gGABA_A_HTC_TC,0.25*ms,5*ms,-80*mV)
#    synHTCTC=generate_syn(HTC,TC,'IsynHTC','',gGABA_A_HTC_TC,0.25*ms,5*ms,-80*mV)
#    if in_mode=='burst':
#        if condition=='mAChR':
#            synHTCTC=generate_syn_delay(HTC,TC,'IsynHTC','',gGABA_A_HTC_TC,0.25*ms,5*ms,-80*mV)
#    #        synHTCTC.delay=delayTC
#            
#        if condition=='mGluR1':
#            synHTCTC=generate_syn_delay(HTC,TC,'IsynHTC','',gGABA_A_HTC_TC,0.25*ms,5*ms,-80*mV)
#    #        synHTCTC.delay=delayTC
#    else :
#        if condition=='mAChR':
#            synHTCTC=generate_syn(HTC,TC,'IsynHTC','',gGABA_A_HTC_TC,0.25*ms,5*ms,-80*mV)
#    #        synHTCTC.delay=delayTC
#            
#        if condition=='mGluR1':
#            synHTCTC=generate_syn(HTC,TC,'IsynHTC','',gGABA_A_HTC_TC,0.25*ms,5*ms,-80*mV)
#    #        synHTCTC.delay=delayTC
        
    GJ_HTC=Synapses(HTC,HTC,'''
                 w : siemens * meter **-2 # gap junction conductance
                 IGJ_post = w * (V_post - V_pre) : amp * meter ** -2 (summed)
                 ''')
    GJ_HTC.connect()
    #GJ_HTC.w = 0.003e-3 * siemens * cm **-2
    GJ_HTC.w = 0.08e-3 * siemens * cm **-2
    
#    gEPSP_RE = 20 * msiemens * cm **-2 #0.42 - 0.70
##    gEPSP_TC = 0.5 * msiemens * cm **-2 #0.42 - 0.70
#    gEPSP_TC = 0 * msiemens * cm **-2 #0.42 - 0.70
#    gEPSP_TC2 = 10 * msiemens * cm **-2 #0.42 - 0.70
#    gIPSP_TC = -1.5 * gEPSP_TC    
#    #Inputs : 
#    inputs_TC=generate_spike_timing(N_TC,100*Hz,0*ms,end_time=3000*ms)
#    G_in = SpikeGeneratorGroup(N_TC, inputs_TC[:,1], inputs_TC[:,0]*second)
#    Syn_in=Synapses(G_in,TC,on_pre='Vinp=Vhigh')
#    Syn_in.connect(j='i')
#    TC.ginp_TC=gEPSP_TC
#    
#    inputs_TC3=generate_spike_timing(N_TC,13*Hz,0*ms,end_time=3000*ms)
#    G_in3 = SpikeGeneratorGroup(N_TC, inputs_TC3[:,1], inputs_TC3[:,0]*second)
#    Syn_in3=Synapses(G_in3,TC,on_pre='Vinp=Vhigh')
#    Syn_in3.connect(j='i')
#    
##    inputs_HTC=generate_spike_timing(N_HTC,13*Hz,0*ms,end_time=3000*ms)
##    G_in4 = SpikeGeneratorGroup(N_HTC, inputs_HTC[:,1], inputs_HTC[:,0]*second)
##    Syn_in4=Synapses(G_in4,HTC,on_pre='Vinp=Vhigh')
##    Syn_in4.connect(j='i')  
##    HTC.ginp_HTC=0e-3 * siemens * cm **-2
#    
#    if theta_phase=='bad':
#        TC.ginp_TC2=gEPSP_TC2
    
    all_neurons=HTC,TC
    all_synapses=(synHTCTC,)
    all_synapses=tuple([y for y in all_synapses if y])
    all_monitors=R1,R2,V1,V2,I2,I3,I4
    all_gap_junctions=(GJ_HTC,)
    
#    if in_mode=='single_spike':
#        HTC.delay_steps = [1]  # delay in time steps per neuron
#        buffer_size = 2  # 1+Maximum delay (in time steps)
#    else :
#        HTC.delay_steps = [3999]  # delay in time steps per neuron
#        buffer_size = 4000  # 1+Maximum delay (in time steps)
#        
#    print(HTC.delay_steps)
#    HTC.variables.add_array('voltage_buffer', dimensions=volt.dim, size=(buffer_size, len(HTC)))
#    
#    update_code = '''buffer_pointer = (buffer_pointer + 1) % buffer_size
#                     voltage_delayed = update_voltage_buffer(V, voltage_buffer, buffer_pointer, delay_steps, buffer_size)'''
#       
#    buffer_updater = HTC.run_regularly(update_code, codeobj_class=NumpyCodeObject)
#        
#    @check_units(V=volt, voltage_buffer=volt, buffer_pointer=1, delay_steps=1, buffer_size=1, result=volt)
#    def update_voltage_buffer(V, voltage_buffer, buffer_pointer, delay_steps, buffer_size):
#        # Write current rate into the buffer
#        voltage_buffer[buffer_pointer, :] = V
#        # Get delayed rates 
#        rows = (buffer_pointer - delay_steps) % buffer_size    
#        return voltage_buffer[rows, arange(len(rows))]
    
    return all_neurons,all_synapses,all_gap_junctions,all_monitors


if __name__=='__main__':
    close('all')
    runtime=1*second
    f=13*Hz #rythmic input frequency
    input_on=False
    N_HTC,N_TC,N_RE= 20,80,100 #Number of neurons of RE, TC, and HTC type
#    N_HTC,N_TC,N_RE= 1,80,100 #Number of neurons of RE, TC, and HTC type

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
        
#    Cm_HTC = 2.5* ufarad * cm ** -2
#    gNa_HTC=90e-3 * siemens * cm **-2
#    ENa_HTC=50*mV
#    gK_HTC=10e-3 * siemens * cm **-2
#    EK_HTC=-100*mV
#    gL_HTC=0.001e-3 * siemens * cm **-2  
#    EL_HTC=-70*mV
    gKL_HTC=0.001e-3 * siemens * cm **-2  
#    EKL_HTC=-100*mV
#    gTLT_HTC= 2.1e-3 * siemens * cm **-2
#    gTHT_HTC= 15e-3 * siemens * cm **-2
#    gAHP_HTC= 45e-3 * siemens * cm **-2
#    EAHP=-95*mV
#    gH_HTC = 0.36e-3 * siemens * cm **-2
#    EH_HTC = -40 * mV
    gapp=0.1*mamp * cmeter ** -2 # in HTC cells
#    gapp=0.001*mamp * cmeter ** -2 # in HTC cells
    
    net=Network()
    all_neurons,all_synapses,all_gap_junctions,all_monitors=create_mdPul_column(N_HTC,N_TC,N_RE,condition,in_mode,theta_phase)
#    R1,R2,R3,V1,V2,V3=all_monitors
    R1,R2,V1,V2,I2,I3,I4=all_monitors
    
    HTC=all_neurons[0]
    if in_mode=='single_spike':
        HTC.delay_steps = [1]  # delay in time steps per neuron
        buffer_size = 2  # 1+Maximum delay (in time steps)
    else :
        HTC.delay_steps = [3999]  # delay in time steps per neuron
        buffer_size = 4000  # 1+Maximum delay (in time steps)
    
    init_array=-70*mV*ones((buffer_size, len(HTC)))
    HTC.variables.add_array('voltage_buffer', dimensions=volt.dim, size=(buffer_size, len(HTC)))#,values=init_array)
    HTC.voltage_buffer=init_array
    
    update_code = '''buffer_pointer = (buffer_pointer + 1) % buffer_size
                     voltage_delayed = update_voltage_buffer(V, voltage_buffer, buffer_pointer, delay_steps, buffer_size)'''
       
    buffer_updater = HTC.run_regularly(update_code, codeobj_class=NumpyCodeObject)
        
    @check_units(V=volt, voltage_buffer=volt, buffer_pointer=1, delay_steps=1, buffer_size=1, result=volt)
    def update_voltage_buffer(V, voltage_buffer, buffer_pointer, delay_steps, buffer_size):
        # Write current rate into the buffer
        voltage_buffer[buffer_pointer, :] = V
        # Get delayed rates 
        rows = (buffer_pointer - delay_steps) % buffer_size    
        return voltage_buffer[rows, arange(len(rows))]
        
#    print(HTC.delay_steps)
        
#    SynHTCTC=all_synapses[5]
    
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_gap_junctions)
    net.add(all_monitors)

#    TC=all_neurons[1]
#    M=StateMonitor(TC,'Isyn', record=True)
#    M2=StateMonitor(HTC,'voltage_delayed', record=True)
#    M3=StateMonitor(SynHTCTC,'s_i',record=[0])
#    net.add(M,M2,M3)
            
    prefs.codegen.target = 'cython'
    net.run(runtime,report='text',report_period=300*second)
    
    
    figure()
    plot(R1.t,R1.i+0,'r.',label='HTC')
    plot(R2.t,R2.i+20,'b.',label='TC')
    xlim(0,runtime/second)
    legend()  
    
    figure()
    plot(V1.t,V1.V[0],label='HTC V')
    plot(V2.t,V2.V[0],label='TC V')
    legend()
    
#    figure()
#    plot(I1.t,I1.ISK[0],label='I_SK HTC')
#    legend()
    
#    figure()
#    plot(I2.t,I2.ITHT[0],label='I_THT HTC')
#    legend()
#    
#    figure()
#    plot(I3.t,I3.ITLT[0],label='I_TLT HTC')
#    legend()
    
    figure()
    plot(I4.t,I4.Itheta[0],label='Itheta HTC')
    legend()
        
#    f,Spectrum_LFP_V1=signal.periodogram(V1.V[0], 100000,'flattop', scaling='spectrum')
#    figure()
#    plot(f,Spectrum_LFP_V1)
#    xlim(0,100)
#    Isyn_calc=(V2.V[0]+80*mV)*0.2 * msiemens * cm **-2*M3.s_i[0]
#    
#    figure()
#    plot(M.t,M.Isyn[0],label='I_syn HTC-TC')
#    plot(M2.t,M2.voltage_delayed[0],label='HTC voltage_delayed')
#    plot(V2.t,V2.V[0],label='TC V')
#    plot(V1.t,V1.V[0],label='HTC V')
#    plot(M3.t,M3.s_i[0],label='s_i syn HTC-TC')
#    plot(M3.t,Isyn_calc,label='Isyn HTC-TC calc')
#    legend()
    

    clear_cache('cython') 
    
