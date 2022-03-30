# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 11:35:36 2019

@author: aaussel
"""

from cells.RE_mdPul import *
from cells.TC_mdPul import *
from cells.HTC_buffer_mdPul_Destxhe_tests import *
from mdPul_one_column_v3 import *
from cells.RS_FEF import *

from scipy import signal


prefs.codegen.target = 'numpy'


defaultclock.dt = 0.01*ms
runtime=1*second

start_scope()

def create_mdPul_2col(N_HTC,N_TC,N_RE,condition,in_mode,theta_phase):
    
    all_neuronsA,all_synapsesA,all_gap_junctionsA,all_monitorsA=create_mdPul_column(N_HTC,N_TC,N_RE,condition,in_mode,theta_phase)
    all_neuronsB,all_synapsesB,all_gap_junctionsB,all_monitorsB=create_mdPul_column(N_HTC,N_TC,N_RE,condition,in_mode,theta_phase)
       
    HTC_A,TC_A=all_neuronsA[0],all_neuronsA[1]
    HTC_B,TC_B=all_neuronsB[0],all_neuronsB[1]
    
    RE_A=NeuronGroup(N_RE,eq_RE_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    RE_A.V = '-70*mvolt+20*mvolt*rand()'
    RE_A.J = '0 * nA * cmeter ** -2'
    RE_A.Ca_i = '1e-7 * mole * metre**-3'
    
    RE_B=NeuronGroup(N_RE,eq_RE_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    RE_B.V = '-70*mvolt+20*mvolt*rand()'
    RE_B.J = '0 * nA * cmeter ** -2'
    RE_B.Ca_i = '1e-7 * mole * metre**-3'
    
    RE_lat_A=NeuronGroup(N_RE,eq_RE_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    RE_lat_A.V = '-70*mvolt+20*mvolt*rand()'
    RE_lat_A.J = '0 * nA * cmeter ** -2'
    RE_lat_A.Ca_i = '1e-7 * mole * metre**-3'
    
    RE_lat_B=NeuronGroup(N_RE,eq_RE_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    RE_lat_B.V = '-70*mvolt+20*mvolt*rand()'
    RE_lat_B.J = '0 * nA * cmeter ** -2'
    RE_lat_B.Ca_i = '1e-7 * mole * metre**-3'

    E_A=NeuronGroup(20,eq_RS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    E_A.V = '-70*mvolt+10*rand()*mvolt'
    E_A.h = '0+0.05*rand()'
    E_A.m = '0+0.05*rand()'
    E_A.mAR = '0.035+0.025*rand()'
    E_A.J='50 * uA * cmeter ** -2'  #article SI=25, code=1

    E_B=NeuronGroup(20,eq_RS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    E_B.V = '-70*mvolt+10*rand()*mvolt'
    E_B.h = '0+0.05*rand()'
    E_B.m = '0+0.05*rand()'
    E_B.mAR = '0.035+0.025*rand()'
    E_B.J='50 * uA * cmeter ** -2'  #article SI=25, code=1

    
    ##Synapses
    eq_syn='''_post=s_i*g_i*(V_post-V_i) : amp * meter ** -2 (summed)
        ds_i/dt=-s_i/taud_i+(1-s_i)/taur_i*0.5*(1+tanh(V_pre/10/mV)) : 1
        g_i : siemens * meter**-2
        V_i : volt
        taud_i : second
        taur_i : second
    '''
    #0.006, 0.03
#    g_HTCRE=0.01 * msiemens * cm **-2 #0.4
#    g_HTCRE=0.2 * msiemens * cm **-2 #0.4
#    g_REHTC=1.5 * msiemens * cm **-2
#    g_REHTC=0.2 * msiemens * cm **-2
#    g_AB=1
    
#    HTC_A.J='1 * nA * cmeter ** -2'
#    HTC_B.J='5 * nA * cmeter ** -2'
    HTC_A.J='0 * nA * cmeter ** -2'
    HTC_B.J='0 * nA * cmeter ** -2'
    
#    HTC_A.V='-10*mvolt+20*mvolt*rand()'
    HTC_A.V='0*mvolt+5*mvolt*rand()'
    TC_A.V='-10*mvolt+5*mvolt*rand()'
    RE_A.V='-40*mvolt+5*mvolt*rand()'
    

    g_HTC_E=0.02 * msiemens * cm **-2
    g_TC_E=0.01 * msiemens * cm **-2
    g_E_RE=1 * msiemens * cm **-2

    synHTCA_EA=generate_syn(HTC_A,E_A,'Isyn_LIP','',g_HTC_E,0.125*ms,1*ms,0*mV) 
    synTCA_EA=generate_syn(TC_A,E_A,'Isyn_mdPul','',g_TC_E,0.125*ms,1*ms,0*mV) 
    synHTCB_EB=generate_syn(HTC_B,E_B,'Isyn_LIP','',g_HTC_E,0.125*ms,1*ms,0*mV) 
    synTCB_EB=generate_syn(TC_B,E_B,'Isyn_mdPul','',g_TC_E,0.125*ms,1*ms,0*mV) 
    synEA_REB=generate_syn(E_A,RE_lat_B,'IsynREA','',g_E_RE,0.125*ms,1*ms,0*mV) 
    synEB_REA=generate_syn(E_B,RE_lat_A,'IsynREB','',g_E_RE,0.125*ms,1*ms,0*mV) 


    gGABA_A_RE_RE = 0.07 * msiemens * cm **-2
    gGABA_A_RE_HTC = 0.035 * msiemens * cm **-2
    gGABA_A_RE_TC = 0.01 * msiemens * cm **-2
    gGABA_B_RE_HTC = 0.005 * msiemens * cm **-2
    gGABA_B_RE_TC = 0.02 * msiemens * cm **-2
#    gAMPA_HTC_RE = 0.01 * msiemens * cm **-2
    gAMPA_HTC_RE = 0.02 * msiemens * cm **-2
#    gAMPA_HTC_RE = 0.0225 * msiemens * cm **-2
    gAMPA_TC_RE = 0.02 * msiemens * cm **-2
    
#    gAMPA_HTC_RE = 0 * msiemens * cm **-2
#    gAMPA_TC_RE = 0 * msiemens * cm **-2
#    gGABA_A_RE_TC = 0.1 * msiemens * cm **-2
#    gGABA_B_RE_TC = 0.0575 * msiemens * cm **-2
#    gGABA_A_RE_HTC = 0.1 * msiemens * cm **-2
#    gGABA_B_RE_HTC = 0.0575 * msiemens * cm **-2
    
    synRERE_A=generate_syn(RE_A,RE_A,'IsynRE','',gGABA_A_RE_RE,0.25*ms,5*ms,-80*mV)
    synREHTC_A=generate_syn(RE_A,HTC_A,'IsynREintGABAA','',gGABA_A_RE_HTC,0.25*ms,5*ms,-80*mV)
    synRETC_A=generate_syn(RE_A,TC_A,'IsynREGABAA','',gGABA_A_RE_TC,0.25*ms,5*ms,-80*mV)
    synBREHTC_A=generate_syn(RE_A,HTC_A,'IsynREintGABAB','',gGABA_B_RE_HTC,25*ms,100*ms,-95*mV)
    synBRETC_A=generate_syn(RE_A,TC_A,'IsynREGABAB','',gGABA_B_RE_TC,25*ms,100*ms,-95*mV)
    synHTCRE_A=generate_syn(HTC_A,RE_A,'IsynHTC','',gAMPA_HTC_RE,0.125*ms,1*ms,0*mV)
    synTCRE_A=generate_syn(TC_A,RE_A,'IsynTC','',gAMPA_TC_RE,0.125*ms,1*ms,0*mV)
     
    synRERE_B=generate_syn(RE_B,RE_B,'IsynRE','',gGABA_A_RE_RE,0.25*ms,5*ms,-80*mV)
    synREHTC_B=generate_syn(RE_B,HTC_B,'IsynREintGABAA','',gGABA_A_RE_HTC,0.25*ms,5*ms,-80*mV)
    synRETC_B=generate_syn(RE_B,TC_B,'IsynREGABAA','',gGABA_A_RE_TC,0.25*ms,5*ms,-80*mV)
    synBREHTC_B=generate_syn(RE_B,HTC_B,'IsynREintGABAB','',gGABA_B_RE_HTC,25*ms,100*ms,-95*mV)
    synBRETC_B=generate_syn(RE_B,TC_B,'IsynREGABAB','',gGABA_B_RE_TC,25*ms,100*ms,-95*mV)
    synHTCRE_B=generate_syn(HTC_B,RE_B,'IsynHTC','',gAMPA_HTC_RE,0.125*ms,1*ms,0*mV)
    synTCRE_B=generate_syn(TC_B,RE_B,'IsynTC','',gAMPA_TC_RE,0.125*ms,1*ms,0*mV)


    gGABA_A_RE_lat_HTC = 1.5 * msiemens * cm **-2
    gGABA_B_RE_lat_HTC = 0.05 * msiemens * cm **-2
    synRElatHTC_A=generate_syn(RE_lat_A,HTC_A,'IsynREA','',gGABA_A_RE_lat_HTC,0.25*ms,5*ms,-80*mV)           
    synBRElatHTC_A=generate_syn(RE_lat_A,HTC_A,'IsynREB','',gGABA_B_RE_lat_HTC,25*ms,100*ms,-95*mV)
    synRElatHTC_B=generate_syn(RE_lat_B,HTC_B,'IsynREA','',gGABA_A_RE_lat_HTC,0.25*ms,5*ms,-80*mV)           
    synBRElatHTC_B=generate_syn(RE_lat_B,HTC_B,'IsynREB','',gGABA_B_RE_lat_HTC,25*ms,100*ms,-95*mV)
    
    
    ##Define monitors
    RA=SpikeMonitor(RE_A,record=True)
    RB=SpikeMonitor(RE_B,record=True)
    RlatA=SpikeMonitor(RE_lat_A,record=True)
    RlatB=SpikeMonitor(RE_lat_B,record=True)
    EA=SpikeMonitor(E_A,record=True)
    EB=SpikeMonitor(E_B,record=True)    
    
    all_neurons=all_neuronsA+all_neuronsB+(RE_A,RE_B,RE_lat_A,RE_lat_B,E_A,E_B)
    all_synapses=all_synapsesA,all_synapsesB+(synRERE_A,synREHTC_A,synRETC_A,synHTCRE_A,synTCRE_A,synBREHTC_A,synBRETC_A, synRERE_B,synREHTC_B,synRETC_B,synHTCRE_B,synTCRE_B,synBREHTC_B,synBRETC_B,synTCA_EA,synHTCA_EA,synTCB_EB,synHTCB_EB,synEA_REB,synEB_REA,synRElatHTC_A,synBRElatHTC_A,synRElatHTC_B,synBRElatHTC_B)
    all_monitors=all_monitorsA+all_monitorsB+(RA,RB,RlatA,RlatB,EA,EB)
    all_gap_junctions=all_gap_junctionsA,all_gap_junctionsB
    
    return all_neurons,all_synapses,all_gap_junctions,all_monitors


if __name__=='__main__':
    close('all')
    prefs.codegen.target = 'numpy'
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
#    in_mode='single_spike'
    in_mode='burst'
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
    gapp=0*mamp * cmeter ** -2 # in HTC cells    
    
    net=Network()
    all_neurons,all_synapses,all_gap_junctions,all_monitors=create_mdPul_2col(N_HTC,N_TC,N_RE,condition,in_mode,theta_phase)
    R1A,R2A,V1A,V2A,I1A,I2A,I3A,R1B,R2B,V1B,V2B,I1B,I2B,I3B,RA,RB,RlatA,RlatB,EA,EB=all_monitors
#    print(all_synapses)

    
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
    
    HTC_A,HTC_B=all_neurons[0],all_neurons[2]
    if in_mode=='single_spike':
        HTC_A.delay_steps = [1]  # delay in time steps per neuron
        HTC_B.delay_steps = [1]  # delay in time steps per neuron
        buffer_size = 2  # 1+Maximum delay (in time steps)
    else :
        HTC_A.delay_steps = [3999]  # delay in time steps per neuron
        HTC_B.delay_steps = [3999]  # delay in time steps per neuron
        buffer_size = 4000  # 1+Maximum delay (in time steps)
        
    HTC_A.variables.add_array('voltage_buffer', dimensions=volt.dim, size=(buffer_size, len(HTC_A)))
    HTC_B.variables.add_array('voltage_buffer', dimensions=volt.dim, size=(buffer_size, len(HTC_B)))
    
    init_array=-70*mV*ones((buffer_size, len(HTC_A)))
    HTC_A.voltage_buffer=init_array
    HTC_B.voltage_buffer=init_array
    
    update_code = '''buffer_pointer = (buffer_pointer + 1) % buffer_size
                     voltage_delayed = update_voltage_buffer(V, voltage_buffer, buffer_pointer, delay_steps, buffer_size)'''
       
    buffer_updater_A = HTC_A.run_regularly(update_code, codeobj_class=NumpyCodeObject)
    buffer_updater_B = HTC_B.run_regularly(update_code, codeobj_class=NumpyCodeObject)
        
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

    REA=all_neurons[4]
    M=StateMonitor(REA,['V','ITRE','INa','IK'],record=[0])
    net.add(M)     
            
    M2=StateMonitor(HTC_A,['Isyn','ITHT','INa','IK','ITLT','IH'],record=[0])
    net.add(M2)  
    
    TC_A=all_neurons[1]
    M3=StateMonitor(TC_A,['Isyn','INa','IK','ITLT','IH'],record=[0])
    net.add(M3)  
    
    prefs.codegen.target = 'cython'
    net.run(runtime,report='text',report_period=300*second)
    
    
    figure()
    plot(R1A.t,R1A.i+0,'r.',label='HTC')
    plot(R2A.t,R2A.i+20,'y.',label='TC')
    plot(RA.t,RA.i+100,'b.',label='RE')
    plot(RlatA.t,RlatA.i+200,'k.',label='RE lat')
    plot(EA.t,EA.i+300,'g.',label='E lat')
    
    plot([0,1],[345,345],'k')

    plot(R1B.t,R1B.i+350,'r.')
    plot(R2B.t,R2B.i+370,'y.')
    plot(RB.t,RB.i+450,'b.')
    plot(RlatB.t,RlatB.i+550,'k.')
    plot(EB.t,EB.i+650,'g.')
    
    xlim(0,runtime/second)  
    yticks([150,500],['Object I','Object II'])
    legend()  
    
    figure()
    plot(V1A.t,V1A.V[0],label='HTC V')
    plot(V2A.t,V2A.V[0],label='TC V')
    plot(M.t,M.V[0],label='RE V')
    legend()
    
    figure()
    plot(M.t,M.INa[0],label='RE INa')
    plot(M.t,M.IK[0],label='RE IK')
    plot(M.t,M.ITRE[0],label='RE ITRE')
    legend()
    
    figure()
    plot(M2.t,M2.INa[0],label='HTC INa')
    plot(M2.t,M2.IK[0],label='HTC IK')
    plot(M2.t,M2.ITHT[0],label='HTC ITHT')
    plot(M2.t,M2.ITLT[0],label='HTC ITLT')
    plot(M2.t,M2.IH[0],label='HTC IH')
    plot(M2.t,M2.Isyn[0],label='HTC Isyn')
    legend()
    
    figure()
    plot(M3.t,M3.INa[0],label='HTC INa')
    plot(M3.t,M3.IK[0],label='HTC IK')
    plot(M3.t,M3.ITLT[0],label='HTC ITLT')
    plot(M3.t,M3.IH[0],label='HTC IH')
    plot(M3.t,M3.Isyn[0],label='HTC Isyn')
    legend()
    
    f,Spectrum_LFP_V1=signal.periodogram(V2A.V[0], 100000,'flattop', scaling='spectrum')
    figure()
    plot(f,Spectrum_LFP_V1)
    xlim(0,100)
    clear_cache('cython') 
    
