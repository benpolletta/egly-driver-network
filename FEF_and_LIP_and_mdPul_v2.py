#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 16:23:57 2020

@author: amelie
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 13:56:08 2020

@author: amelie
"""

from brian2 import *

from scipy import signal

from FEF_full import *
from LIP_full import *
from mdPul_two_columns import *

from itertools import *


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

if __name__=='__main__':
    prefs.codegen.target = 'numpy'
    theta_phase='bad'
    target_on=True    
    
    start_scope()
    close('all')

    runtime=2*second
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    
    NN=1 #multiplicative factor on the number of neurons
    N_RS,N_FS,N_SI,N_IB= NN*80,NN*20,NN*20,NN*20 #Number of neurons of RE, TC, and HTC type
    N_SI,N_RS_gran,N_SI_gran=20,20,20
    N_RS_vis,N_FS_vis,N_RS_mot,N_SI_mot,N_dSI_vm,N_RS_vm,N_gSI_vm=[20]*7
    N_HTC,N_TC,N_RE= 20,80,100 #Number of neurons of RE, TC, and HTC type
    
    all_SIdFSg=[2*msiemens * cm **-2] #1
    all_FSgRSg=[1* msiemens * cm **-2]
    all_RSgFSg=[1*msiemens * cm **-2]
    all_RSgRSg=[0.3*msiemens * cm **-2]
    all_FSgFSg=[0.3* msiemens * cm **-2]
    all_RSgRSs=[2*msiemens * cm **-2]
    all_RSgFSs=[0.1*msiemens * cm **-2]
    all_FSgRSs=[0.1* msiemens * cm **-2]
    all_J_RSg=['30 * uA * cmeter ** -2']
    all_J_FSg=['5 * uA * cmeter ** -2']
    all_thal=[10* msiemens * cm **-2]
    thal=all_thal[0]
    
    all_syn_cond=list(product(all_SIdFSg,all_FSgRSg,all_RSgFSg,all_RSgRSg,all_FSgFSg,all_RSgRSs,all_RSgFSs,all_FSgRSs))
    all_J=list(product(all_J_RSg,all_J_FSg))
    syn_cond=all_syn_cond[0]
    J=all_J[0]
    
    if theta_phase=='bad':
        input_beta2_IB=False
        input_beta2_RS=False
        input_beta2_FS_SI=True
        input_thalamus_gran=True
        gFS=0* msiemens * cm **-2
        ginp_SI=0* msiemens * cm **-2
        ginpSIdeep=0* msiemens * cm **-2
        thal_cond=2* msiemens * cm **-2
        kainate='low'
        
    if theta_phase=='good':
#        input_beta2_IB=True
        input_beta2_IB=False
        ginp_IB=500* msiemens * cm **-2
        ginpSIdeep=500* msiemens * cm **-2
        input_beta2_RS=False
        input_beta2_FS_SI=False
        input_thalamus_gran=True
        thal_cond=thal
        kainate='low'
        
    if theta_phase=='mixed':
        input_mixed=True
        ginp_IB=500* msiemens * cm **-2
        ginpSIdeep=500* msiemens * cm **-2
        input_beta2_IB=False
        input_beta2_RS=False
        input_beta2_RS=False
        input_beta2_FS_SI=False
        input_thalamus_gran=False
        kainate='low'
        
    Vrev_inp2=0*mV
    taurinp2=0.1*ms
    taudinp2=0.5*ms
    tauinp2=taudinp2
    Vhigh2=0*mV
    Vlow2=-80*mV
        
    condition='mAChR'
    in_mode='single_spike'
#    in_mode='burst'
    theta_phase='mixed'

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
    
    print('Network setup')
    
    net=Network()
    
    all_neurons_FEF,all_synapses_FEF,all_monitors_FEF=create_network_no_motor2(N_RS_vis,N_FS_vis,N_RS_mot,N_SI_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime)
    R1FEF,R2FEF,R3FEF,V1FEF,V2FEF,V3FEF,R4FEF,R5FEF,V4FEF,V5FEF,mon_FS=all_monitors_FEF
    RSvm_FEF,SIvm_FEF=all_neurons_FEF[1],all_neurons_FEF[2]
    SI2vm_FEF=all_neurons_FEF[0]
    
    all_neurons_LIP, all_synapses_LIP, all_gap_junctions_LIP, all_monitors_LIP=make_full_network(syn_cond,J,thal,theta_phase)
    V1,V2,V3,R1,R2,R3,I1,I2,I3,V4,R4,I4s,I4a,I4ad,I4bd,R5,R6,R7,V5,V6,V7,inpmon,inpIBmon=all_monitors_LIP
    RS_sup_LIP,IB_LIP,SI_deep_LIP=all_neurons_LIP[0],all_neurons_LIP[5],all_neurons_LIP[9]
    RS_gran_LIP,FS_gran_LIP=all_neurons_LIP[7],all_neurons_LIP[8]
    
    all_neurons_mdPul,all_synapses_mdPul,all_gap_junctions_mdPul,all_monitors_mdPul=create_mdPul(N_HTC,N_TC,N_RE,condition,in_mode,theta_phase)
    R1A,R2A,R3A,V1A,V2A,V3A,I1A,I2A,R1B,R2B,R3B,V1B,V2B,V3B,I1B,I2B,RA,RB=all_monitors_mdPul
    TC_B,RE_B=all_neurons_mdPul[1],all_neurons_mdPul[2]     
    
    IB_LIP.ginp_IB=0* msiemens * cm **-2 #the input to RS_sup_LIP is provided with synapses from FEF 
    SI_deep_LIP.ginp_SI=0* msiemens * cm **-2
    SI2vm_FEF.ginp_VIP_good=0* msiemens * cm **-2
    SI2vm_FEF.ginp_VIP_bad=0* msiemens * cm **-2
    RS_gran_LIP.ginp_RS_good=0* msiemens * cm **-2 #5
    FS_gran_LIP.ginp_FS_good=0* msiemens * cm **-2 #5
    RS_gran_LIP.ginp_RS_bad=0* msiemens * cm **-2 #5
    FS_gran_LIP.ginp_FS_bad=0* msiemens * cm **-2 #5
#    if theta_phase=='good' or theta_phase=='mixed':
#        RSvm_FEF.ginp_RS=10* msiemens * cm **-2
#        SIvm_FEF.ginp_SI=10* msiemens * cm **-2

    RSvm_FEF.ginp_RS=0* msiemens * cm **-2 #10
    SIvm_FEF.ginp_SI=0* msiemens * cm **-2 #10
    
    HTC_A,HTC_B=all_neurons_mdPul[0],all_neurons_mdPul[6]
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
    
    net.add(all_neurons_FEF)
    net.add(all_synapses_FEF)
    net.add(all_monitors_FEF)    
    
    net.add(all_neurons_LIP)
    net.add(all_synapses_LIP)
    net.add(all_gap_junctions_LIP)
    net.add(all_monitors_LIP)
    
    net.add(all_neurons_mdPul)
    net.add(all_synapses_mdPul)
    net.add(all_gap_junctions_mdPul)
    net.add(all_monitors_mdPul)
    
    S_FEF_IB_LIP=generate_syn(RSvm_FEF,IB_LIP,'Isyn_FEF','',0*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_FEF_SIdeep_LIP=generate_syn(RSvm_FEF,SI_deep_LIP,'Isyn_FEF','',0.03*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_LIP_RS_FEF=generate_syn(RS_sup_LIP,RSvm_FEF,'Isyn_LIP','',0.004*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    S_LIP_FS_FEF=generate_syn(RS_sup_LIP,SIvm_FEF,'Isyn_LIP','',0.004*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    
    S_FEF_mdPul=generate_syn(RSvm_FEF,TC_B,'Isyn_FEF','',0*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_LIP_mdPul=generate_syn(IB_LIP,RE_B,'Isyn_LIP','',0.001*msiemens * cm **-2,0.125*ms,1*ms,0*mV)    

    S_mdPul_FEF_VIP=generate_syn(TC_B,SI2vm_FEF,'Isyn_mdPul','',0.01*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_mdPul_LIP_RSg=generate_syn(TC_B,RS_gran_LIP,'Isyn_mdPul','',0.01*msiemens * cm **-2,0.125*ms,1*ms,0*mV)    
    S_mdPul_LIP_FSg=generate_syn(TC_B,FS_gran_LIP,'Isyn_mdPul','',0.01*msiemens * cm **-2,0.125*ms,1*ms,0*mV)    
  
    
    net.add(S_FEF_IB_LIP)
    net.add(S_FEF_SIdeep_LIP)
    net.add(S_LIP_RS_FEF)
    net.add(S_LIP_FS_FEF)
    net.add([S_FEF_mdPul,S_LIP_mdPul])
    net.add([S_mdPul_FEF_VIP,S_mdPul_LIP_RSg,S_mdPul_LIP_FSg])
    
    print('Compiling with cython')
    
    prefs.codegen.target = 'cython' #cython=faster, numpy = default python
    
    net.run(runtime,report='text',report_period=300*second)
    
    
    
    # LIP Plots
    figure()
    plot(R1.t,R1.i+140,'r.',label='RS cells')
    plot(R2.t,R2.i+120,'m.',label='FS cells')
    plot(R3.t,R3.i+100,'y.',label='SI cells')
    plot(R5.t,R5.i+70,'g.',label='Granular RS')
    plot(R6.t,R6.i+50,'c.',label='Granular FS')
    plot(R4.t,R4.i+20,'b.',label='IB cells')
    plot(R7.t,R7.i,'k.',label='Deep SI')
    xlim(0,runtime/second)
    legend(loc='upper left')
    
#    min_t=int(50*ms*100000*Hz)
#    LFP_V_RS=1/N_RS*sum(V1.V,axis=0)[min_t:]
#    LFP_V_FS=1/N_FS*sum(V2.V,axis=0)[min_t:]
#    LFP_V_SI=1/N_SI*sum(V3.V,axis=0)[min_t:]
#    LFP_V_IB=1/N_IB*sum(V4.V,axis=0)[min_t:]
#    LFP_V_RSg=1/N_FS*sum(V5.V,axis=0)[min_t:]
#    LFP_V_FSg=1/N_FS*sum(V6.V,axis=0)[min_t:]
#    LFP_V_SId=1/N_SI*sum(V7.V,axis=0)[min_t:]
#    
#    f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_FS=signal.periodogram(LFP_V_FS, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_SI=signal.periodogram(LFP_V_SI, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_IB=signal.periodogram(LFP_V_IB, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_RSg=signal.periodogram(LFP_V_RSg, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_FSg=signal.periodogram(LFP_V_FSg, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_SId=signal.periodogram(LFP_V_SId, 100000,'flattop', scaling='spectrum')
#    
#    figure(figsize=(10,8))    
#    subplot(421)
#    plot(f,Spectrum_LFP_V_RS)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('RS cell')
#    subplot(422)
#    plot(f,Spectrum_LFP_V_FS)
#    yticks([],[])
#    xlim(0,100)
#    title('FS cell')
#    subplot(423)
#    plot(f,Spectrum_LFP_V_SI)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('SI cell')
#    subplot(425)
#    plot(f,Spectrum_LFP_V_RSg)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('gran RS cell')
#    subplot(426)
#    plot(f,Spectrum_LFP_V_FSg)
#    yticks([],[])
#    xlim(0,100)
#    title('gran FS cell')
#    subplot(427)
#    plot(f,Spectrum_LFP_V_IB)
#    xlabel('Frequency (Hz)')
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('IB cell')
#    subplot(428)
#    plot(f,Spectrum_LFP_V_SId)
#    yticks([],[])
#    xlim(0,100)
#    xlabel('Frequency (Hz)')
#    title('deep SI cell')
#    
#    tight_layout()
    
    
    
    #FEF Plots
    figure(figsize=(10,4))
    subplot(121)
    title('Visual Neurons')
    plot(R4FEF.t,R4FEF.i+20,'r.',label='RS')
    plot(R5FEF.t,R5FEF.i+0,'b.',label='FS')
    xlim(0,runtime/second)
    legend(loc='upper left')   
    
    subplot(122)
    title('Visual-Motor Neurons')
    plot(R3FEF.t,R3FEF.i+0,'c.',label='VIP')
    plot(R1FEF.t,R1FEF.i+60,'r.',label='RS')
    plot(R2FEF.t,R2FEF.i+40,'b.',label='SI')
    xlim(0,runtime/second)
    legend(loc='upper left') 
    
#    subplot(133)
#    title('Motor Neurons')
#    plot(R6.t,R6.i+60,'r.',label='RS')
#    plot(R7.t,R7.i+40,'b.',label='SI')
#    plot(R8.t,R8.i+0,'c.',label='Fix')
#    xlim(0,runtime/second)
#    legend(loc='upper left') 
    
    
#    min_t=int(50*ms*100000*Hz)
#    LFP_V1=1/20*sum(V1FEF.V,axis=0)[min_t:]
#    LFP_V2=1/20*sum(V2FEF.V,axis=0)[min_t:]
#    LFP_V3=1/20*sum(V3FEF.V,axis=0)[min_t:]
#    LFP_V4=1/20*sum(V4FEF.V,axis=0)[min_t:]
#    LFP_V5=1/20*sum(V5FEF.V,axis=0)[min_t:]
##    LFP_V6=1/20*sum(V6.V,axis=0)[min_t:]
##    LFP_V7=1/20*sum(V7.V,axis=0)[min_t:]
#    
#    f,Spectrum_LFP_V1=signal.periodogram(LFP_V1, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V2=signal.periodogram(LFP_V2, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V3=signal.periodogram(LFP_V3, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V4=signal.periodogram(LFP_V4, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V5=signal.periodogram(LFP_V5, 100000,'flattop', scaling='spectrum')
##    f,Spectrum_LFP_V6=signal.periodogram(LFP_V6, 100000,'flattop', scaling='spectrum')
##    f,Spectrum_LFP_V7=signal.periodogram(LFP_V7, 100000,'flattop', scaling='spectrum')
#
#    figure(figsize=(10,4))
#    subplot(321)
#    plot(f,Spectrum_LFP_V4)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('visual RS')
#    subplot(323)
#    plot(f,Spectrum_LFP_V5)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('visual FS')  
#    
#    subplot(322)
#    plot(f,Spectrum_LFP_V1)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('visual-motor gran RS')
#    subplot(324)
#    plot(f,Spectrum_LFP_V2)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('visual-motor gran SI')  
#    subplot(326)
#    plot(f,Spectrum_LFP_V3)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('visual-motor deep SI')    
    
    #mdPul plots
    
    figure()
    plot(R1A.t,R1A.i+0,'r.',label='HTC')
    plot(R2A.t,R2A.i+20,'y.',label='TC')
    plot(R3A.t,R3A.i+100,'g.',label='RE int')
    plot(RA.t,RA.i+200,'b.',label='RE lat')
    
    plot([0,1],[325,325],'k')

    plot(R1B.t,R1B.i+350,'r.')
    plot(R2B.t,R2B.i+370,'y.')
    plot(R3B.t,R3B.i+450,'g.')
    plot(RB.t,RB.i+550,'b.')
    
    xlim(0,runtime/second)  
    yticks([150,500],['Object I','Object II'])
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
