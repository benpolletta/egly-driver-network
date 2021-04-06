#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19

@author: amelie
"""

from mdPul_two_columns_v3 import *
from FEF_and_LIP import *

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

def create_FEF_and_LIP_and_mdPul(N_RS_vis,N_FS_vis,N_RS_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime,target_time, syn_cond,J,thal, N_HTC,N_TC,N_RE,J_RS_inter):
    all_neurons_FEF,all_synapses_FEF,all_monitors_FEF=create_FEF_full2(N_RS_vis,N_FS_vis,N_RS_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime,target_time)
    R8,R9,R10,V_RS,V_FS,V_SI,R11,R12,R13,R14,mon_RS=all_monitors_FEF
    RSvm_FEF,SIvm_FEF,RSv_FEF,SIv_FEF,VIPv_FEF=all_neurons_FEF[1],all_neurons_FEF[2],all_neurons_FEF[6],all_neurons_FEF[9],all_neurons_FEF[8]
    
    all_neurons_LIP, all_synapses_LIP, all_gap_junctions_LIP, all_monitors_LIP=make_full_network(syn_cond,J,thal,theta_phase)
    V1,V2,V3,R1,R2,R3,I1,I2,I3,V4,R4,I4s,I4a,I4ad,I4bd,R5,R6,R7,V5,V6,V7,inpmon,inpIBmon=all_monitors_LIP
    RS_sup_LIP,IB_LIP,SI_deep_LIP=all_neurons_LIP[0],all_neurons_LIP[5],all_neurons_LIP[9]
    RS_gran_LIP,FS_gran_LIP=all_neurons_LIP[7],all_neurons_LIP[8]

    IB_LIP.ginp_IB=0* msiemens * cm **-2 #the input to RS_sup_LIP is provided with synapses from FEF 
    SI_deep_LIP.ginp_SI=0* msiemens * cm **-2
    RSvm_FEF.ginp_RS=0* msiemens * cm **-2
    SIvm_FEF.ginp_SI=0* msiemens * cm **-2
    RSv_FEF.ginp_RS=0* msiemens * cm **-2
    SIv_FEF.ginp_SI=0* msiemens * cm **-2
    VIPv_FEF.ginp_VIP_good=0* msiemens * cm **-2
    VIPv_FEF.ginp_VIP_bad=0* msiemens * cm **-2
    
    if theta_phase=='good':
        VIP_FEF=all_neurons_FEF[0]
        VIP_FEF.ginp_VIP_good=10* msiemens * cm **-2
        RS_gran_LIP.ginp_RS_good=15* msiemens * cm **-2
        FS_gran_LIP.ginp_FS_good=15* msiemens * cm **-2
        VIP_FEF.ginp_VIP_bad=10* msiemens * cm **-2
        RS_gran_LIP.ginp_RS_bad=15* msiemens * cm **-2
        FS_gran_LIP.ginp_FS_bad=15* msiemens * cm **-2
    if theta_phase=='mixed':
        VIP_FEF=all_neurons_FEF[0]
#        VIP_FEF.ginp_VIP_good=10* msiemens * cm **-2
#        RS_gran_LIP.ginp_RS_good=5* msiemens * cm **-2
#        FS_gran_LIP.ginp_FS_good=5* msiemens * cm **-2
#        RS_gran_LIP.ginp_RS_bad=5* msiemens * cm **-2
#        FS_gran_LIP.ginp_FS_bad=5* msiemens * cm **-2
#        VIP_FEF.ginp_VIP_bad=10* msiemens * cm **-2
        VIP_FEF.ginp_VIP_good=0* msiemens * cm **-2
        RS_gran_LIP.ginp_RS_good=0* msiemens * cm **-2
        FS_gran_LIP.ginp_FS_good=0* msiemens * cm **-2
        RS_gran_LIP.ginp_RS_bad=0* msiemens * cm **-2
        FS_gran_LIP.ginp_FS_bad=0* msiemens * cm **-2
        VIP_FEF.ginp_VIP_bad=0* msiemens * cm **-2
        
    S_FEF_IB_LIP=generate_syn(RSvm_FEF,IB_LIP,'Isyn_FEF','',0*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_FEF_SIdeep_LIP=generate_syn(RSvm_FEF,SI_deep_LIP,'Isyn_FEF','',0.05*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_LIP_RS_FEF=generate_syn(RS_sup_LIP,RSvm_FEF,'Isyn_LIP','',0.009*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
#    S_LIP_FS_FEF=generate_syn(RS_sup_LIP,SIvm_FEF,'Isyn_LIP','',0.009*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
#    S_LIP_RS_FEF=generate_syn(RS_sup_LIP,RSvm_FEF,'Isyn_LIP','',0.015*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    S_LIP_FS_FEF=generate_syn(RS_sup_LIP,SIvm_FEF,'Isyn_LIP','',0.009*msiemens * cm **-2,0.125*ms,1*ms,0*mV)  
#    S_LIP_FS_FEF=generate_syn(RS_sup_LIP,SIvm_FEF,'Isyn_LIP','',0.015*msiemens * cm **-2,0.125*ms,1*ms,0*mV)  
    
#    S_LIP_RSv_FEF=generate_syn(RS_sup_LIP,RSv_FEF,'Isyn_LIP','',0.010*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
#    S_LIP_SIv_FEF=generate_syn(RS_sup_LIP,SIv_FEF,'Isyn_LIP','',0.020*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
#    S_LIP_VIPv_FEF=generate_syn(RS_sup_LIP,VIPv_FEF,'Isyn_LIP','',0.005*msiemens * cm **-2,0.125*ms,1*ms,0*mV)

    #Creating the intermediate between LIP and FEF vis in mdPul
    RS_inter=NeuronGroup(20,eq_RS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    RS_inter.V = '-70*mvolt+10*rand()*mvolt'
    RS_inter.h = '0+0.05*rand()'
    RS_inter.m = '0+0.05*rand()'
    RS_inter.mAR = '0.035+0.025*rand()'
#    RS_inter.J='80 * uA * cmeter ** -2'  #if low LIP->FEF conductance (controlled by SC)
#    RS_inter.J='20 * uA * cmeter ** -2'  #if high LIP->FEF conductance (controlled by SC)
#    RS_inter.run_regularly('J = 80 * (1-0.75*t/second) * uA * cmeter ** -2', dt=0.1*ms)
    RS_inter.J=J_RS_inter
    Rinter=SpikeMonitor(RS_inter,record=True)

    S_LIP_RS_inter_mdPul=generate_syn(RS_sup_LIP,RS_inter,'Isyn_LIP','',0.015*msiemens * cm **-2,0.125*ms,1*ms,0*mV)

#    S_RS_inter_RSv_FEF=generate_syn(RS_inter,RSv_FEF,'Isyn_LIP','',0.010*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
#    S_RS_inter_SIv_FEF=generate_syn(RS_inter,SIv_FEF,'Isyn_LIP','',0.020*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
#    S_RS_inter_VIPv_FEF=generate_syn(RS_inter,VIPv_FEF,'Isyn_LIP','',0.005*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_RS_inter_RSv_FEF=generate_syn(RS_inter,RSv_FEF,'Isyn_LIP','',0.05*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    S_RS_inter_SIv_FEF=generate_syn(RS_inter,SIv_FEF,'Isyn_LIP','',0.05*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    S_RS_inter_VIPv_FEF=generate_syn(RS_inter,VIPv_FEF,'Isyn_LIP','',0.01*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
   
    if target_on:
        RSv_FEF.ginp_RS2=2.5* msiemens * cm **-2
        SIv_FEF.ginp_SI2=2.5* msiemens * cm **-2
        VIPv_FEF.ginp_VIP2=2.5* msiemens * cm **-2
        
    condition='mAChR'
    in_mode='single_spike'
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
    
    all_neurons_mdPul,all_synapses_mdPul,all_gap_junctions_mdPul,all_monitors_mdPul=create_mdPul(N_HTC,N_TC,N_RE,condition,in_mode,theta_phase)

    HTC_A,HTC_B=all_neurons_mdPul[0],all_neurons_mdPul[2]
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

    #Connecting mdPul TC cells to FEF and LIP
    TC_A,TC_B=all_neurons_mdPul[1],all_neurons_mdPul[3]
    S_TC_A_RSg_LIP=generate_syn(TC_A,RS_gran_LIP,'Isyn_FEF','',0.05*msiemens * cm **-2,0.125*ms,1*ms,0*mV)    
    S_TC_A_FSg_LIP=generate_syn(TC_A,FS_gran_LIP,'Isyn_FEF','',0.05*msiemens * cm **-2,0.125*ms,1*ms,0*mV)    
    S_TC_A_VIP_FEF=generate_syn(TC_A,VIP_FEF,'Isyn_mdPul','',0.05*msiemens * cm **-2,0.125*ms,1*ms,0*mV)    


    all_neurons=all_neurons_FEF+all_neurons_LIP+all_neurons_mdPul+(RS_inter,)
    all_synapses=all_synapses_FEF+all_synapses_LIP+all_synapses_mdPul+(S_FEF_IB_LIP,S_FEF_SIdeep_LIP,S_LIP_RS_FEF,S_LIP_FS_FEF,S_LIP_RS_inter_mdPul,S_RS_inter_RSv_FEF,S_RS_inter_SIv_FEF,S_RS_inter_VIPv_FEF)+(S_TC_A_RSg_LIP,S_TC_A_FSg_LIP,S_TC_A_VIP_FEF)
    all_monitors=all_monitors_FEF+all_monitors_LIP+all_monitors_mdPul+(Rinter,)
    all_gap_junctions=all_gap_junctions_LIP,all_gap_junctions_mdPul
    
    return all_neurons,all_synapses,all_gap_junctions,all_monitors

def save_raster(name,raster_i,raster_t,path):
    raster_file=open(path+'/raster_'+name+'_i.txt','w')
    for elem in raster_i:
        raster_file.write(str(elem)+',')
    raster_file.close()
    raster_file=open(path+'/raster_'+name+'_t.txt','w')
    for elem in raster_t:
        raster_file.write(str(elem)+',')
    raster_file.close()
    return

def run_simu(simu,path):
    prefs.codegen.target = 'numpy'
    
    target_time,J_RS_inter,N_simu=simu
    
    new_path=path+"/results_"+str(N_simu)
    os.mkdir(new_path)
    
    theta_phase='mixed'
    target_on=True  
#    target_on=False 
#    target_time = 500*msecond
#    target_time = 625*msecond
#    target_time = 575*msecond #good to bad
#    target_time = 450*msecond #bad to good
    target_time = 800*msecond
    
    start_scope()
#    set_device('genn')

    runtime=1*second
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
#    taurinp=1*ms
#    taudinp=5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    Vhigh2=0*mV
    Vlow2=-80*mV
    
    NN=1 #multiplicative factor on the number of neurons
    N_RS,N_FS,N_SI,N_IB= NN*80,NN*20,NN*20,NN*20 #Number of neurons of RE, TC, and HTC type
    N_SI,N_RS_gran,N_SI_gran=20,20,20
    N_RS_vis,N_FS_vis,N_RS_mot,N_SI_mot,N_dSI_vm,N_RS_vm,N_gSI_vm=[20]*7

    
    all_SIdFSg=[2*msiemens * cm **-2] #1
    all_FSgRSg=[1* msiemens * cm **-2]
    all_RSgFSg=[1*msiemens * cm **-2]
    all_RSgRSg=[0.5*msiemens * cm **-2]
    all_FSgFSg=[0.3* msiemens * cm **-2]
    all_RSgRSs=[2*msiemens * cm **-2]
    all_RSgFSs=[0.1*msiemens * cm **-2]
    all_FSgRSs=[0.1* msiemens * cm **-2]
    all_J_RSg=['10 * uA * cmeter ** -2']
    all_J_FSg=['-5 * uA * cmeter ** -2']
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
        
    N_HTC, N_TC,N_RE= 20,80,100
    condition='mAChR'
    in_mode='single_spike'
    if in_mode=='single_spike':
        buffer_size = 2  # 1+Maximum delay (in time steps)
    else :
        buffer_size = 4000  # 1+Maximum delay (in time steps)
    @check_units(V=volt, voltage_buffer=volt, buffer_pointer=1, delay_steps=1, buffer_size=1, result=volt)
    def update_voltage_buffer(V, voltage_buffer, buffer_pointer, delay_steps, buffer_size):
        # Write current rate into the buffer
        voltage_buffer[buffer_pointer, :] = V
        # Get delayed rates 
        rows = (buffer_pointer - delay_steps) % buffer_size    
        return voltage_buffer[rows, arange(len(rows))]
            
    
    print('Network setup')
    
    net=Network()    

    all_neurons,all_synapses,all_gap_junctions,all_monitors=create_FEF_and_LIP_and_mdPul(N_RS_vis,N_FS_vis,N_RS_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime,target_time, syn_cond,J,thal, N_HTC,N_TC,N_RE,J_RS_inter)    
    
    print('Compiling with cython')
    
    prefs.codegen.target = 'cython' #cython=faster, numpy = default python
#    set_device('genn')
#    defaultclock.dt = 0.01*ms
    
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
#    taurinp2=2*ms
#    taudinp2=10*ms
#    tauinp2=taudinp2
    taurinp2=0.1*ms
    taudinp2=0.5*ms
    tauinp2=taudinp2
    
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_gap_junctions)
    net.add(all_monitors)
    net.run(runtime,report='text',report_period=300*second)
    
    R8,R9,R10,V_RS,V_FS,V_SI,R11,R12,R13,R14,mon_RS, V1,V2,V3,R1,R2,R3,I1,I2,I3,V4,R4,I4s,I4a,I4ad,I4bd,R5,R6,R7,V5,V6,V7,inpmon,inpIBmon,R1A,R2A,V1A,V2A,I1A,I2A,I3A,R1B,R2B,V1B,V2B,I1B,I2B,I3B,RA,RB,Rinter=all_monitors
    
    save_raster('LIP RS',R1.i,R1.t,new_path)
    save_raster('LIP FS',R2.i,R2.t,new_path)
    save_raster('LIP SI',R3.i,R3.t,new_path)
    save_raster('LIP IB',R4.i,R4.t,new_path)
    save_raster('LIP RS gran',R5.i,R5.t,new_path)
    save_raster('LIP FS gran',R6.i,R6.t,new_path)
    save_raster('LIP SI deep',R7.i,R7.t,new_path)
    save_raster('FEF RS vm',R8.i,R8.t,new_path)
    save_raster('FEF SI2 vm',R9.i,R9.t,new_path)
    save_raster('FEF SI1 vm',R10.i,R10.t,new_path)
    save_raster('FEF RS v',R11.i,R11.t,new_path)
    save_raster('FEF FS v',R12.i,R12.t,new_path)
    save_raster('FEF VIP v',R13.i,R13.t,new_path)
    save_raster('FEF SI v',R14.i,R14.t,new_path)
    save_raster('FEF RS m',mon_RS.i,mon_RS.t,new_path)

    save_raster('HTC A',R1A.i,R1A.t,new_path)
    save_raster('TC A',R2A.i,R2A.t,new_path)
    save_raster('RE A',RA.i,RA.t,new_path)
    save_raster('HTC B',R1B.i,R1B.t,new_path)
    save_raster('TC B',R2B.i,R2B.t,new_path)
    save_raster('RE B',RB.i,RB.t,new_path)
    save_raster('R_inter',Rinter.i,Rinter.t,new_path) 

#    # LIP Plots
#    figure(figsize=(6,9))
##    subplot(411)
#    up=100
#    plot(R1.t,R1.i+140+up,'r.',label='RS cells')
#    plot(R2.t,R2.i+120+up,'b.',label='FS cells')
#    plot(R3.t,R3.i+100+up,'g.',label='SI cells')
#    plot([0.2,runtime/second],[95+up,95+up],'k--')
#    plot(R5.t,R5.i+70+up,'r.',label='Granular RS')
#    plot(R6.t,R6.i+50+up,'b.',label='Granular FS')
#    plot([0.2,runtime/second],[45+up,45+up],'k--')
#    plot(R4.t,R4.i+20+up,'m.',label='IB cells')
#    plot(R7.t,R7.i+up,'g.',label='Deep SI')
#    xlim(0.2,runtime/second)
#    plot([0.2,runtime/second],[up-10,up-10],'k')
#
##    subplot(412)
#    up=0
##    title('FEF Visual-Motor Neurons')
#    plot(R10.t,R10.i+0+up,'k.',label='VIP')
#    plot(R8.t,R8.i+60+up,'r.',label='RS')
#    plot(R9.t,R9.i+40+up,'g.',label='SOM')
#    xlim(0.2,runtime/second)
##    plot([0.2,runtime/second],[up-10,up-10],'k')
##    xticks([],[])
##    legend(loc='upper left') 
#    #xlabel('Time (s)')
##    ylabel('Neuron index')
#    xlabel('Time (s)')
#    ylabel('Neuron index')
#
#    #FEF Plots    
#    figure(figsize=(10,4))
#    subplot(131)
#    title('Visual Neurons')
#    plot(R11.t,R11.i+0,'r.',label='RS')
#    plot(R12.t,R12.i+20,'b.',label='FS')
#    plot(R13.t,R13.i+40,'k.',label='VIP')
#    plot(R14.t,R14.i+60,'g.',label='SOM')
#    xlim(0.2,runtime/second)
#    legend(loc='upper left')   
#    xlabel('Time (s)')
#    ylabel('Neuron index')
#    
#    subplot(132)
#    title('Visual-Motor Neurons')
#    plot(R10.t,R10.i+0,'g.',label='VIP')
#    plot(R8.t,R8.i+60,'r.',label='RS')
#    plot(R9.t,R9.i+40,'.',label='SOM',color='lime')
#    xlim(0.2,runtime/second)
#    legend(loc='upper left') 
#    xlabel('Time (s)')
#    ylabel('Neuron index')
#    
#    subplot(133)
#    title('Decision cells')
#    plot(mon_RS.t,mon_RS.i+0,'r.',label='RS')
#    xlim(0.2,runtime/second)
#    legend(loc='upper left') 
#    xlabel('Time (s)')
#    ylabel('Neuron index')
#    
#    tight_layout()
#    
#    
#    figure(figsize=(5,5))
##    subplot(411)
#    up=0
#    plot(R1.t,R1.i+140+up,'r.',label='RS cells')
#    plot(R2.t,R2.i+120+up,'b.',label='FS cells')
#    plot(R3.t,R3.i+100+up,'g.',label='SI cells')
#    plot([0.2,runtime/second],[95+up,95+up],'k--')
#    plot(R5.t,R5.i+70+up,'r.',label='Granular RS')
#    plot(R6.t,R6.i+50+up,'b.',label='Granular FS')
#    plot([0.2,runtime/second],[45+up,45+up],'k--')
#    plot(R4.t,R4.i+20+up,'m.',label='IB cells')
#    plot(R7.t,R7.i+up,'g.',label='Deep SI')
#    xlim(0.2,runtime/second)
#
#    xlabel('Time (s)')
#    ylabel('Neuron index')
#    
#    figure()
#    plot(R1A.t,R1A.i+0,'r.',label='HTC')
#    plot(R2A.t,R2A.i+20,'y.',label='TC')
#    plot(RA.t,RA.i+100,'b.',label='RE lat')
#    
#    plot([0,1],[225,225],'k')
#
#    plot(R1B.t,R1B.i+250,'r.')
#    plot(R2B.t,R2B.i+270,'y.')
#    plot(RB.t,RB.i+350,'b.')
#    
#    xlim(0,runtime/second)  
#    yticks([150,500],['Object I','Object II'])
#    legend()  
#    
#    figure()
#    plot(Rinter.t,Rinter.i+0,'r.',label='Rinter')
#    xlim(0,runtime/second)
#    legend()  
#    
#    clear_cache('cython')
