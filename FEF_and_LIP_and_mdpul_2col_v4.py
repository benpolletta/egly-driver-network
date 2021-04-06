#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 09:07:50 2021

@author: amelie
"""

from FEF_and_LIP import *
from mdPul_two_columns_v4 import *
from itertools import *

eq_syn='''_post=s_i*g_i*(V_post-V_i) : amp * meter ** -2 (summed)
    ds_i/dt=-s_i/taud_i+(1-s_i)/taur_i*0.5*(1+tanh(V_pre/10/mV)) : 1
    g_i : siemens * meter**-2
    V_i : volt
    taud_i : second
    taur_i : second
'''

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

def create_FEF_and_LIP_and_mdPul(N_RS_vis,N_FS_vis,N_RS_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime,target_time, syn_cond,J,thal, N_HTC,N_TC,N_RE):
    all_neurons_FEF,all_synapses_FEF,all_monitors_FEF=create_FEF_full2(N_RS_vis,N_FS_vis,N_RS_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime,target_time)
    R8,R9,R10,V_RS,V_FS,V_SI,R11,R12,R13,R14,mon_RS=all_monitors_FEF
    RSvm_FEF,SIvm_FEF,RSv_FEF,SIv_FEF,VIPv_FEF,RS_decision=all_neurons_FEF[1],all_neurons_FEF[2],all_neurons_FEF[6],all_neurons_FEF[9],all_neurons_FEF[8],all_neurons_FEF[-1]
    
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
    
    TC_inter=NeuronGroup(20,eq_TC_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    TC_inter.V = '-45*mvolt+20*mvolt*rand()'
    TC_inter.J='-0.25 * uA * cmeter ** -2' #5 or 7=high LIP-FEF conductance
    TC_inter.Ca_i = '1e-7 * mole * metre**-3'
    TC_inter.run_regularly('J = -0.25 * (0.5+0.5*t/second) * uA * cmeter ** -2', dt=0.1*ms)
    
    TCinter=SpikeMonitor(TC_inter,record=True)
    
    
    #IB cells in SC
    IB_SC_soma=NeuronGroup(N_IB,eq_IB_soma,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    IB_SC_soma.V = '-100*mvolt+10*rand()*mvolt'
    IB_SC_soma.h = '0+0.05*rand()'
    IB_SC_soma.m = '0+0.05*rand()'
#    IB_SC_soma.J='-4.5 * uA * cmeter ** -2' #article SI=-3.5, code=-4.5
    IB_SC_soma.J='0 * uA * cmeter ** -2' #article SI=-3.5, code=-4.5
    
    IB_SC_axon=NeuronGroup(N_IB,eq_IB_axon,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    IB_SC_axon.V = '-100*mvolt+10*rand()*mvolt'
    IB_SC_axon.h = '0+0.05*rand()'
    IB_SC_axon.m = '0+0.05*rand()'
    IB_SC_axon.mKM = '0+0.05*rand()'
#    IB_SC_axon.J='-0.4 * uA * cmeter ** -2' #article SI=+0.1, code=-0.4
    IB_SC_axon.J='0 * uA * cmeter ** -2' #article SI=+0.1, code=-0.4
    
    IB_SC_ad=NeuronGroup(N_IB,eq_IB_ad,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    IB_SC_ad.V = '-100*mvolt+10*rand()*mvolt'
    IB_SC_ad.h = '0+0.05*rand()'
    IB_SC_ad.m = '0+0.05*rand()'
    IB_SC_ad.mAR = '0+0.001*rand()'
    IB_SC_ad.mKM = '0+0.05*rand()'
    IB_SC_ad.mCaH = '0+0.01*rand()'
    IB_SC_ad.J='25.5 * uA * cmeter ** -2'  #article SI=27.5, code=25.5
    
    IB_SC_bd=NeuronGroup(N_IB,eq_IB_bd,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    IB_SC_bd.V = '-100*mvolt+10*rand()*mvolt'
    IB_SC_bd.h = '0+0.05*rand()'
    IB_SC_bd.m = '0+0.05*rand()'
    IB_SC_bd.mAR = '0+0.001*rand()'
    IB_SC_bd.mKM = '0+0.05*rand()'
    IB_SC_bd.mCaH = '0+0.01*rand()'
    IB_SC_bd.J='42.5 * uA * cmeter ** -2' #article SI=44.5, code=42.5
    
    S_IBIB=Synapses(IB_SC_axon,IB_SC_bd,model='IsynIB_LIP'+eq_syn,method='exact')
    S_IBIB.connect()
    S_IBIB.g_i=1/500* msiemens * cm **-2
    S_IBIB.taur_i=0.25*ms  #0.5 in Mark
    S_IBIB.taud_i=100*ms
    S_IBIB.V_i=0*mV
    
    
    eq_gap='''_post=g_i*(V_post-V_pre) : amp * meter ** -2 (summed)
        g_i : siemens * meter**-2
    '''
    
    gapIB_SomaAd=Synapses(IB_SC_soma,IB_SC_ad,model='Igap_soma'+eq_gap,method='exact')
    gapIB_SomaAd.connect(j='i')
    gapIB_SomaAd.g_i=0.2* msiemens * cm **-2
    
    gapIB_SomaBd=Synapses(IB_SC_soma,IB_SC_bd,model='Igap_soma'+eq_gap,method='exact')
    gapIB_SomaBd.connect(j='i')
    gapIB_SomaBd.g_i=0.2* msiemens * cm **-2
    
    gapIB_SomaAxon=Synapses(IB_SC_soma,IB_SC_axon,model='Igap_soma'+eq_gap,method='exact')
    gapIB_SomaAxon.connect(j='i')
    gapIB_SomaAxon.g_i=0.3* msiemens * cm **-2
    
    gapIB_AdSoma=Synapses(IB_SC_ad,IB_SC_soma,model='Igap_ad'+eq_gap,method='exact')
    gapIB_AdSoma.connect(j='i')
    gapIB_AdSoma.g_i=0.4* msiemens * cm **-2
    
    gapIB_BdSoma=Synapses(IB_SC_bd,IB_SC_soma,model='Igap_bd'+eq_gap,method='exact')
    gapIB_BdSoma.connect(j='i')
    gapIB_BdSoma.g_i=0.4* msiemens * cm **-2
    
    gapIB_AxonSoma=Synapses(IB_SC_axon,IB_SC_soma,model='Igap_axon'+eq_gap,method='exact')
    gapIB_AxonSoma.connect(j='i')
    gapIB_AxonSoma.g_i=0.3* msiemens * cm **-2
    
    gap_IBIB=Synapses(IB_SC_axon,IB_SC_axon,model='Igap_axon'+eq_gap,method='exact')
    gap_IBIB.connect()
    gap_IBIB.g_i=0.0025* msiemens * cm **-2
    
    IBSC=SpikeMonitor(IB_SC_soma,record=True)
    
#    RS_inter=NeuronGroup(20,eq_RS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
#    RS_inter.V = '-70*mvolt+10*rand()*mvolt'
#    RS_inter.h = '0+0.05*rand()'
#    RS_inter.m = '0+0.05*rand()'
#    RS_inter.mAR = '0.035+0.025*rand()'
#    RS_inter.J='60 * uA * cmeter ** -2'  #if low LIP->FEF conductance (controlled by SC)
#    RS_inter.J='20 * uA * cmeter ** -2'  #if high LIP->FEF conductance (controlled by SC)
#    RS_inter.run_regularly('J = 80 * (1-0.75*t/second) * uA * cmeter ** -2', dt=0.1*ms)

   
    if target_on:
        RSv_FEF.ginp_RS2=2.5* msiemens * cm **-2
        SIv_FEF.ginp_SI2=2.5* msiemens * cm **-2
        VIPv_FEF.ginp_VIP2=2.5* msiemens * cm **-2 
    
    all_neurons_mdPul,all_synapses_mdPul,all_gap_junctions_mdPul,all_monitors_mdPul=create_mdPul_2col(N_HTC,N_TC,N_RE,condition,in_mode,theta_phase)


#    HTC_A,HTC_B=all_neurons_mdPul[0],all_neurons_mdPul[2]
##    RE_A,RE_B=all_neurons_mdPul[4],all_neurons_mdPul[5]
##    HTC_A.V='-10*mvolt+5*mvolt*rand()'
##    RE_A.V='-20*mvolt+5*mvolt*rand()'
###    RE_A.V='-50*mvolt+5*mvolt*rand()'
##    HTC_B.V='-60*mvolt+5*mvolt*rand()'
##    RE_B.V='-60*mvolt+5*mvolt*rand()'
#    if in_mode=='single_spike':
#        HTC_A.delay_steps = [1]  # delay in time steps per neuron
#        HTC_B.delay_steps = [1]  # delay in time steps per neuron
#        buffer_size = 2  # 1+Maximum delay (in time steps)
#    else :
#        HTC_A.delay_steps = [3999]  # delay in time steps per neuron
#        HTC_B.delay_steps = [3999]  # delay in time steps per neuron
#        buffer_size = 4000  # 1+Maximum delay (in time steps)
#        
#    HTC_A.variables.add_array('voltage_buffer', dimensions=volt.dim, size=(buffer_size, len(HTC_A)))
#    HTC_B.variables.add_array('voltage_buffer', dimensions=volt.dim, size=(buffer_size, len(HTC_B)))
#    
#    init_array=-70*mV*ones((buffer_size, len(HTC_A)))
#    HTC_A.voltage_buffer=init_array
#    HTC_B.voltage_buffer=init_array
#    
#    
#    
#    update_code = '''buffer_pointer = (buffer_pointer + 1) % buffer_size
#                     voltage_delayed = update_voltage_buffer(V, voltage_buffer, buffer_pointer, delay_steps, buffer_size)'''
#       
#    buffer_updater_A = HTC_A.run_regularly(update_code, codeobj_class=NumpyCodeObject)
#    buffer_updater_B = HTC_B.run_regularly(update_code, codeobj_class=NumpyCodeObject)
#        
#    @check_units(V=volt, voltage_buffer=volt, buffer_pointer=1, delay_steps=1, buffer_size=1, result=volt)
#    def update_voltage_buffer(V, voltage_buffer, buffer_pointer, delay_steps, buffer_size):
#        # Write current rate into the buffer
#        voltage_buffer[buffer_pointer, :] = V
#        # Get delayed rates 
#        rows = (buffer_pointer - delay_steps) % buffer_size    
#        return voltage_buffer[rows, arange(len(rows))]

#    Connecting mdPul TC cells to FEF and LIP
    HTC_A,HTC_B=all_neurons_mdPul[0],all_neurons_mdPul[2]
    TC_A,TC_B=all_neurons_mdPul[1],all_neurons_mdPul[3]
    TC_A.V='-60*mvolt+5*mvolt*rand()'
    TC_B.V='-60*mvolt+5*mvolt*rand()'
    S_TC_A_RSg_LIP=generate_syn(TC_A,RS_gran_LIP,'Isyn_FEF','',0.05*msiemens * cm **-2,0.125*ms,1*ms,0*mV)    
    S_TC_A_FSg_LIP=generate_syn(TC_A,FS_gran_LIP,'Isyn_FEF','',0.05*msiemens * cm **-2,0.125*ms,1*ms,0*mV)    
    S_TC_A_VIP_FEF=generate_syn(TC_A,VIP_FEF,'Isyn_mdPul','',0.05*msiemens * cm **-2,0.125*ms,1*ms,0*mV)    


#    S_LIP_RS_inter_mdPul=generate_syn(RS_sup_LIP,TC_inter,'Isyn_LIP','',0.05*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
#    S_LIP_RS_inter_mdPul=generate_syn(RS_sup_LIP,TC_inter,'Isyn_LIP','',0.025*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_LIP_RS_inter_mdPul=generate_syn(IB_LIP,TC_inter,'Isyn_LIP','',0.025*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
#    S_RS_inter_RSv_FEF=generate_syn(RS_inter,RSv_FEF,'Isyn_LIP','',0.010*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
#    S_RS_inter_SIv_FEF=generate_syn(RS_inter,SIv_FEF,'Isyn_LIP','',0.020*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
#    S_RS_inter_VIPv_FEF=generate_syn(RS_inter,VIPv_FEF,'Isyn_LIP','',0.005*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_RS_inter_RSv_FEF=generate_syn(TC_inter,RSv_FEF,'Isyn_LIP','',0.05*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    S_RS_inter_SIv_FEF=generate_syn(TC_inter,SIv_FEF,'Isyn_LIP','',0.05*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   
    S_RS_inter_VIPv_FEF=generate_syn(TC_inter,VIPv_FEF,'Isyn_LIP','',0.005*msiemens * cm **-2,0.125*ms,1*ms,0*mV)

    gGABA_A_HTC_TC_inter = 0.2 * msiemens * cm **-2
    
    synHTCTCinter=generate_syn_delay(HTC_A,TC_inter,'IsynHTC','',gGABA_A_HTC_TC_inter,0.25*ms,5*ms,-80*mV)

    S_RS_decision_IB_SC=generate_syn(RS_decision,IB_SC_bd,'Isyn_FEF','',0.5*msiemens * cm **-2,0.125*ms,1*ms,0*mV)   


    all_neurons=all_neurons_FEF+all_neurons_LIP+all_neurons_mdPul+(TC_inter,IB_SC_soma,IB_SC_axon,IB_SC_ad,IB_SC_bd)
#    print(len(all_neurons_FEF),len(all_neurons_LIP))
    all_synapses=all_synapses_FEF+all_synapses_LIP+all_synapses_mdPul+(S_FEF_IB_LIP,S_FEF_SIdeep_LIP,S_LIP_RS_FEF,S_LIP_FS_FEF,S_LIP_RS_inter_mdPul,S_RS_inter_RSv_FEF,S_RS_inter_SIv_FEF,S_RS_inter_VIPv_FEF)+(S_TC_A_RSg_LIP,S_TC_A_FSg_LIP,S_TC_A_VIP_FEF,synHTCTCinter,S_RS_decision_IB_SC,S_IBIB)
    all_monitors=all_monitors_FEF+all_monitors_LIP+all_monitors_mdPul+(TCinter,IBSC)
    all_gap_junctions=all_gap_junctions_LIP,all_gap_junctions_mdPul+(gapIB_SomaAd,gapIB_SomaBd,gapIB_SomaAxon,gapIB_AdSoma,gapIB_BdSoma,gapIB_AxonSoma,gap_IBIB)
    
    return all_neurons,all_synapses,all_gap_junctions,all_monitors
#    return all_neurons_mdPul,all_synapses_mdPul,all_gap_junctions_mdPul,all_monitors_mdPul

if __name__=='__main__':
    close('all')
    prefs.codegen.target = 'numpy'
    
    theta_phase='mixed'
    target_on=True  
#    target_on=False 
    target_time = 500*msecond
#    target_time = 625*msecond
#    target_time = 575*msecond #good to bad
#    target_time = 450*msecond #bad to good
#    target_time = 2000*msecond
    
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
#    in_mode='single_spike'
    in_mode='burst'
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

    all_neurons,all_synapses,all_gap_junctions,all_monitors=create_FEF_and_LIP_and_mdPul(N_RS_vis,N_FS_vis,N_RS_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime,target_time, syn_cond,J,thal, N_HTC,N_TC,N_RE)     
    
#    HTC_A,HTC_B=all_neurons[0],all_neurons[2]
    HTC_A,HTC_B=all_neurons[26],all_neurons[28]
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

    condition='mAChR'
#    in_mode='single_spike'
    in_mode='burst'
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
    
    for neur in all_neurons:
        net.add(neur)
#    net.add(all_neurons)
    for syn in all_synapses:
        net.add(syn)
#    net.add(all_synapses)
    for gj in all_gap_junctions:
        net.add(gj)
#    net.add(all_gap_junctions)
    for moni in all_monitors:
        net.add(moni)
#    net.add(all_monitors)
    net.run(runtime,report='text',report_period=300*second)
    
    R8,R9,R10,V_RS,V_FS,V_SI,R11,R12,R13,R14,mon_RS, V1,V2,V3,R1,R2,R3,I1,I2,I3,V4,R4,I4s,I4a,I4ad,I4bd,R5,R6,R7,V5,V6,V7,inpmon,inpIBmon,R1A,R2A,V1A,V2A,I1A,I2A,I3A,R1B,R2B,V1B,V2B,I1B,I2B,I3B,RA,RB,Rinter,IBSC=all_monitors
#    R1A,R2A,V1A,V2A,I1A,I2A,I3A,R1B,R2B,V1B,V2B,I1B,I2B,I3B,RA,RB=all_monitors
    # LIP Plots
    figure(figsize=(6,9))
#    subplot(411)
    up=100
    plot(R1.t,R1.i+140+up,'r.',label='RS cells')
    plot(R2.t,R2.i+120+up,'b.',label='FS cells')
    plot(R3.t,R3.i+100+up,'g.',label='SI cells')
    plot([0.2,runtime/second],[95+up,95+up],'k--')
    plot(R5.t,R5.i+70+up,'r.',label='Granular RS')
    plot(R6.t,R6.i+50+up,'b.',label='Granular FS')
    plot([0.2,runtime/second],[45+up,45+up],'k--')
    plot(R4.t,R4.i+20+up,'m.',label='IB cells')
    plot(R7.t,R7.i+up,'g.',label='Deep SI')
    xlim(0.2,runtime/second)
    plot([0.2,runtime/second],[up-10,up-10],'k')

#    subplot(412)
    up=0
#    title('FEF Visual-Motor Neurons')
    plot(R10.t,R10.i+0+up,'k.',label='VIP')
    plot(R8.t,R8.i+60+up,'r.',label='RS')
    plot(R9.t,R9.i+40+up,'g.',label='SOM')
    xlim(0.2,runtime/second)
#    plot([0.2,runtime/second],[up-10,up-10],'k')
#    xticks([],[])
#    legend(loc='upper left') 
    #xlabel('Time (s)')
#    ylabel('Neuron index')
    xlabel('Time (s)')
    ylabel('Neuron index')

    #FEF Plots    
    figure(figsize=(10,4))
    subplot(131)
    title('Visual Neurons')
    plot(R11.t,R11.i+0,'r.',label='RS')
    plot(R12.t,R12.i+20,'b.',label='FS')
    plot(R13.t,R13.i+40,'k.',label='VIP')
    plot(R14.t,R14.i+60,'g.',label='SOM')
    xlim(0.2,runtime/second)
    legend(loc='upper left')   
    xlabel('Time (s)')
    ylabel('Neuron index')
    
    subplot(132)
    title('Visual-Motor Neurons')
    plot(R10.t,R10.i+0,'g.',label='VIP')
    plot(R8.t,R8.i+60,'r.',label='RS')
    plot(R9.t,R9.i+40,'.',label='SOM',color='lime')
    xlim(0.2,runtime/second)
    legend(loc='upper left') 
    xlabel('Time (s)')
    ylabel('Neuron index')
    
    subplot(133)
    title('Decision cells')
    plot(mon_RS.t,mon_RS.i+0,'r.',label='RS')
    xlim(0.2,runtime/second)
    legend(loc='upper left') 
    xlabel('Time (s)')
    ylabel('Neuron index')
    
    tight_layout()
    
    
    figure(figsize=(5,5))
#    subplot(411)
    up=0
    plot(R1.t,R1.i+140+up,'r.',label='RS cells')
    plot(R2.t,R2.i+120+up,'b.',label='FS cells')
    plot(R3.t,R3.i+100+up,'g.',label='SI cells')
    plot([0.2,runtime/second],[95+up,95+up],'k--')
    plot(R5.t,R5.i+70+up,'r.',label='Granular RS')
    plot(R6.t,R6.i+50+up,'b.',label='Granular FS')
    plot([0.2,runtime/second],[45+up,45+up],'k--')
    plot(R4.t,R4.i+20+up,'m.',label='IB cells')
    plot(R7.t,R7.i+up,'g.',label='Deep SI')
    xlim(0.2,runtime/second)

    xlabel('Time (s)')
    ylabel('Neuron index')
    
    figure()
    plot(R1A.t,R1A.i+0,'r.',label='HTC')
    plot(R2A.t,R2A.i+20,'y.',label='TC')
    plot(RA.t,RA.i+100,'b.',label='RE lat')
    
    plot([0,1],[225,225],'k')

    plot(R1B.t,R1B.i+250,'r.')
    plot(R2B.t,R2B.i+270,'y.')
    plot(RB.t,RB.i+350,'b.')
    
    xlim(0,runtime/second)  
    yticks([150,500],['Object I','Object II'])
    legend()  
    
    figure()
    plot(V1A.t,V1A.V[0],label='HTC V')
    plot(V2A.t,V2A.V[0],label='TC V')
#    plot(M.t,M.V[0],label='RE V')
    legend()
    
    figure()
    plot(Rinter.t,Rinter.i+0,'r.',label='Rinter')
    xlim(0,runtime/second)
    legend()  
    
    figure()
    plot(IBSC.t,IBSC.i+0,'r.',label='IB_SC')
    xlim(0,runtime/second)
    legend()  
    
    clear_cache('cython')
