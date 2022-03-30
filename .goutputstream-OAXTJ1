#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 11:24:49 2020

@author: amelie
"""

from brian2 import *

from scipy import signal
from cells.RS_FEF import *
from cells.FS_FEF import *
from cells.SI_FEF import *
from cells.VIP_FEF import *

from FEF_visuomotor import *

runtime=3*second
    
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


def generate_spike_timing(N,f,start_time,end_time=runtime):
    list_time_and_i=[]
    for i in range(N):
        list_time=[(start_time,i)]
        next_spike=list_time[-1][0]+(1+0.01*rand())/f
        while next_spike<end_time:
            list_time.append((next_spike,i))
            next_spike=list_time[-1][0]+(1+0.01*rand())/f
        list_time_and_i+=list_time
    return array(list_time_and_i)


    
def create_visual_neurons(N_RS,N_FS,target_on=False):
    ##Define neuron groups
    RS_vis=NeuronGroup(N_RS,eq_RS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    RS_vis.V = '-70*mvolt+10*rand()*mvolt'
    RS_vis.h = '0+0.05*rand()'
    RS_vis.m = '0+0.05*rand()'
    RS_vis.mAR = '0.035+0.025*rand()'
    RS_vis.J='30 * uA * cmeter ** -2'

    FS_vis=NeuronGroup(N_FS,eq_FS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    FS_vis.V = '-110*mvolt+10*rand()*mvolt'
    FS_vis.h = '0+0.05*rand()'
    FS_vis.m = '0+0.05*rand()'
    FS_vis.J='10 * uA * cmeter ** -2'

    ##Synapses    
    S_RSvRSv=generate_syn(RS_vis,RS_vis,'IsynRS_FEF_V','',0.3*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_RSvFSv=generate_syn(RS_vis,FS_vis,'IsynRS_FEF_V','',0.6*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_FSvRSv=generate_syn(FS_vis,RS_vis,'IsynFS_FEF_V','',0.6* msiemens * cm **-2,0.25*ms,5*ms,-80*mV)
    S_FSvFSv=generate_syn(FS_vis,FS_vis,'IsynFS_FEF_V','',0.3*msiemens * cm **-2,0.25*ms,5*ms,-75*mV)
   
    ##Input : gamma when target appears
    if target_on:    
        FS_vis.ginp_FS=5*msiemens * cm **-2
        RS_vis.ginp_RS=5*msiemens * cm **-2
    else:    
        FS_vis.ginp_FS=0*msiemens * cm **-2
        RS_vis.ginp_RS=0*msiemens * cm **-2
        
    fvis=50*Hz
    visual=generate_spike_timing(N_FS,fvis,1500*ms,end_time=3000*ms)
    visual_input = SpikeGeneratorGroup(N_FS, visual[:,1], visual[:,0]*second)
    visual_in=Synapses(visual_input,FS_vis,on_pre='Vinp=Vhigh')
    visual_in.connect(j='i')
    visual_in2=Synapses(visual_input,RS_vis,on_pre='Vinp=Vhigh')
    visual_in2.connect(j='i')    
        
    ##Monitors
    R1=SpikeMonitor(RS_vis,record=True)
    R2=SpikeMonitor(FS_vis,record=True)
    V1=StateMonitor(RS_vis,'V',record=True)
    V2=StateMonitor(FS_vis,'V',record=True)
    
    all_neurons=RS_vis,FS_vis,visual_input
    all_synapses=S_RSvRSv,S_RSvFSv,S_FSvRSv,S_FSvFSv,visual_in,visual_in2
    all_monitors=R1,R2,V1,V2
    
    return all_neurons,all_synapses,all_monitors

#
#def create_motor_neurons(N_RS,N_SI):
#    RS_mot=NeuronGroup(N_RS,eq_RS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
#    RS_mot.V = '-70*mvolt+10*rand()*mvolt'
#    RS_mot.h = '0+0.05*rand()'
#    RS_mot.m = '0+0.05*rand()'
#    RS_mot.mAR = '0.035+0.025*rand()'
#    RS_mot.J='20 * uA * cmeter ** -2'   
#    
#    SI_mot=NeuronGroup(N_SI,eq_SI_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
#    SI_mot.V = '-110*mvolt+10*rand()*mvolt'
#    SI_mot.h = '0+0.05*rand()'
#    SI_mot.m = '0+0.05*rand()'
#    SI_mot.J='5 * uA * cmeter ** -2' #article=code=35
##    SI_mot.ginpSIdeep='0 * msiemens * meter **-2'
#    
#    Fix=NeuronGroup(N_SI,eq_SI_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
#    Fix.V = '-110*mvolt+10*rand()*mvolt'
#    Fix.h = '0+0.05*rand()'
#    Fix.m = '0+0.05*rand()'
#    Fix.J='-20 * uA * cmeter ** -2' #article=code=35
#    
#    S_RSmRSm=generate_syn(RS_mot,RS_mot,'IsynRS','',0.3*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
#    S_RSmSIm=generate_syn(RS_mot,SI_mot,'IsynRS','',0.8*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
#    S_SImRSm=generate_syn(SI_mot,RS_mot,'IsynSIdeep','',0.8* msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
#    S_SImSIm=generate_syn(SI_mot,SI_mot,'IsynSIdeep','',0.3*msiemens * cm **-2,0.25*ms,20*ms,-75*mV)
##    S_SImRSm=generate_syn(SI_mot,RS_mot,'IsynSIdeep','',0* msiemens * cm **-2,0.25*ms,5*ms,-80*mV)
##    S_SImSIm=generate_syn(SI_mot,SI_mot,'IsynSIdeep','',0*msiemens * cm **-2,0.25*ms,5*ms,-75*mV)
#    S_FixRSm=generate_syn(Fix,RS_mot,'IsynFSgran','',1.5* msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
#    S_FixSIm=generate_syn(Fix,SI_mot,'IsynFSgran','',0* msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
#    S_SImFix=generate_syn(SI_mot,Fix,'IsynSIdeep','',1.5* msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
#    
#    R1=SpikeMonitor(RS_mot,record=True)
#    V1=StateMonitor(RS_mot,'V',record=True)
#    R2=SpikeMonitor(SI_mot,record=True)
#    V2=StateMonitor(SI_mot,'V',record=True)
#    R3=SpikeMonitor(Fix,record=True)
#    
#    all_neurons=RS_mot,SI_mot,Fix
#    all_synapses=S_RSmRSm,S_RSmSIm,S_SImRSm,S_SImSIm,S_FixRSm,S_FixSIm,S_SImFix
#    all_monitors=R1,R2,V1,V2,R3
#    
#    return all_neurons,all_synapses,all_monitors

#def create_network(N_RS_vis,N_FS_vis,N_RS_mot,N_SI_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime):
#    start_scope()
#    net=Network()
#    
#    #create each functional group of neurons individually
#    all_neurons_vm,all_synapses_vm,all_monitors_vm=generate_deepSI_and_gran_layers(theta_phase,N_dSI_vm,N_RS_vm,N_gSI_vm,runtime)
#    RS_vm,gSI_vm=all_neurons_vm[1],all_neurons_vm[2]
#    
#    all_neurons_v,all_synapses_v,all_monitors_v=create_visual_neurons(N_RS_vis,N_FS_vis,target_on)
#    RS_vis,FS_vis=all_neurons_v[0],all_neurons_v[1]
#    print(FS_vis.J)
#    
#    all_neurons_m,all_synapses_m,all_monitors_m=create_motor_neurons(N_RS_mot,N_SI_mot)
#    RS_mot,SI_mot=all_neurons_m[0],all_neurons_m[1]
#    
#    #other column's VM neurons:
#    if theta_phase=='good':
#        all_neurons_vm2,all_synapses_vm2,all_monitors_vm2=generate_deepSI_and_gran_layers('bad',N_dSI_vm,N_RS_vm,N_gSI_vm,runtime)
#        RS_vm2,gSI_vm2=all_neurons_vm2[1],all_neurons_vm2[2]
#    else :
#        all_neurons_vm2,all_synapses_vm2,all_monitors_vm2=generate_deepSI_and_gran_layers('good',N_dSI_vm,N_RS_vm,N_gSI_vm,runtime)
#        RS_vm2,gSI_vm2=all_neurons_vm2[1],all_neurons_vm2[2]
#    
#    #add the synapses between groups
#    #From visual-motor to visual
##    S_gSIvmRSv=generate_syn(gSI_vm,RS_vis,'IsynSI','',1*msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
#    S_gSIvmFSv=generate_syn(gSI_vm,FS_vis,'IsynSI','',0.5*msiemens * cm **-2,0.25*ms,20*ms,-80*mV)    
#    
#    #From visual to visual-motor
#    S_RSvRSvm=generate_syn(RS_vis,RS_vm,'IsynRS','',0*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
#    S_RSvgSIvm=generate_syn(RS_vis,gSI_vm,'IsynRS','',0*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
#
#    #From visual-motor to motor       
##    S_RSvmRSm=generate_syn(RS_vm,RS_mot,'IsynEgran','',1*msiemens * cm **-2,0.125*ms,1*ms,0*mV) 
##    S_RSvmSIm=generate_syn(RS_vm,SI_mot,'IsynEgran','',1*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
#
#    S_SIvmSIm=generate_syn(gSI_vm,SI_mot,'IsynSI','',0.3*msiemens * cm **-2,0.25*ms,20*ms,-80*mV)  
##    S_SIvmSIm=generate_syn(RS_vm,SI_mot,'IsynEgran','',0.5*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
##    S_SIvmRSm=generate_syn(RS_vm,RS_mot,'IsynEgran','',0.5*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
#
#    S_SIvm2SIm=generate_syn(RS_vm2,SI_mot,'IsynIB','',0.3*msiemens * cm **-2,0.25*ms,20*ms,-80*mV) 
##    S_SIvm2SIm=generate_syn(RS_vm2,SI_mot,'IsynIB','',0.5*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
##    S_SIvm2RSm=generate_syn(RS_vm2,RS_mot,'IsynIB','',0.5*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
#
#
#    #From visual to motor
#    S_RSvRSm=generate_syn(RS_vis,RS_mot,'IsynEgran','',0.5*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
#    S_RSvSIm=generate_syn(RS_vis,SI_mot,'IsynEgran','',0.5*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
##    S_RSvSIm=generate_syn(FS_vis,SI_mot,'IsynFS','',1*msiemens * cm **-2,0.25*ms,5*ms,-80*mV)
#        
#     
#    for elem in [all_neurons_vm,all_synapses_vm,all_monitors_vm,all_neurons_v,all_synapses_v,all_monitors_v,all_neurons_m,all_synapses_m,all_monitors_m,all_neurons_vm2,all_synapses_vm2]:
#        net.add(elem)
#
#        
##    for elem in [S_gSIvmRSv,S_gSIvmFSv,S_RSvRSvm,S_RSvgSIvm,S_RSvmRSm]:
##    for elem in [S_gSIvmFSv,S_RSvRSvm,S_RSvgSIvm,S_RSvmRSm,S_RSvmSIm]:
#    for elem in [S_gSIvmFSv,S_RSvRSvm,S_RSvgSIvm,S_SIvmSIm,S_RSvRSm,S_RSvSIm,S_SIvm2SIm]:
##    for elem in [S_gSIvmFSv,S_RSvRSvm,S_RSvgSIvm,S_SIvmSIm,S_SIvmRSm,S_RSvSIm,S_SIvm2SIm,S_SIvm2RSm]:
#        net.add(elem)
#        
#    mon_RS=StateMonitor(RS_vis,'Isyn',record=True)
#    net.add(mon_RS)
#        
#    all_monitors=all_monitors_vm+all_monitors_v+all_monitors_m+(mon_RS,)
#    
#    return net,all_monitors

#def create_network_no_motor(N_RS_vis,N_FS_vis,N_RS_mot,N_SI_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime):
#    start_scope()
#    net=Network()
#    
#    #create each functional group of neurons individually
#    all_neurons_vm,all_synapses_vm,all_monitors_vm=generate_deepSI_and_gran_layers(theta_phase,N_dSI_vm,N_RS_vm,N_gSI_vm,runtime)
#    RS_vm,gSI_vm=all_neurons_vm[1],all_neurons_vm[2]
#    
#    all_neurons_v,all_synapses_v,all_monitors_v=create_visual_neurons(N_RS_vis,N_FS_vis,target_on)
#    RS_vis,FS_vis=all_neurons_v[0],all_neurons_v[1]
#    
#    
#    #add the synapses between groups
#    #From visual-motor to visual
##    S_gSIvmRSv=generate_syn(gSI_vm,RS_vis,'IsynSI','',1*msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
#    S_gSIvmFSv=generate_syn(gSI_vm,FS_vis,'IsynSI_FEF_VM','',1*msiemens * cm **-2,0.25*ms,20*ms,-80*mV)    
#    
#    #From visual to visual-motor
#    S_RSvRSvm=generate_syn(RS_vis,RS_vm,'IsynRS_FEF_V','',0*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
#    S_RSvgSIvm=generate_syn(RS_vis,gSI_vm,'IsynRS_FEF_V','',0*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
#
#
#     
#    for elem in [all_neurons_vm,all_synapses_vm,all_monitors_vm,all_neurons_v,all_synapses_v,all_monitors_v]:
#        net.add(elem)
#
#        
##    for elem in [S_gSIvmRSv,S_gSIvmFSv,S_RSvRSvm,S_RSvgSIvm,S_RSvmRSm]:
##    for elem in [S_gSIvmFSv,S_RSvRSvm,S_RSvgSIvm,S_RSvmRSm,S_RSvmSIm]:
#    for elem in [S_gSIvmFSv,S_RSvRSvm,S_RSvgSIvm]:
#        net.add(elem)
#        
#    mon_FS=StateMonitor(FS_vis,'Isyn',record=True)
#    net.add(mon_FS)
#        
#    all_monitors=all_monitors_vm+all_monitors_v+(mon_FS,)
#    
#    return net,all_monitors


def create_network_no_motor2(N_RS_vis,N_FS_vis,N_RS_mot,N_SI_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime):

    
    #create each functional group of neurons individually
    all_neurons_vm,all_synapses_vm,all_monitors_vm=generate_deepSI_and_gran_layers(theta_phase,N_dSI_vm,N_RS_vm,N_gSI_vm,runtime)
    RS_vm,gSI_vm=all_neurons_vm[1],all_neurons_vm[2]
    
    all_neurons_v,all_synapses_v,all_monitors_v=create_visual_neurons(N_RS_vis,N_FS_vis,target_on)
    RS_vis,FS_vis=all_neurons_v[0],all_neurons_v[1]
    
    
    #add the synapses between groups
    #From visual-motor to visual
#    S_gSIvmRSv=generate_syn(gSI_vm,RS_vis,'IsynSI','',1*msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
    S_gSIvmFSv=generate_syn(gSI_vm,FS_vis,'IsynSI_FEF_VM','',1*msiemens * cm **-2,0.25*ms,20*ms,-80*mV)    
    
    #From visual to visual-motor
    S_RSvRSvm=generate_syn(RS_vis,RS_vm,'IsynRS_FEF_V','',0*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_RSvgSIvm=generate_syn(RS_vis,gSI_vm,'IsynRS_FEF_V','',0*msiemens * cm **-2,0.125*ms,1*ms,0*mV)


     
#    for elem in [all_neurons_vm,all_synapses_vm,all_monitors_vm,all_neurons_v,all_synapses_v,all_monitors_v]:
#        net.add(elem)
#
#        
##    for elem in [S_gSIvmRSv,S_gSIvmFSv,S_RSvRSvm,S_RSvgSIvm,S_RSvmRSm]:
##    for elem in [S_gSIvmFSv,S_RSvRSvm,S_RSvgSIvm,S_RSvmRSm,S_RSvmSIm]:
#    for elem in [S_gSIvmFSv,S_RSvRSvm,S_RSvgSIvm]:
#        net.add(elem)
       
    mon_FS=StateMonitor(FS_vis,'Isyn',record=True)
#    net.add(mon_FS)
        
    all_monitors=all_monitors_vm+all_monitors_v+(mon_FS,)
    all_neurons=all_neurons_vm+all_neurons_v
    all_synapses=all_synapses_vm+all_synapses_v+(S_gSIvmFSv,S_RSvRSvm,S_RSvgSIvm)
    
    
    return all_neurons,all_synapses,all_monitors




if __name__=='__main__':
    close('all')
    prefs.codegen.target = 'numpy'
    defaultclock.dt = 0.01*ms
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp=0* msiemens * cm **-2
    ginp_SI=0* msiemens * cm **-2
    
    print('Creating the network')
    N_RS_vis,N_FS_vis,N_RS_mot,N_SI_mot,N_dSI_vm,N_RS_vm,N_gSI_vm=[20]*7
    
    theta_phase='mixed'
    target_on=True
    runtime=1*second
    
    net=Network()
    
#    net,all_monitors=create_network(N_RS_vis,N_FS_vis,N_RS_mot,N_SI_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime)
    all_neurons,all_synapses,all_monitors=create_network_no_motor2(N_RS_vis,N_FS_vis,N_RS_mot,N_SI_mot,N_dSI_vm,N_RS_vm,N_gSI_vm,theta_phase,target_on,runtime)
    
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_monitors)
    
    print('Compiling with cython')
    prefs.codegen.target = 'cython'
    net.run(runtime,report='text',report_period=300*second)
    
#    R1,R2,R3,V1,V2,V3,R4,R5,V4,V5,R6,R7,V6,V7,R8,mon_RS=all_monitors
    R1,R2,R3,V1,V2,V3,R4,R5,V4,V5,mon_FS=all_monitors
    
    figure(figsize=(10,4))
    subplot(131)
    title('Visual Neurons')
    plot(R4.t,R4.i+20,'r.',label='RS')
    plot(R5.t,R5.i+0,'b.',label='FS')
    xlim(0,runtime/second)
    legend(loc='upper left')   
    
    subplot(132)
    title('Visual-Motor Neurons')
    plot(R3.t,R3.i+0,'c.',label='deep SI')
    plot(R1.t,R1.i+60,'r.',label='gran RS')
    plot(R2.t,R2.i+40,'b.',label='gran SI')
    xlim(0,runtime/second)
    legend(loc='upper left') 
    
#    subplot(133)
#    title('Motor Neurons')
#    plot(R6.t,R6.i+60,'r.',label='RS')
#    plot(R7.t,R7.i+40,'b.',label='SI')
#    plot(R8.t,R8.i+0,'c.',label='Fix')
#    xlim(0,runtime/second)
#    legend(loc='upper left') 
    
    
    min_t=int(50*ms*100000*Hz)
    LFP_V1=1/20*sum(V1.V,axis=0)[min_t:]
    LFP_V2=1/20*sum(V2.V,axis=0)[min_t:]
    LFP_V3=1/20*sum(V3.V,axis=0)[min_t:]
    LFP_V4=1/20*sum(V4.V,axis=0)[min_t:]
    LFP_V5=1/20*sum(V5.V,axis=0)[min_t:]
#    LFP_V6=1/20*sum(V6.V,axis=0)[min_t:]
#    LFP_V7=1/20*sum(V7.V,axis=0)[min_t:]
    
    f,Spectrum_LFP_V1=signal.periodogram(LFP_V1, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V2=signal.periodogram(LFP_V2, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V3=signal.periodogram(LFP_V3, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V4=signal.periodogram(LFP_V4, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V5=signal.periodogram(LFP_V5, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V6=signal.periodogram(LFP_V6, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V7=signal.periodogram(LFP_V7, 100000,'flattop', scaling='spectrum')

    figure(figsize=(10,4))
    subplot(331)
    plot(f,Spectrum_LFP_V4)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,100)
    title('visual RS')
    subplot(334)
    plot(f,Spectrum_LFP_V5)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,100)
    title('visual FS')  
    
    subplot(332)
    plot(f,Spectrum_LFP_V1)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,100)
    title('visual-motor gran RS')
    subplot(335)
    plot(f,Spectrum_LFP_V2)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,100)
    title('visual-motor gran SI')  
    subplot(338)
    plot(f,Spectrum_LFP_V3)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,100)
    title('visual-motor deep SI')  
    
#    subplot(333)
#    plot(f,Spectrum_LFP_V6)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('motor RS')
#    subplot(336)
#    plot(f,Spectrum_LFP_V7)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('motor SI')
#    tight_layout()
    
#    figure()
#    plot(mon_RS.t,mon_RS.Isyn[0])
#    plot(mon_RS.t,mon_RS.Isyn[10])
    
    clear_cache('cython')
    