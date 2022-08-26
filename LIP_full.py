# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 14:28:44 2019

@author: aaussel
"""

#This is used for simulations that change the time constants of all interneurons of one type at the same time.

from brian2 import *

start_scope()

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
from LIP_beta1 import *

import os

os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_DYNAMIC'] = 'FALSE'

import time
from itertools import *

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
    
def make_full_network(syn_cond,J,thal,t_SI,t_FS,theta_phase,target_time):  
    
    NN=1 #multiplicative factor on the number of neurons
    N_RS,N_FS,N_SI,N_IB= NN*80,NN*20,NN*20,NN*20 #Number of neurons of RE, TC, and HTC type
    
    gSIdFSg,gFSgRSg,gRSgFSg,gRSgRSg,gFSgFSg,gRSgRSs,gRSgFSs,gFSgRSs=syn_cond
    J_RSg,J_FSg=J
    
    runtime=3*second
    
    all_neurons, all_synapses, all_gap_junctions, all_monitors=create_Mark_Alex_network(kainate,version,Nf=NN)
    V1,V2,V3,V4,R1,R2,R3,R4,I1,I2,I3,I4,V5,R5,I5s,I5a,I5ad,I5IBbd=all_monitors
    RS, FS, SI, VIP, IB_soma, IB_axon, IB_bd, IB_ad =all_neurons

    prefs.codegen.target = 'numpy'
    defaultclock.dt = 0.01*ms
    
    #Single column network
    
    ##Define neuron groups
    E_gran=NeuronGroup(N_FS,eq_RS_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    E_gran.V = '-70*mvolt+10*rand()*mvolt'
    E_gran.h = '0+0.05*rand()'
    E_gran.m = '0+0.05*rand()'
    E_gran.mAR = '0.035+0.025*rand()'
    E_gran.J=J_RSg
    
    FS_gran=NeuronGroup(N_FS,eq_FS_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    FS_gran.V = '-110*mvolt+10*rand()*mvolt'
    FS_gran.h = '0+0.05*rand()'
    FS_gran.m = '0+0.05*rand()'
    FS_gran.J=J_FSg
    
    SI_deep=NeuronGroup(N_SI,eq_SI_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    SI_deep.V = '-100*mvolt+10*rand()*mvolt'
    SI_deep.h = '0+0.05*rand()'
    SI_deep.m = '0+0.05*rand()'
    SI_deep.mAR = '0.02+0.04*rand()'
    SI_deep.J='35* uA * cmeter ** -2' #article SI=50, code=35, Mark = 45
 
    if theta_phase=='bad':
        SI_deep.ginp_SI=0* msiemens * cm **-2 #FEF input to deep SOM cells is zero
        mdpul_input_amplitude=0* msiemens * cm **-2 #mdPul input to granular layer is zero
        
    if theta_phase=='good':
        SI_deep.ginp_SI=5* msiemens * cm **-2 #FEF input to deep SOM cells
        mdpul_input_amplitude=thal #mdPul input to granular layer
        
    if theta_phase=='mixed':
        SI_deep.ginp_SI=5* msiemens * cm **-2 #FEF input to deep SOM cells
        mdpul_input_amplitude=thal #mdPul input to granular layer
            
    
    
    Vlow=-80*mV
    SI_deep.Vinp=Vlow
        
    ##Synapses
    eq_syn='''_post=s_i*g_i*(V_post-V_i) : amp * meter ** -2 (summed)
        ds_i/dt=-s_i/taud_i+(1-s_i)/taur_i*0.5*(1+tanh(V_pre/10/mV)) : 1
        g_i : siemens * meter**-2
        V_i : volt
        taud_i : second
        taur_i : second
    '''
    
    def generate_syn(source,target,syntype,connection_pattern,g_i,taur_i,taud_i,V_i):
        S=Synapses(source,target,model=syntype+eq_syn,method='exact')
        if connection_pattern=='':
            S.connect()
        else :
            S.connect(condition=connection_pattern, skip_if_invalid=True)
        S.g_i=g_i
        S.taur_i=taur_i
        S.taud_i=taud_i
        S.V_i=V_i  
        return S
    
    #From E (granular layer) cells
    S_EgranFS=generate_syn(E_gran,FS,'IsynRS_LIP_gran','i//10==j//10',2*gRSgFSs,0.125*ms,1*ms,0*mV)
    S_EgranEgran=generate_syn(E_gran,E_gran,'IsynRS_LIP_gran','i//10==j//10',2*gRSgRSg,0.125*ms,1*ms,0*mV)
    S_EgranFSgran=generate_syn(E_gran,FS_gran,'IsynRS_LIP_gran','i//10==j//10',2*gRSgFSg,0.125*ms,1*ms,0*mV)
    S_EgranRS=generate_syn(E_gran,RS,'IsynRS_LIP_gran','i//10==j//40',2*gRSgRSs,0.125*ms,1*ms,0*mV)
    S_EgranIB=generate_syn(E_gran,IB_ad,'IsynRS_LIP_gran','i//10==j//10',2*0.212*usiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_EgranVIP=generate_syn(E_gran,VIP,'IsynRS_LIP_gran','i//10==j//10',2*gRSgFSg,0.125*ms,1*ms,0*mV)
    
    #From FS (granular layer) cells
    S_FSgranEgran=generate_syn(FS_gran,E_gran,'IsynFS_LIP_gran','i//10==j//10',2*gFSgRSg,0.25*ms,t_FS,-80*mV)
    S_FSgranFSgran=generate_syn(FS_gran,FS_gran,'IsynFS_LIP_gran','i//10==j//10',2*gFSgFSg,0.25*ms,t_FS,-75*mV)
    S_FSgranRS=generate_syn(FS_gran,RS,'IsynFS_LIP_gran','i//10==j//40',2*gFSgRSs,0.25*ms,t_FS,-80*mV)
    
    #From IB cells
    S_IBSIdeep=generate_syn(IB_axon,SI_deep,'IsynIB_LIP','i//10==j//10',2*0.01* msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    
    #From deep SOM cells    
    S_SIdeepIB=generate_syn(SI_deep,IB_bd,'IsynSI_LIP_deep','i//10==j//10',2*10* msiemens * cm **-2,0.25*ms,t_SI,-80*mV)
    S_SIdeepFSgran=generate_syn(SI_deep,FS_gran,'IsynSI_LIP_deep','i//10==j//10',2*gSIdFSg,0.25*ms,t_SI,-80*mV)
    
    
    def generate_spike_timing(N,f,start_time,end_time=runtime):
        list_time_and_i=[]
        for i in range(N):
            list_time=[(start_time,i)]
            next_spike=list_time[-1][0]+(1+0.001*rand())/f
            while next_spike<end_time:
                list_time.append((next_spike,i))
                next_spike=list_time[-1][0]+(1+0.001*rand())/f
            list_time_and_i+=list_time
        return array(list_time_and_i)
    
    def generate_spike_timing_theta(N,f,start_time,end_time=runtime,f_theta=4*Hz):
        list_time_and_i=[]
        for i in range(N):
            list_time=[(start_time,i)]
#            next_spike=list_time[-1][0]+(1+0.1*rand())/f
            next_spike=list_time[-1][0]+1/f
            while next_spike<end_time:
                if int(sin(2*pi*next_spike*f_theta)>0)==1:
                    list_time.append((next_spike,i))
                next_spike=next_spike+1/f
            list_time_and_i+=list_time
        return array(list_time_and_i)
    
    G_topdown,G_topdown2,G_topdown3,G_lateral,G_lateral2,Poisson_input,Poisson_input2=[None]*7
    topdown_in,topdown_in2,topdown_in3,lateral_in,lateral_in2,bottomup_in,bottomup_in2=[None]*7


    #defining input to the granular layer (from mdpul)
    E_gran.ginp_RS_good=mdpul_input_amplitude
    FS_gran.ginp_FS_good=mdpul_input_amplitude
    E_gran.ginp_RS_bad=mdpul_input_amplitude
    FS_gran.ginp_FS_bad=mdpul_input_amplitude
    
    inputs_mdpul=generate_spike_timing(N_FS,13*Hz,0*ms,end_time=10000*ms)
    
    if target_on:
        # if theta_phase=='good':
        #     fFEF=25*Hz
        # else :
        #     fFEF=0*Hz
            
        # gamma_background=generate_spike_timing(N_FS,fFEF,0*ms,end_time=3000*ms)
        
        # if theta_phase=='mixed':
        #     t0=0*ms
        #     t1=125*ms
        #     gamma_background=generate_spike_timing(N_FS,fFEF,t0,end_time=t1)
        #     while t0+125*ms<runtime:
        #         fFEF=25*Hz*int(fFEF==0*Hz)+0*Hz*int(fFEF==25*Hz)
        #         t0,t1=t0+125*ms,t1+125*ms
        #         gamma_background=vstack((gamma_background,generate_spike_timing(N_FS,fLIP,t0,end_time=t1)))
            
        gamma_target=generate_spike_timing(10,50*Hz,target_time,end_time=target_time+100*ms)
        
        # Poisson_background = SpikeGeneratorGroup(N_FS, gamma_background[:,1], gamma_background[:,0]*second)
        Poisson_target = SpikeGeneratorGroup(10, gamma_target[:,1], gamma_target[:,0]*second)
        
        # S_in_bg_RS=Synapses(Poisson_background,RS,on_pre='Vinp=Vhigh')
        # S_in_bg_RS.connect(j='i')
        # S_in_bg_SI=Synapses(Poisson_background,SI,on_pre='Vinp=Vhigh')
        # S_in_bg_SI.connect(j='i')
        # S_in_bg_VIP=Synapses(Poisson_background,VIP,on_pre='Vinp=Vhigh')
        # S_in_bg_VIP.connect(j='i')
        
        S_in_target_VIP=Synapses(Poisson_target,VIP,on_pre='Vinp2=Vhigh')
        S_in_target_VIP.connect(j='i')
        S_in_target_SI=Synapses(Poisson_target,SI,on_pre='Vinp2=Vhigh')
        S_in_target_SI.connect(j='i')
        S_in_target_RS=Synapses(Poisson_target,RS,on_pre='Vinp2=Vhigh')
        S_in_target_RS.connect(j='i')
        S_in_target_Eg=Synapses(Poisson_target,E_gran,on_pre='Vinp2=Vhigh')
        S_in_target_Eg.connect(j='i')
        # S_in_target_FSg=Synapses(Poisson_target,FS_gran,on_pre='Vinp2=Vhigh')
        # S_in_target_FSg.connect(j='i')
        
        target_multiplier = 1
        SI.ginp_SI2=target_multiplier*2.5* msiemens * cm **-2
        VIP.ginp_VIP2=target_multiplier*3* msiemens * cm **-2
        RS.ginp_RS2=target_multiplier*2.5* msiemens * cm **-2
        E_gran.ginp_RS2=target_multiplier*2.5* msiemens * cm **-2

    if theta_phase=='mixed':
        t0=0*ms
        t1=125*ms
        inputs_mdpul=generate_spike_timing(N_FS,13*Hz,t0,end_time=t1)
        while t0+250*ms<runtime:
            t0,t1=t0+250*ms,t1+250*ms
            inputs_mdpul=vstack((inputs_mdpul,generate_spike_timing(N_SI,13*Hz,t0,end_time=t1)))
                          
        
    Poisson_input = SpikeGeneratorGroup(N_FS, inputs_mdpul[:,1], inputs_mdpul[:,0]*second)
    bottomup_in = Synapses(Poisson_input,FS_gran, on_pre='Vinp=Vhigh')
    bottomup_in.connect(j='i')

    Poisson_input2 = SpikeGeneratorGroup(N_FS, inputs_mdpul[:,1], inputs_mdpul[:,0]*second)
    bottomup_in2 = Synapses(Poisson_input2,E_gran, on_pre='Vinp=Vhigh')
    bottomup_in2.connect(j='i')
    
    #defining input to the deep layer layer (from FEFvm)
    if theta_phase=='good':
        inputs_topdown3=generate_spike_timing(N_SI,25*Hz,0*ms,end_time=5100*ms)
        G_topdown3 = SpikeGeneratorGroup(N_SI, inputs_topdown3[:,1], inputs_topdown3[:,0]*second)
        topdown_in3=Synapses(G_topdown3,SI_deep,on_pre='Vinp=Vhigh')
        topdown_in3.connect(j='i')
    
    if theta_phase=='mixed':
        inputs_topdown3=generate_spike_timing_theta(N_SI,25*Hz,0*ms,end_time=5100*ms)
        G_topdown3 = SpikeGeneratorGroup(N_SI, inputs_topdown3[:,1], inputs_topdown3[:,0]*second)
        topdown_in3=Synapses(G_topdown3,SI_deep,on_pre='Vinp=Vhigh')
        topdown_in3.connect(j='i')
    
    
    g_inputs=[G_topdown2,G_topdown3,G_lateral,G_lateral2,Poisson_input,Poisson_input2,Poisson_target]
    g_inputs=[y for y in g_inputs if y]
    syn_inputs=[topdown_in2,topdown_in3,lateral_in,lateral_in2,bottomup_in,bottomup_in2,S_in_target_VIP,S_in_target_SI,S_in_target_RS,S_in_target_Eg,S_EgranVIP]
    syn_inputs=[y for y in syn_inputs if y]
    

    
    #Define monitors and run network :
    R6=SpikeMonitor(E_gran,record=True)
    R7=SpikeMonitor(FS_gran,record=True)
    R8=SpikeMonitor(SI_deep,record=True)
    
    V6=StateMonitor(E_gran,'V',record=True)
    V7=StateMonitor(FS_gran,'V',record=True)
    V8=StateMonitor(SI_deep,'V',record=True)
    
#    inpmon=StateMonitor(E_gran,'Iinp1',record=True)
    inpmon=StateMonitor(E_gran,'sinp',record=True)
    #graninpmon=StateMonitor(FS,'IsynEgran',record=[0])
    inptarget=SpikeMonitor(Poisson_target,record=True)
    all_neurons=all_neurons+(E_gran,FS_gran,SI_deep)+tuple(g_inputs)
    all_synapses=all_synapses+(S_EgranFS,S_EgranEgran,S_EgranFSgran,S_EgranRS,S_EgranIB,S_FSgranEgran,S_FSgranFSgran,S_FSgranRS,S_IBSIdeep,S_SIdeepIB,S_SIdeepFSgran)+tuple(syn_inputs)
    all_monitors=all_monitors+(R6,R7,R8,V6,V7,V8,inpmon,inptarget)
    return all_neurons, all_synapses, all_gap_junctions, all_monitors


def run_one_LIP_simulation(simu,path,plot_raster=False):
#    print(simu,len(simu))
    start_scope()

    target_time,N_simu,t_SI,t_FS,theta_phase,g_LIP_FEF_v,target_on,runtime=simu[0],simu[1],simu[2],simu[3],simu[4],simu[5],simu[6],simu[7]

    if not plot_raster :
        new_path=path+"/results_"+str(N_simu)
        os.mkdir(new_path)
    else :
        new_path=path
        
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp_IB=0* msiemens * cm **-2
    ginp_SI=0* msiemens * cm **-2
    ginp=0* msiemens * cm **-2

    NN=1 #multiplicative factor on the number of neurons
    N_RS,N_FS,N_SI,N_IB= NN*80,NN*20,NN*20,NN*20 #Number of neurons of RE, TC, and HTC type
    
    
    all_SIdFSg=2*msiemens * cm **-2
    all_FSgRSg=1* msiemens * cm **-2
    all_RSgFSg=1*msiemens * cm **-2
    all_RSgRSg=0.5*msiemens * cm **-2
    all_FSgFSg=0.3* msiemens * cm **-2
    all_RSgRSs=2*msiemens * cm **-2
    all_RSgFSs=0.1*msiemens * cm **-2
    all_FSgRSs=0.1* msiemens * cm **-2
    all_J_RSg='10 * uA * cmeter ** -2'
    all_J_FSg='-5 * uA * cmeter ** -2'
    
    thal=10* msiemens * cm **-2
    
    syn_cond=(all_SIdFSg,all_FSgRSg,all_RSgFSg,all_RSgRSg,all_FSgFSg,all_RSgRSs,all_RSgFSs,all_FSgRSs)
    J=(all_J_RSg,all_J_FSg)
    
    net = Network(collect())
    
    print('Network setup')
    all_neurons, all_synapses, all_gap_junctions, all_monitors=make_full_network(syn_cond,J,thal,t_SI,t_FS,theta_phase)
    V1,V2,V3,R1,R2,R3,I1,I2,I3,V4,R4,I4s,I4a,I4ad,I4bd,R5,R6,R7,V5,V6,V7=all_monitors
    
    
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_gap_junctions)
#    net.add(all_monitors)
    net.add((V1,R1,R2,R3,R4,R5,R6,R7))
    
#    taurinp=0.1*ms
#    taudinp=0.5*ms    
    taurinp=2*ms
    taudinp=10*ms
#    taurinp=10*ms
#    taudinp=50*ms
    tauinp=taudinp
    
    print('Compiling with cython')
    
    prefs.codegen.target = 'cython' #cython=faster, numpy = default python
    
    net.run(runtime,report='text',report_period=300*second)
    
    if plot_raster :
    
        min_t=int(50*ms*100000*Hz)
        LFP_V_RS=1/N_RS*sum(V1.V,axis=0)[min_t:]
        f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='density')
    
        figure()
        plot(LFP_V_RS)
        
        
        figure()
        plot(f,Spectrum_LFP_V_RS)
        ylabel('Spectrum')
        xlabel('Frequency (Hz)')
        xlim(0,100)
        
        figure()
        plot(f,Spectrum_LFP_V_RS)
        ylabel('Spectrum')
        xlabel('Frequency (Hz)')
        xlim(0,50)
        
        figure()
        plot(R1.t,R1.i+140,'r.',label='RS cells')
        plot(R2.t,R2.i+120,'b.',label='FS cells')
        plot(R3.t,R3.i+100,'g.',label='SOM cells')
        plot([0.2,runtime/second],[95,95],'k--')
        plot(R5.t,R5.i+70,'r.')
        plot(R6.t,R6.i+50,'b.')
        plot([0.2,runtime/second],[45,45],'k--')
        plot(R4.t,R4.i+20,'m.',label='IB cells')
        plot(R7.t,R7.i,'g.')
        xlim(0.2,runtime/second)
        ylim(0,220)
    #    legend(loc='upper left')
        xlabel('Time (s)')
        ylabel('Neuron index')
    

    save_raster('LIP_RS',R1.i,R1.t,new_path)
    save_raster('LIP_FS',R2.i,R2.t,new_path)
    save_raster('LIP_SI',R3.i,R3.t,new_path)
    save_raster('LIP_RS_gran',R5.i,R5.t,new_path)
    save_raster('LIP_FS_gran',R6.i,R6.t,new_path)
    save_raster('LIP_IB',R4.i,R4.t,new_path)
    save_raster('LIP_SI_deep',R7.i,R7.t,new_path)
    
    return



def run_one_simulation(simu,path,index_var):
#    print(simu,len(simu))
    start_scope()
    close('all')

    runtime=1*second
#    runtime=5*second
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp_IB=0* msiemens * cm **-2
    ginp_SI=0* msiemens * cm **-2
    ginp=0* msiemens * cm **-2

    NN=1 #multiplicative factor on the number of neurons
    N_RS,N_FS,N_SI,N_IB= NN*80,NN*20,NN*20,NN*20 #Number of neurons of RE, TC, and HTC type
    
    syn_cond,J,thal,theta_phase,target_time,index=simu
    print('Simulation '+str(index))
    

    net = Network(collect())
    
    print('Network setup')

    all_neurons, all_synapses, all_gap_junctions, all_monitors=make_full_network(syn_cond,J,thal,theta_phase,target_time)
    V1,V2,V3,V4,R1,R2,R3,R4,I1,I2,I3,I4,V5,R5,Is,I5a,I5ad,I5bd,R6,R7,R8,V6,V7,V8,inpmon,inptarget=all_monitors
    
    
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_gap_junctions)
#    net.add(all_monitors)
    net.add((V1,R1,R2,R3,R4,R5,R6,R7,R8,inpmon,inptarget))
    
#    taurinp=0.1*ms
#    taudinp=0.5*ms    
    taurinp=2*ms
    taudinp=10*ms
#    taurinp=10*ms
#    taudinp=50*ms
    tauinp=taudinp
    
    print('Compiling with cython')
    
    prefs.codegen.target = 'cython' #cython=faster, numpy = default python
    
    net.run(runtime,report='text',report_period=300*second)
    
    figure()
    plot(R1.t,R1.i+160,'r.',label='RS cells')
    plot(R2.t,R2.i+140,'m.',label='FS cells')
    plot(R3.t,R3.i+120,'y.',label='SI cells')
    plot(R4.t,R4.i+100,'.',label='VIP cells',color='orange')
    plot(R6.t,R6.i+70,'g.',label='Granular RS')
    plot(R7.t,R7.i+50,'c.',label='Granular FS')
    plot(R5.t,R5.i+20,'b.',label='IB cells')
    plot(R8.t,R8.i,'k.',label='Deep SI')
    xlim(0,runtime/second)
    legend(loc='upper left')
    
    figure()
    plot(inpmon.t,inpmon.sinp[0])
    xlabel('Time (s)')
    ylabel('Input synapse opening variable $s_{input}$')
    tight_layout()
    
    figure()
    plot(inptarget.t,inptarget.i,'k.')
    xlabel('Time (s)')
    ylabel('Target Input')
    
    min_t=int(50*ms*100000*Hz)
    LFP_V_RS=1/N_RS*sum(V1.V,axis=0)[min_t:]

    f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='density')

    
    figure()
    plot(LFP_V_RS)
    
    figure()
    plot(R1.t,R1.i+160,'r.',label='RS cells')
    plot(R2.t,R2.i+140,'b.',label='FS cells')
    plot(R3.t,R3.i+120,'g.',label='SI cells')
    plot(R4.t,R4.i+100,'y.',label='VIP cells')
    plot(R6.t,R6.i+70,'.',label='Granular RS',color='C1')
    plot(R7.t,R7.i+50,'c.',label='Granular FS')
    plot(R5.t,R5.i+20,'m.',label='IB cells')
    plot(R8.t,R8.i,'.',label='Deep SI',color='lime')
    xlim(0,runtime/second)
    ylim(0,240)
    legend(loc='upper left')
    xlabel('Time (s)')
    ylabel('Neuron index')
    
    figure()
    plot(f,Spectrum_LFP_V_RS)
    ylabel('Spectrum')
    xlabel('Frequency (Hz)')
    xlim(0,100)
    
    figure()
    plot(f,Spectrum_LFP_V_RS)
    ylabel('Spectrum')
    xlabel('Frequency (Hz)')
    xlim(0,50)
    
    figure()
    plot(R1.t,R1.i+160,'r.',label='RS cells')
    plot(R2.t,R2.i+140,'b.',label='FS cells')
    plot(R3.t,R3.i+120,'g.',label='SOM cells')
    plot(R4.t,R4.i+100,'k.',label='VIP cells')
    plot([0.2,runtime/second],[95,95],'k--')
    plot(R6.t,R6.i+70,'r.')
    plot(R7.t,R7.i+50,'b.')
    plot([0.2,runtime/second],[45,45],'k--')
    plot(R5.t,R5.i+20,'m.',label='IB cells')
    plot(R8.t,R8.i,'g.')
    xlim(0.2,runtime/second)
    ylim(0,240)
    xticks(fontsize=12)
    yticks(fontsize=12)
#    legend(loc='upper left')
    xlabel('Time (s)',fontsize=12)
    ylabel('Neuron index',fontsize=12)
    
    
    ##save figures
    new_path=path+"/results_"+str(index)
    os.mkdir(new_path)

    for n in get_fignums():
        current_fig=figure(n)
        current_fig.savefig(new_path+'/figure'+str(n)+'.png')
        
    save_raster('LIP_RS',R1.i,R1.t,new_path)
    save_raster('LIP_FS',R2.i,R2.t,new_path)
    save_raster('LIP_SI',R3.i,R3.t,new_path)
    save_raster('LIP_VIP',R4.i,R4.t,new_path)
    save_raster('LIP_RS_gran',R6.i,R6.t,new_path)
    save_raster('LIP_FS_gran',R7.i,R7.t,new_path)
    save_raster('LIP_IB',R5.i,R5.t,new_path)
    save_raster('LIP_SI_deep',R8.i,R8.t,new_path)

    
if __name__=='__main__':
    
    runtime=1*second
    all_theta=['good','bad'] 
    all_t_SOM=[20*msecond]
    all_t_FS=[5*msecond]

    simu=(target_presentation_time,0,t_SOM,t_FS,theta_phase,g_LIP_FEF_v,target_presence,runtime)
   
    path="./results_"+str(datetime.datetime.now())
    os.mkdir(path)
        
    all_sim=list(product(all_t_SOM,all_t_FS,all_theta))
#    index_var=[-1]
    
    all_sim=[[800*second,0]+list(all_sim[i])+[0,True,runtime,i] for i in range(len(all_sim))]
    
    #saving simulation parameters as a txt file
    param_file=open(path+'/parameters.txt','w')
    for simu in all_sim:
        param_file.write(str(simu))
        param_file.write('\n\n')
    param_file.close()
    
    for simu in all_sim:
        new_path=path+'/'+str(simu[-1])
        os.mkdir(new_path)
        run_one_LIP_simulation(simu,new_path,plot_raster=True)
    
        print('Saving figures')
        os.mkdir(new_path+'/figures')
        for i in get_fignums():
            current_fig=figure(i)
            current_fig.savefig(new_path+'/figures/figure'+str(i)+'.png')

    
    clear_cache('cython')    