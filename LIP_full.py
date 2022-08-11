# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 14:28:44 2019

@author: aaussel
"""
#import matplotlib
#matplotlib.use('Agg')

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

#from joblib import Parallel, delayed
#from joblib import parallel_backend
#import multiprocessing
import os

os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_DYNAMIC'] = 'FALSE'

import time
from itertools import *

    
##Custom
#input_beta2_IB=False
#input_beta2_RS=False
#input_beta2_FS_SI=False
#input_thalamus_gran=False
#thal_cond=8* msiemens * cm **-2
#kainate='low'

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

def zeros_ones_monitor(spikemon,record_dt,runtime):
    L=int(runtime/record_dt)
    zeros_ones=[0]*L
    for time in spikemon.t:
        zeros_ones[int(time/record_dt)]+=1
    return zeros_ones
    
def make_full_network(syn_cond,J,thal,theta_phase,target_time):
    
#    print(syn_cond,J,thal,theta_phase)    
    
    NN=1 #multiplicative factor on the number of neurons
    N_RS,N_FS,N_SI,N_IB= NN*80,NN*20,NN*20,NN*20 #Number of neurons of RE, TC, and HTC type
    
    gSIdFSg,gFSgRSg,gRSgFSg,gRSgRSg,gFSgFSg,gRSgRSs,gRSgFSs,gFSgRSs=syn_cond
    J_RSg,J_FSg=J
    
 
    FLee=(0.05*mS/cm**2)/(0.4*uS/cm**2)*0.5
    version = 'Alex'
    runtime=1*second
    kainate='low'
    
    all_neurons, all_synapses, all_gap_junctions, all_monitors=create_Mark_Alex_network(kainate,version,Nf=NN)
    V1,V2,V3,V4,R1,R2,R3,R4,I1,I2,I3,I4,V5,R5,I5s,I5a,I5ad,I5IBbd=all_monitors
    RS, FS, SI, VIP, IB_soma, IB_axon, IB_bd, IB_ad =all_neurons
   
    if theta_phase=='bad':
        input_beta2_IB=False
        input_beta2_RS=False
        input_beta2_FS_SI=True
        input_thalamus_gran=True
        gFS=0* msiemens * cm **-2
        SI.ginp_SI=0* msiemens * cm **-2
        thal_cond=2* msiemens * cm **-2
        input_mixed=False
        target_on=True
        
        
    if theta_phase=='good':
        input_beta2_IB=True
        IB_bd.ginp_IB=500* msiemens * cm **-2
#        IB_bd.ginp_IB=0* msiemens * cm **-2
        input_beta2_RS=False
        input_beta2_FS_SI=False
        input_thalamus_gran=True
        thal_cond=thal
#        thal_cond=thal*2754.660086037123/12782.0904181147
#        thal_cond=thal*2754.660086037123/139.46773954954165
        input_mixed=False
        target_on=True
        
    if theta_phase=='mixed':
        input_mixed=True
        IB_bd.ginp_IB=500* msiemens * cm **-2
        input_beta2_IB=False
        input_beta2_RS=False
        input_beta2_RS=False
        input_beta2_FS_SI=False
        input_thalamus_gran=False
        thal_cond=thal
        target_on=True
        
#    print(input_mixed,input_beta2_IB,input_beta2_RS,input_beta2_FS_SI,input_thalamus_gran)
        
    prefs.codegen.target = 'numpy'
    defaultclock.dt = 0.01*ms
    
    #Single column network
    
    ##Define neuron groups
    E_gran=NeuronGroup(N_FS,eq_RS_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    E_gran.V = '-70*mvolt+10*rand()*mvolt'
    E_gran.h = '0+0.05*rand()'
    E_gran.m = '0+0.05*rand()'
    E_gran.mAR = '0.035+0.025*rand()'
    if kainate=='low':
    #    E_gran.J='1 * uA * cmeter ** -2'  #article SI=25, code=1
        E_gran.J=J_RSg
    elif kainate=='high':
        E_gran.J='-10 * uA * cmeter ** -2'  #article SI=25, code=1
#        E_gran.J='-5 * uA * cmeter ** -2'  #article SI=25, code=1
        
    FS_gran=NeuronGroup(N_FS,eq_FS_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    FS_gran.V = '-110*mvolt+10*rand()*mvolt'
    FS_gran.h = '0+0.05*rand()'
    FS_gran.m = '0+0.05*rand()'
    if kainate=='low':
    #    FS_gran.J='5 * uA * cmeter ** -2' #article=code=35
        FS_gran.J=J_FSg
    elif kainate=='high':
        FS_gran.J='16 * uA * cmeter ** -2'
    
    SI_deep=NeuronGroup(N_SI,eq_SI_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    SI_deep.V = '-100*mvolt+10*rand()*mvolt'
    SI_deep.h = '0+0.05*rand()'
    SI_deep.m = '0+0.05*rand()'
    SI_deep.mAR = '0.02+0.04*rand()'
    if version == 'Alex':
        if kainate=='low':
            SI_deep.J='35* uA * cmeter ** -2' #article SI=50, code=35, Mark = 45
        elif kainate=='high':
            SI_deep.J='30* uA * cmeter ** -2' #article SI=50, code=35, Mark = 45
    elif version=='Mark':
        if kainate=='low':
            SI_deep.J='45* uA * cmeter ** -2' #article SI=50, code=35, Mark = 45
        elif kainate=='high':
            SI_deep.J='40* uA * cmeter ** -2' #article SI=50, code=35, Mark = 45
    
    if theta_phase=='good' or theta_phase=='mixed':
#        SI_deep.ginp_SI=50* msiemens * cm **-2
        SI_deep.ginp_SI=5* msiemens * cm **-2
#        SI_deep.ginp_SI=0* msiemens * cm **-2
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
            S.connect(j=connection_pattern, skip_if_invalid=True)
        S.g_i=g_i
        S.taur_i=taur_i
        S.taud_i=taud_i
        S.V_i=V_i  
        return S
    
    #From E (granular layer) cells
    #S_EgranFS=generate_syn(E_gran,FS,'IsynEgran','',0.2*usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)
    #S_EgranFS=generate_syn(E_gran,FS,'IsynEgran','',0.1*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_EgranFS=generate_syn(E_gran,FS,'IsynRS_LIP_gran','i//10==j//10',gRSgFSs,0.125*ms,1*ms,0*mV)
    #S_EgranEgran=generate_syn(E_gran,E_gran,'IsynEgran','',0.4*usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)
    S_EgranEgran=generate_syn(E_gran,E_gran,'IsynRS_LIP_gran','i//10==j//10',gRSgRSg,0.125*ms,1*ms,0*mV)
    #S_EgranFSgran=generate_syn(E_gran,FS_gran,'IsynEgran','',0.2*usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)
    S_EgranFSgran=generate_syn(E_gran,FS_gran,'IsynRS_LIP_gran','i//10==j//10',gRSgFSg,0.125*ms,1*ms,0*mV)
    #S_EgranRS=generate_syn(E_gran,RS,'IsynEgran','',0.2*usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)
    #S_EgranRS=generate_syn(E_gran,RS,'IsynEgran','',1*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_EgranRS=generate_syn(E_gran,RS,'IsynRS_LIP_gran','i//10==j//10',gRSgRSs,0.125*ms,1*ms,0*mV)
    
    S_EgranIB=generate_syn(E_gran,IB_ad,'IsynRS_LIP_gran','i//10==j//10',0.212*usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)
    
    #From FS (granular layer) cells, normal timescale
    #S_FSgranEgran=generate_syn(FS_gran,E_gran,'IsynFSgran','',1* usiemens * cm **-2*FLee,0.25*ms,5*ms,-80*mV)
    S_FSgranEgran=generate_syn(FS_gran,E_gran,'IsynFS_LIP_gran','i//10==j//10',gFSgRSg,0.25*ms,5*ms,-80*mV)
    S_FSgranFSgran=generate_syn(FS_gran,FS_gran,'IsynFS_LIP_gran','i//10==j//10',gFSgFSg,0.25*ms,5*ms,-75*mV)
    #S_FSgranRS=generate_syn(FS_gran,RS,'IsynFSgran','',0.02* usiemens * cm **-2*FLee,0.25*ms,5*ms,-80*mV)
    S_FSgranRS=generate_syn(FS_gran,RS,'IsynFS_LIP_gran','i//10==j//10',gFSgRSs,0.25*ms,5*ms,-80*mV)
    #S_FSgranRS=generate_syn(FS_gran,RS,'IsynFSgran','',0* msiemens * cm **-2,0.25*ms,5*ms,-80*mV)

#    #From FS (granular layer) cells, beta timescale
#    S_FSgranEgran=generate_syn(FS_gran,E_gran,'IsynFS_LIP_gran','',gFSgRSg,0.25*ms,20*ms,-80*mV)
#    S_FSgranFSgran=generate_syn(FS_gran,FS_gran,'IsynFS_LIP_gran','',gFSgFSg,0.25*ms,20*ms,-75*mV)
#    S_FSgranRS=generate_syn(FS_gran,RS,'IsynFS_LIP_gran','',gFSgRSs,0.25*ms,20*ms,-80*mV)
    
    #From IB cells
    #S_IBSIdeep=generate_syn(IB_axon,SI_deep,'IsynIB','',0.12* usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)
#    S_IBSIdeep=generate_syn(IB_axon,SI_deep,'IsynIB','',0.2* msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_IBSIdeep=generate_syn(IB_axon,SI_deep,'IsynIB_LIP','i//10==j//10',0.01* msiemens * cm **-2,0.125*ms,1*ms,0*mV)
        
    
    #From deep SI cells    
    S_SIdeepIB=generate_syn(SI_deep,IB_bd,'IsynSI_LIP_deep','i//10==j//10',10* msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
    #S_SIdeepFSgran=generate_syn(SI_deep,FS_gran,'IsynSIdeep','',0.4* usiemens * cm **-2*FLee,0.25*ms,20*ms,-80*mV)
    S_SIdeepFSgran=generate_syn(SI_deep,FS_gran,'IsynSI_LIP_deep','i//10==j//10',gSIdFSg,0.25*ms,20*ms,-80*mV)
    
    
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
    
    if input_beta2_IB:
#        inputs_topdown=generate_spike_timing(N_IB,25*Hz,0*ms,end_time=3000*ms)
#        G_topdown = SpikeGeneratorGroup(N_IB, inputs_topdown[:,1], inputs_topdown[:,0]*second)
#        topdown_in=Synapses(G_topdown,IB_bd,on_pre='Vinp=Vhigh')
#        topdown_in.connect(j='i')
        
#        inputs_topdown3=generate_spike_timing(N_SI,50*Hz,0*ms,end_time=8000*ms)
        inputs_topdown3=generate_spike_timing(N_SI,25*Hz,0*ms,end_time=10000*ms)
        G_topdown3 = SpikeGeneratorGroup(N_SI, inputs_topdown3[:,1], inputs_topdown3[:,0]*second)
        topdown_in3=Synapses(G_topdown3,SI_deep,on_pre='Vinp=Vhigh')
        topdown_in3.connect(j='i')
        
    if input_beta2_RS:    
        RS.ginp_RS_good=4* msiemens * cm **-2
        RS.ginp_RS_bad=4* msiemens * cm **-2
        inputs_topdown2=generate_spike_timing(N_RS,25*Hz,0*ms,end_time=10000*ms)
        G_topdown2 = SpikeGeneratorGroup(N_RS, inputs_topdown2[:,1], inputs_topdown2[:,0]*second)
        topdown_in2=Synapses(G_topdown2,RS,on_pre='Vinp=Vhigh')
        topdown_in2.connect(j='i')
        
    if input_beta2_FS_SI:
        FS.ginp_FS_good=gFS
        FS.ginp_FS_bad=gFS
        inputs_lateral=generate_spike_timing(N_FS,25*Hz,0*ms,end_time=10000*ms)
        G_lateral = SpikeGeneratorGroup(N_FS, inputs_lateral[:,1], inputs_lateral[:,0]*second)
        lateral_in=Synapses(G_lateral,FS,on_pre='Vinp=Vhigh')
        lateral_in.connect(j='i')
        
        inputs_lateral2=generate_spike_timing(N_SI,25*Hz,0*ms,end_time=10000*ms)
        G_lateral2 = SpikeGeneratorGroup(N_SI, inputs_lateral2[:,1], inputs_lateral2[:,0]*second)
        lateral_in2=Synapses(G_lateral2,SI,on_pre='Vinp=Vhigh')
        lateral_in2.connect(j='i')
        
    if input_thalamus_gran:
        E_gran.ginp_RS_good=thal_cond
        FS_gran.ginp_FS_good=thal_cond
        E_gran.ginp_RS_bad=thal_cond
        FS_gran.ginp_FS_bad=thal_cond
#        print(E_gran.ginp_RS_good,FS_gran.ginp_FS_good,E_gran.ginp_RS_bad,FS_gran.ginp_FS_bad)
#        FS_gran.ginp_FS=thal_cond
        inputs_mdpul=generate_spike_timing(N_FS,13*Hz,0*ms,end_time=10000*ms)
        Poisson_input = SpikeGeneratorGroup(N_FS, inputs_mdpul[:,1], inputs_mdpul[:,0]*second)
#        Poisson_input = PoissonGroup(N_FS,100*Hz)
        bottomup_in = Synapses(Poisson_input,FS_gran, on_pre='Vinp=Vhigh')
        bottomup_in.connect(j='i')
#        E_gran.ginp_RS=thal_cond
        Poisson_input2 = SpikeGeneratorGroup(N_FS, inputs_mdpul[:,1], inputs_mdpul[:,0]*second)
#        Poisson_input2 = PoissonGroup(N_FS,100*Hz)
        bottomup_in2 = Synapses(Poisson_input2,E_gran, on_pre='Vinp=Vhigh')
        bottomup_in2.connect(j='i')
    #    print(bottomup_in,bottomup_in2)
        
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
        SI.ginp_SI2=2.5* msiemens * cm **-2
        VIP.ginp_VIP2=2.5* msiemens * cm **-2
        RS.ginp_RS2=2.5* msiemens * cm **-2
    
    if input_mixed:
        E_gran.ginp_RS_good=thal_cond
        FS_gran.ginp_FS_good=thal_cond
#        E_gran.ginp_RS_bad=2* msiemens * cm **-2
#        FS_gran.ginp_FS_bad=2* msiemens * cm **-2
        E_gran.ginp_RS_bad=thal_cond
        FS_gran.ginp_FS_bad=thal_cond
#        FS_gran.ginp_FS='15* msiemens * cm **-2 * int(sin(2*pi*t*4*Hz)>0) + 2* msiemens * cm **-2 * int(sin(2*pi*t*4*Hz)<0)'
#        inputs_mdpul=generate_spike_timing_theta(N_FS,13*Hz,0*ms,end_time=3000*ms)

        #4Hz theta
        t0=0*ms
        t1=125*ms
        inputs_mdpul=generate_spike_timing(N_SI,13*Hz,t0,end_time=t1)
        while t0+250*ms<runtime:
            t0,t1=t0+250*ms,t1+250*ms
            inputs_mdpul=vstack((inputs_mdpul,generate_spike_timing(N_SI,13*Hz,t0,end_time=t1)))

#        #8Hz theta
#        t0=0*ms
#        t1=62.5*ms
#        inputs_mdpul=generate_spike_timing(N_SI,13*Hz,t0,end_time=t1)
#        while t0+125*ms<runtime:
#            t0,t1=t0+125*ms,t1+125*ms
#            inputs_mdpul=vstack((inputs_mdpul,generate_spike_timing(N_SI,13*Hz,t0,end_time=t1)))

#        #2Hz theta
#        t0=0*ms
#        t1=250*ms
#        inputs_mdpul=generate_spike_timing(N_SI,13*Hz,t0,end_time=t1)
#        while t0+500*ms<runtime:
#            t0,t1=t0+500*ms,t1+500*ms
#            inputs_mdpul=vstack((inputs_mdpul,generate_spike_timing(N_SI,13*Hz,t0,end_time=t1)))

#        #16Hz theta
#        t0=0*ms
#        t1=31.25*ms
#        inputs_mdpul=generate_spike_timing(N_SI,13*Hz,t0,end_time=t1)
#        while t0+62.5*ms<runtime:
#            t0,t1=t0+62.5*ms,t1+62.5*ms
#            inputs_mdpul=vstack((inputs_mdpul,generate_spike_timing(N_SI,13*Hz,t0,end_time=t1)))
                                                            
        
        Poisson_input = SpikeGeneratorGroup(N_FS, inputs_mdpul[:,1], inputs_mdpul[:,0]*second)
#        Poisson_input = PoissonGroup(N_FS,100*Hz)
        bottomup_in = Synapses(Poisson_input,FS_gran, on_pre='Vinp=Vhigh')
        bottomup_in.connect(j='i')
#        E_gran.ginp_RS='15* msiemens * cm **-2 * int(sin(2*pi*t*4*Hz)>0) + 2* msiemens * cm **-2 * int(sin(2*pi*t*4*Hz)<0)'
        Poisson_input2 = SpikeGeneratorGroup(N_FS, inputs_mdpul[:,1], inputs_mdpul[:,0]*second)
#        Poisson_input2 = PoissonGroup(N_FS,100*Hz)
        bottomup_in2 = Synapses(Poisson_input2,E_gran, on_pre='Vinp=Vhigh')
        bottomup_in2.connect(j='i')
        
        inputs_topdown=generate_spike_timing_theta(N_IB,25*Hz,0*ms,end_time=5100*ms)
        G_topdown = SpikeGeneratorGroup(N_IB, inputs_topdown[:,1], inputs_topdown[:,0]*second)
        topdown_in=Synapses(G_topdown,IB_bd,on_pre='Vinp=Vhigh')
        topdown_in.connect(j='i')
        
        #theta=4Hz
        inputs_topdown3=generate_spike_timing_theta(N_SI,25*Hz,0*ms,end_time=5100*ms)
#        #theta=8Hz
#        inputs_topdown3=generate_spike_timing_theta(N_SI,25*Hz,0*ms,end_time=5100*ms,f_theta=8*Hz)
#        #theta=2Hz
#        inputs_topdown3=generate_spike_timing_theta(N_SI,25*Hz,0*ms,end_time=5100*ms,f_theta=2*Hz)
#        #theta=16Hz
#        inputs_topdown3=generate_spike_timing_theta(N_SI,25*Hz,0*ms,end_time=5100*ms,f_theta=16*Hz)
#        inputs_topdown3=generate_spike_timing_theta(N_SI,50*Hz,0*ms,end_time=5100*ms)
        G_topdown3 = SpikeGeneratorGroup(N_SI, inputs_topdown3[:,1], inputs_topdown3[:,0]*second)
        topdown_in3=Synapses(G_topdown3,SI_deep,on_pre='Vinp=Vhigh')
        topdown_in3.connect(j='i')
    
#    g_inputs=[G_topdown,G_topdown2,G_topdown3,G_lateral,G_lateral2,Poisson_input,Poisson_input2]
#    g_inputs=[y for y in g_inputs if y]
#    syn_inputs=[topdown_in,topdown_in2,topdown_in3,lateral_in,lateral_in2,bottomup_in,bottomup_in2]
#    syn_inputs=[y for y in syn_inputs if y]
    numpy.set_printoptions(threshold=sys.maxsize)
    
#    print(E_gran.ginp_RS_good)
#    print(Poisson_input2)
    
    g_inputs=[G_topdown2,G_topdown3,G_lateral,G_lateral2,Poisson_input,Poisson_input2]
    g_inputs=[y for y in g_inputs if y]
    syn_inputs=[topdown_in2,topdown_in3,lateral_in,lateral_in2,bottomup_in,bottomup_in2]
    syn_inputs=[y for y in syn_inputs if y]
    
#    print(len(g_inputs))
#    print(len(syn_inputs))
#    print(Poisson_input2 in g_inputs)
    
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
    inpIBmon=StateMonitor(SI_deep,'Iapp',record=[0])

    all_neurons=all_neurons+(E_gran,FS_gran,SI_deep)+tuple(g_inputs)
    all_synapses=all_synapses+(S_EgranFS,S_EgranEgran,S_EgranFSgran,S_EgranRS,S_EgranIB,S_FSgranEgran,S_FSgranFSgran,S_FSgranRS,S_IBSIdeep,S_SIdeepIB,S_SIdeepFSgran)+tuple(syn_inputs)
    all_monitors=all_monitors+(R6,R7,R8,V6,V7,V8,inpmon,inpIBmon)
    return all_neurons, all_synapses, all_gap_junctions, all_monitors


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
    
    syn_cond,J,thal,theta_phase,index,target_time=simu
    print('Simulation '+str(index))
    
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
        input_mixed=False
        
    if theta_phase=='good':
#        input_beta2_IB=True
        input_beta2_IB=False
        ginp_IB=500* msiemens * cm **-2
#        ginpSIdeep=500* msiemens * cm **-2
        ginpSIdeep=0* msiemens * cm **-2
        input_beta2_RS=False
        input_beta2_FS_SI=False
        input_thalamus_gran=True
        thal_cond=thal
        kainate='low'
        input_mixed=False
        
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
    
    net = Network(collect())
    
    print('Network setup')
    all_neurons, all_synapses, all_gap_junctions, all_monitors=make_full_network(syn_cond,J,thal,theta_phase)
    V1,V2,V3,V4,R1,R2,R3,R4,I1,I2,I3,I4,V5,R5,Is,I5a,I5ad,I5bd,R6,R7,R8,V6,V7,V8,inpmon,inpIBmon=all_monitors
    
    
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_gap_junctions)
#    net.add(all_monitors)
    net.add((V1,R1,R2,R3,R4,R5,R6,R7,R8,inpmon,inpIBmon))
    
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
    legend(loc='upper left',fontsize=12)
    xticks(fontsize=12)
    yticks(fontsize=12)
    
    figure()
#    plot(inpmon.t,inpmon.Iinp1[0])
    plot(inpmon.t,inpmon.sinp[0])
    xlabel('Time (s)')
    ylabel('Input synapse opening variable $s_{input}$')
    tight_layout()
    
    print(sum(inpmon.sinp[0]))
#    print(inpmon.t)
#    print(inpmon.Iinp1[0])
#    
#    print(inpIBmon.t)
#    print(inpIBmon.Iapp[0])
    
    figure()
    plot(inpIBmon.t,inpIBmon.Iapp[0])
    
    min_t=int(50*ms*100000*Hz)
#    min_t=int(150*ms*100000*Hz)
    LFP_V_RS=1/N_RS*sum(V1.V,axis=0)[min_t:]
#    LFP_V_FS=1/N_FS*sum(V2.V,axis=0)[min_t:]
#    LFP_V_SI=1/N_SI62.903225806451616*sum(V3.V,axis=0)[min_t:]
#    LFP_V_IB=1/N_IB*sum(V5.V,axis=0)[min_t:]
#    LFP_V_RSg=1/N_FS*sum(V5.V,axis=0)[min_t:]
#    LFP_V_FSg=1/N_FS*sum(V7.V,axis=0)[min_t:]
#    LFP_V_SId=1/N_SI*sum(V8.V,axis=0)[min_t:]
    
#    f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='density')
#    f,Spectrum_LFP_V_FS=signal.periodogram(LFP_V_FS, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_SI=signal.periodogram(LFP_V_SI, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_IB=signal.periodogram(LFP_V_IB, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_RSg=signal.periodogram(LFP_V_RSg, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_FSg=signal.periodogram(LFP_V_FSg, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_SId=signal.periodogram(LFP_V_SId, 100000,'flattop', scaling='spectrum')
    
    figure()
    plot(LFP_V_RS)
    xticks(fontsize=12)
    yticks(fontsize=12)
    
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
    
#    tight_layout()
    
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
    ylim(0,220)
    legend(loc='upper left')
    xlabel('Time (s)')
    ylabel('Neuron index')
    
    figure()
    plot(f,Spectrum_LFP_V_RS)
    ylabel('Spectrum',fontsize=12)
    xlabel('Frequency (Hz)',fontsize=12)
    xticks(fontsize=12)
    yticks(fontsize=12)
    xlim(0,100)
    
    figure()
    plot(f,Spectrum_LFP_V_RS)
    ylabel('Spectrum',fontsize=12)
    xlabel('Frequency (Hz)',fontsize=12)
    xticks(fontsize=12)
    yticks(fontsize=12)
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
    ylim(0,220)
    xticks(fontsize=12)
    yticks(fontsize=12)
#    legend(loc='upper left')
    xlabel('Time (s)',fontsize=12)
    ylabel('Neuron index',fontsize=12)
    
#    N_RS_spikes=zeros_ones_monitor(R1,defaultclock.dt,runtime)
    N_RS_spikes=zeros_ones_monitor(R4,defaultclock.dt,runtime)
    figure()
    plot(N_RS_spikes)
    f, t, Sxx = signal.spectrogram(array(N_RS_spikes), 100000*Hz,nperseg=20000,noverlap=15000)
    figure()
    pcolormesh(t, f, Sxx)#, cmap=)
    colorbar(format='%.1e')
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    ylim(0,100)
    title('Power ($V^2$)')
    
    f,Spectrum_spike_RS=signal.periodogram(array(N_RS_spikes), 100000,'flattop', scaling='density')
    figure()
    plot(f,Spectrum_spike_RS)
    xlabel('Frequency (Hz)')   
    ylabel('Spectrum')
    xlim(0,100)
    
    L=int(runtime/defaultclock.dt)
    zeros_ones=[0]*L
    ind=0
    for time in R1.t:
        if R1.i[ind]==1:
            zeros_ones[int(time/defaultclock.dt)]+=1
        ind+=1
    f,Spectrum_spike_RS=signal.periodogram(array(zeros_ones), 100000,'flattop', scaling='density')
    figure()
    plot(f,Spectrum_spike_RS)
    xlabel('Frequency (Hz)')   
    ylabel('Spectrum')
    xlim(0,100)    
    

    
#    min_t=int(50*ms*100000*Hz)
#    LFP_V_RS=1/N_RS*sum(V1.V,axis=0)[min_t:]
#    LFP_V_FS=1/N_FS*sum(V2.V,axis=0)[min_t:]
#    LFP_V_SI=1/N_SI*sum(V3.V,axis=0)[min_t:]
#    LFP_V_IB=1/N_IB*sum(V5.V,axis=0)[min_t:]
#    
#    f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_FS=signal.periodogram(LFP_V_FS, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_SI=signal.periodogram(LFP_V_SI, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_IB=signal.periodogram(LFP_V_IB, 100000,'flattop', scaling='spectrum')
#    
#    figure()
#    subplot(421)
#    plot((V1.t/second)[min_t:],LFP_V_RS)
#    ylabel('LFP')
#    title('RS cell')
#    subplot(423)
#    plot((V1.t/second)[min_t:],LFP_V_FS)
#    ylabel('LFP')
#    title('FS cell')
#    subplot(425)
#    plot((V1.t/second)[min_t:],LFP_V_SI)
#    ylabel('LFP')
#    title('SI cell')
#    subplot(427)
#    plot((V1.t/second)[min_t:],LFP_V_IB)
#    xlabel('Time (s)')
#    ylabel('LFP')
#    title('IB cell')
#    
#    subplot(422)
#    plot(f,Spectrum_LFP_V_RS)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('RS cell')
#    subplot(424)
#    plot(f,Spectrum_LFP_V_FS)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('FS cell')
#    subplot(426)
#    plot(f,Spectrum_LFP_V_SI)
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('SI cell')
#    subplot(428)
#    plot(f,Spectrum_LFP_V_IB)
#    xlabel('Frequency (Hz)')
#    ylabel('Spectrum')
#    yticks([],[])
#    xlim(0,100)
#    title('IB cell')
    
    ##save figures
    new_path=path+"/results_"+str(index)
    os.mkdir(new_path) 

    for n in get_fignums():
        current_fig=figure(n)
        current_fig.savefig(new_path+'/figure'+str(n)+'.png')
        current_fig.savefig(new_path+'/figure'+str(n)+'.eps')
        
    save_raster('LIP_RS',R1.i,R1.t,new_path)
    save_raster('LIP_FS',R2.i,R2.t,new_path)
    save_raster('LIP_SI',R3.i,R3.t,new_path)
    save_raster('LIP_VIP',R4.i,R4.t,new_path)
    save_raster('LIP_RS_gran',R6.i,R6.t,new_path)
    save_raster('LIP_FS_gran',R7.i,R7.t,new_path)
    save_raster('LIP_IB',R5.i,R5.t,new_path)
    save_raster('LIP_SI_deep',R8.i,R8.t,new_path)

    
if __name__=='__main__':
        
    Vlow=-80*mV
    FLee=(0.05*mS/cm**2)/(0.4*uS/cm**2)*0.5   
    all_SIdFSg=[2*msiemens * cm **-2]
    all_FSgRSg=[1* msiemens * cm **-2]
    all_RSgFSg=[1*msiemens * cm **-2]
    all_RSgRSg=[0.5*msiemens * cm **-2]
    all_FSgFSg=[0.3* msiemens * cm **-2]
    all_RSgRSs=[2*msiemens * cm **-2]
    all_RSgFSs=[0.1*msiemens * cm **-2]
    all_FSgRSs=[0.1* msiemens * cm **-2]
    all_J_RSg=['10 * uA * cmeter ** -2'] 
    all_J_FSg=['-5 * uA * cmeter ** -2']

##   A
#    all_J_RSg=['-5 * uA * cmeter ** -2'] #['20 * uA * cmeter ** -2']
#    all_J_FSg=['-10 * uA * cmeter ** -2']

##   B
#    all_J_RSg=['-5 * uA * cmeter ** -2'] #['20 * uA * cmeter ** -2']
#    all_J_FSg=['5 * uA * cmeter ** -2']

##   C
#    all_J_RSg=['20 * uA * cmeter ** -2'] #['20 * uA * cmeter ** -2']
#    all_J_FSg=['-5 * uA * cmeter ** -2']

#    all_J_RSg=['-10 * uA * cmeter ** -2','-5 * uA * cmeter ** -2','0 * uA * cmeter ** -2','5 * uA * cmeter ** -2','10 * uA * cmeter ** -2','15 * uA * cmeter ** -2','20 * uA * cmeter ** -2']
#    all_J_FSg=['-9 * uA * cmeter ** -2','-8 * uA * cmeter ** -2','-7 * uA * cmeter ** -2','-6 * uA * cmeter ** -2','1 * uA * cmeter ** -2','2 * uA * cmeter ** -2','3 * uA * cmeter ** -2','4 * uA * cmeter ** -2']
#    all_thal=[10* msiemens * cm **-2]
    all_thal=[10* msiemens * cm **-2]
#    all_thal=[0* msiemens * cm **-2]
#    all_theta=['mixed']
    all_theta=['good']
#    all_theta=['mixed','mixed','mixed','mixed','mixed']
    
    all_target_time=800*ms
    
    #FLee=(0.05*mS/cm**2)/(0.4*uS/cm**2)*0.5   
    #all_SIdFSg=[1*msiemens * cm **-2]
    #all_FSgRSg=[1* msiemens * cm **-2]
    #all_RSgFSg=[0.5*msiemens * cm **-2]
    #all_RSgRSg=[0.3*msiemens * cm **-2]
    #all_FSgFSg=[0.3* msiemens * cm **-2]
    #all_RSgRSs=[0.25*msiemens * cm **-2]
    #all_RSgFSs=[0.01*msiemens * cm **-2]
    #all_FSgRSs=[0.01* msiemens * cm **-2]
    #all_J_RSg=['30 * uA * cmeter ** -2']
    #all_J_FSg=['30 * uA * cmeter ** -2']
    #all_thal=[15* msiemens * cm **-2]
    #all_theta=['good','bad']
    
    
    all_syn_cond=list(product(all_SIdFSg,all_FSgRSg,all_RSgFSg,all_RSgRSg,all_FSgFSg,all_RSgRSs,all_RSgFSs,all_FSgRSs))
    all_J=list(product(all_J_RSg,all_J_FSg))
    
    path="./results_"+str(datetime.datetime.now())
    os.mkdir(path)
        
    all_sim=list(product(all_syn_cond,all_J,all_thal,all_theta,all_target_time))
    index_var=[-1]
    
    all_sim=[list(all_sim[i])+[i] for i in range(len(all_sim))]
    
    #saving simulation parameters as a txt file
    param_file=open(path+'/parameters.txt','w')
    for simu in all_sim:
        param_file.write(str(simu))
        param_file.write('\n\n')
    param_file.close()
    
    #close('all')
    #all_results=[]
    #num_cores = multiprocessing.cpu_count()
    #num_cores = 2
    
    #with parallel_backend('multiprocessing'):
    #    Parallel(n_jobs=num_cores)(delayed(run_one_simulation)(simu,path,index_var) for simu in all_sim)
    
    for simu in all_sim:
        run_one_simulation(simu,path,index_var)
    
    clear_cache('cython')    