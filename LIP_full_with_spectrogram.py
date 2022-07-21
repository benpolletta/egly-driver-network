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
#from cells.SI_LIP_deep import *
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
    
def make_full_network(syn_cond,J,thal,theta_phase):
    
#    print(syn_cond,J,thal,theta_phase)    
    
    NN=1 #multiplicative factor on the number of neurons
    N_RS,N_FS,N_SI,N_IB= NN*80,NN*20,NN*20,NN*20 #Number of neurons of RE, TC, and HTC type
    
    gSIdFSg,gFSgRSg,gRSgFSg,gRSgRSg,gFSgFSg,gRSgRSs,gRSgFSs,gFSgRSs=syn_cond
    J_RSg,J_FSg=J
    
    FLee=(0.05*mS/cm**2)/(0.4*uS/cm**2)*0.5
    version = 'Alex'
    runtime=3*second
    kainate='low'
    
    timeslots=zeros((int(around(runtime*10/second)),1))
    sinp_SI=TimedArray(timeslots, dt=100*ms)
    
    all_neurons, all_synapses, all_gap_junctions, all_monitors=create_Mark_Alex_network(kainate,version,Nf=NN)
    V1,V2,V3,R1,R2,R3,I1,I2,I3,V4,R4,I4s,I4a,I4ad,I4bd=all_monitors
    RS, FS, SI, IB_soma,IB_axon,IB_bd,IB_ad =all_neurons
    
    SI.ginp_SI=0* msiemens * cm **-2
   
    if theta_phase=='bad':
        input_beta2_IB=False
        input_beta2_RS=False
        input_beta2_FS_SI=True
        input_thalamus_gran=True
        gFS=0* msiemens * cm **-2
#        SI.ginp_SI1=0* msiemens * cm **-2
        SI.ginp_SI=0* msiemens * cm **-2
        thal_cond=2* msiemens * cm **-2
        input_mixed=False
        
    if theta_phase=='good':
        input_beta2_IB=True
        IB_bd.ginp_IB=500* msiemens * cm **-2
        input_beta2_RS=False
        input_beta2_FS_SI=False
        input_thalamus_gran=True
        thal_cond=thal
#        thal_cond=thal#*2754.660086037123/12782.0904181147
#        thal_cond=thal*2754.660086037123/139.46773954954165
        input_mixed=False
        
    if theta_phase=='mixed':
        input_mixed=True
        IB_bd.ginp_IB=500* msiemens * cm **-2
        input_beta2_IB=False
        input_beta2_RS=False
        input_beta2_RS=False
        input_beta2_FS_SI=False
        input_thalamus_gran=False
        thal_cond=thal
        
#    print(input_mixed,input_beta2_IB,input_beta2_RS,input_beta2_FS_SI,input_thalamus_gran)
        
    prefs.codegen.target = 'numpy'
    defaultclock.dt = 0.01*ms
    
    #Single column network
    
    ##Define neuron groups
    E_gran=NeuronGroup(N_FS,eq_RS_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4',name='RSgranLIP')
    E_gran.V = '-70*mvolt+10*rand()*mvolt'
    E_gran.h = '0+0.05*rand()'
    E_gran.m = '0+0.05*rand()'
    E_gran.mAR = '0.035+0.025*rand()'
    if kainate=='low':
    #    E_gran.J='1 * uA * cmeter ** -2'  #article SI=25, code=1
        E_gran.J=J_RSg
    elif kainate=='high':
        E_gran.J='-10 * uA * cmeter ** -2'  #article SI=25, code=1
        
    FS_gran=NeuronGroup(N_FS,eq_FS_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4',name='FSgranLIP')
    FS_gran.V = '-110*mvolt+10*rand()*mvolt'
    FS_gran.h = '0+0.05*rand()'
    FS_gran.m = '0+0.05*rand()'
    if kainate=='low':
    #    FS_gran.J='5 * uA * cmeter ** -2' #article=code=35
        FS_gran.J=J_FSg
    elif kainate=='high':
        FS_gran.J='16 * uA * cmeter ** -2'
    
    SI_deep=NeuronGroup(N_SI,eq_SI_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4',name='SIdeepLIP')
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
    S_EgranFS=generate_syn(E_gran,FS,'IsynRS_LIP_gran','',gRSgFSs,0.125*ms,1*ms,0*mV)
    #S_EgranEgran=generate_syn(E_gran,E_gran,'IsynEgran','',0.4*usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)
    S_EgranEgran=generate_syn(E_gran,E_gran,'IsynRS_LIP_gran','',gRSgRSg,0.125*ms,1*ms,0*mV)
    #S_EgranFSgran=generate_syn(E_gran,FS_gran,'IsynEgran','',0.2*usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)
    S_EgranFSgran=generate_syn(E_gran,FS_gran,'IsynRS_LIP_gran','',gRSgFSg,0.125*ms,1*ms,0*mV)
    #S_EgranRS=generate_syn(E_gran,RS,'IsynEgran','',0.2*usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)
    #S_EgranRS=generate_syn(E_gran,RS,'IsynEgran','',1*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_EgranRS=generate_syn(E_gran,RS,'IsynRS_LIP_gran','',gRSgRSs,0.125*ms,1*ms,0*mV)
    
    S_EgranIB=generate_syn(E_gran,IB_ad,'IsynRS_LIP_gran','',0.212*usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)
    
    #From FS (granular layer) cells, normal timescale
    #S_FSgranEgran=generate_syn(FS_gran,E_gran,'IsynFSgran','',1* usiemens * cm **-2*FLee,0.25*ms,5*ms,-80*mV)
    S_FSgranEgran=generate_syn(FS_gran,E_gran,'IsynFS_LIP_gran','',gFSgRSg,0.25*ms,5*ms,-80*mV)
    S_FSgranFSgran=generate_syn(FS_gran,FS_gran,'IsynFS_LIP_gran','',gFSgFSg,0.25*ms,5*ms,-75*mV)
    #S_FSgranRS=generate_syn(FS_gran,RS,'IsynFSgran','',0.02* usiemens * cm **-2*FLee,0.25*ms,5*ms,-80*mV)
    S_FSgranRS=generate_syn(FS_gran,RS,'IsynFS_LIP_gran','',gFSgRSs,0.25*ms,5*ms,-80*mV)
    #S_FSgranRS=generate_syn(FS_gran,RS,'IsynFSgran','',0* msiemens * cm **-2,0.25*ms,5*ms,-80*mV)

#    #From FS (granular layer) cells, beta timescale
#    S_FSgranEgran=generate_syn(FS_gran,E_gran,'IsynFS_LIP_gran','',gFSgRSg,0.25*ms,20*ms,-80*mV)
#    S_FSgranFSgran=generate_syn(FS_gran,FS_gran,'IsynFS_LIP_gran','',gFSgFSg,0.25*ms,20*ms,-75*mV)
#    S_FSgranRS=generate_syn(FS_gran,RS,'IsynFS_LIP_gran','',gFSgRSs,0.25*ms,20*ms,-80*mV)
    
    #From IB cells
    #S_IBSIdeep=generate_syn(IB_axon,SI_deep,'IsynIB','',0.12* usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)
#    S_IBSIdeep=generate_syn(IB_axon,SI_deep,'IsynIB','',0.2* msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_IBSIdeep=generate_syn(IB_axon,SI_deep,'IsynIB_LIP','',0.01* msiemens * cm **-2,0.125*ms,1*ms,0*mV)
        
    
    #From deep SI cells    
    S_SIdeepIB=generate_syn(SI_deep,IB_bd,'IsynSI_LIP_deep','',10* msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
    #S_SIdeepFSgran=generate_syn(SI_deep,FS_gran,'IsynSIdeep','',0.4* usiemens * cm **-2*FLee,0.25*ms,20*ms,-80*mV)
    S_SIdeepFSgran=generate_syn(SI_deep,FS_gran,'IsynSI_LIP_deep','',gSIdFSg,0.25*ms,20*ms,-80*mV)
    
    
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
    
    beta2freq = 25*Hz
    
    if input_beta2_IB:
#        inputs_topdown=generate_spike_timing(N_IB,25*Hz,0*ms,end_time=3000*ms)
#        G_topdown = SpikeGeneratorGroup(N_IB, inputs_topdown[:,1], inputs_topdown[:,0]*second)
#        topdown_in=Synapses(G_topdown,IB_bd,on_pre='Vinp=Vhigh')
#        topdown_in.connect(j='i')
        
#        inputs_topdown3=generate_spike_timing(N_SI,50*Hz,0*ms,end_time=8000*ms)
        inputs_topdown3=generate_spike_timing(N_SI,beta2freq,0*ms,end_time=10000*ms)
        G_topdown3 = SpikeGeneratorGroup(N_SI, inputs_topdown3[:,1], inputs_topdown3[:,0]*second)
        topdown_in3=Synapses(G_topdown3,SI_deep,on_pre='Vinp=Vhigh')
        topdown_in3.connect(j='i')
        
    if input_beta2_RS:    
        RS.ginp_RS_good=4* msiemens * cm **-2
        RS.ginp_RS_bad=4* msiemens * cm **-2
        inputs_topdown2=generate_spike_timing(N_RS,beta2freq,0*ms,end_time=10000*ms)
        G_topdown2 = SpikeGeneratorGroup(N_RS, inputs_topdown2[:,1], inputs_topdown2[:,0]*second)
        topdown_in2=Synapses(G_topdown2,RS,on_pre='Vinp=Vhigh')
        topdown_in2.connect(j='i')
        
    if input_beta2_FS_SI:
        FS.ginp_FS_good=gFS
        FS.ginp_FS_bad=gFS
        inputs_lateral=generate_spike_timing(N_FS,beta2freq,0*ms,end_time=10000*ms)
        G_lateral = SpikeGeneratorGroup(N_FS, inputs_lateral[:,1], inputs_lateral[:,0]*second)
        lateral_in=Synapses(G_lateral,FS,on_pre='Vinp=Vhigh')
        lateral_in.connect(j='i')
        
        inputs_lateral2=generate_spike_timing(N_SI,beta2freq,0*ms,end_time=10000*ms)
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
        
        inputs_topdown=generate_spike_timing_theta(N_IB,beta2freq,0*ms,end_time=5100*ms)
        G_topdown = SpikeGeneratorGroup(N_IB, inputs_topdown[:,1], inputs_topdown[:,0]*second)
        topdown_in=Synapses(G_topdown,IB_bd,on_pre='Vinp=Vhigh')
        topdown_in.connect(j='i')
        
        #theta=4Hz
        inputs_topdown3=generate_spike_timing_theta(N_SI,beta2freq,0*ms,end_time=5100*ms)
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
    R5=SpikeMonitor(E_gran,record=True)
    R6=SpikeMonitor(FS_gran,record=True)
    R7=SpikeMonitor(SI_deep,record=True)
    
    V5=StateMonitor(E_gran,'V',record=True)
    V6=StateMonitor(FS_gran,'V',record=True)
    V7=StateMonitor(SI_deep,'V',record=True)
    
#    inpmon=StateMonitor(E_gran,'Iinp1',record=True)
    inpmon=StateMonitor(E_gran,'sinp',record=True)
    #graninpmon=StateMonitor(FS,'IsynEgran',record=[0])
    inpIBmon=StateMonitor(SI_deep,'Iapp',record=[0])
    
    targets=[]
    for i in range(len(all_synapses)):
        targets.append(all_synapses[i].target.name)
    #RSsupTargeting = [1 for x in targets if x=='RSsupLIP']
    RSsupTargeting = [i for i, x in enumerate(targets) if x=='RSsupLIP']
    
    #for i in range(len(RSsupTargeting)):
    #    all_monitors+=StateMonitor(all_synapses[RSsupTargeting[i]])
    #for property, value in vars(all_monitors[6].variables).items():
    #    print(property, ": ", value)

    all_neurons=all_neurons+(E_gran,FS_gran,SI_deep)+tuple(g_inputs)
    all_synapses=all_synapses+(S_EgranFS,S_EgranEgran,S_EgranFSgran,S_EgranRS,S_EgranIB,S_FSgranEgran,S_FSgranFSgran,S_FSgranRS,S_IBSIdeep,S_SIdeepIB,S_SIdeepFSgran)+tuple(syn_inputs)
    all_monitors=all_monitors+(R5,R6,R7,V5,V6,V7,inpmon,inpIBmon)
    return all_neurons, all_synapses, all_gap_junctions, all_monitors


def run_one_simulation(simu,path,index_var):
#    print(simu,len(simu))
    start_scope()
    close('all')

    runtime=2*second
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
    
    syn_cond,J,thal,theta_phase,index=simu
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
        ginpSIdeep=500* msiemens * cm **-2
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
    V1,V2,V3,R1,R2,R3,I1,I2,I3,V4,R4,I4s,I4a,I4ad,I4bd,R5,R6,R7,V5,V6,V7,inpmon,inpIBmon=all_monitors
    
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_gap_junctions)
#    net.add(all_monitors)
    net.add((V1,I1,R1,R2,R3,R4,R5,R6,R7,inpmon,inpIBmon))
#    net.add((V1,I1))
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
    
#    figure()
#    plot(R1.t,R1.i+140,'r.',label='RS cells')
#    plot(R2.t,R2.i+120,'m.',label='FS cells')
#    plot(R3.t,R3.i+100,'y.',label='SI cells')
#    plot(R5.t,R5.i+70,'g.',label='Granular RS')
#    plot(R6.t,R6.i+50,'c.',label='Granular FS')
#    plot(R4.t,R4.i+20,'b.',label='IB cells')
#    plot(R7.t,R7.i,'k.',label='Deep SI')
#    xlim(0,runtime/second)
#    legend(loc='upper left')
    
#    figure()
##    plot(inpmon.t,inpmon.Iinp1[0])
#    plot(inpmon.t,inpmon.sinp[0])
#    xlabel('Time (s)')
#    ylabel('Input synapse opening variable $s_{input}$')
#    tight_layout()
    
#    print(sum(inpmon.sinp[0]))
#    print(inpmon.t)
#    print(inpmon.Iinp1[0])
#    
#    print(inpIBmon.t)
#    print(inpIBmon.Iapp[0])
    
#    figure()
#    plot(inpIBmon.t,inpIBmon.Iapp[0])
    
    min_t=int(50*ms*100000*Hz)
#    min_t=int(150*ms*100000*Hz)
    LFP_V_RS = 1/N_RS*sum(V1.V,axis=0)#[min_t:]
    LFP_I_RS = 1/N_RS*sum(I1.Isyn, axis=0)
#    LFP_V_FS=1/N_FS*sum(V2.V,axis=0)[min_t:]
#    LFP_V_SI=1/N_SI62.903225806451616*sum(V3.V,axis=0)[min_t:]
#    LFP_V_IB=1/N_IB*sum(V4.V,axis=0)[min_t:]
#    LFP_V_RSg=1/N_FS*sum(V5.V,axis=0)[min_t:]
#    LFP_V_FSg=1/N_FS*sum(V6.V,axis=0)[min_t:]
#    LFP_V_SId=1/N_SI*sum(V7.V,axis=0)[min_t:]
    
    
    def zeros_ones_monitor(spikemon,record_dt,runtime):
        L=math.ceil(runtime/record_dt)
        zeros_ones=[0]*L
        for time in spikemon.t:
            zeros_ones[int(time/record_dt)]+=1
        return zeros_ones
    
    RS_spikes=array(zeros_ones_monitor(R1,defaultclock.dt,runtime))
    IB_spikes=array(zeros_ones_monitor(R4,defaultclock.dt,runtime))
    spike_shape = atleast_2d(shape(RS_spikes))
    print("Size of RS_spikes = "+str(spike_shape[0]))#+","+str(spike_shape[1])+"]")
    
    figure()
    subplot(311)
    plot(LFP_V_RS)
    subplot(312)
    plot(LFP_I_RS)
    subplot(313)
    plot(RS_spikes)

    
##    f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='spectrum')
#    f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS[min_t:], 100000,'flattop', scaling='density')
#    f,Spectrum_LFP_I_RS=signal.periodogram(LFP_I_RS[min_t:], 100000,'flattop', scaling='density')
###    f,Spectrum_LFP_V_FS=signal.periodogram(LFP_V_FS, 100000,'flattop', scaling='spectrum')
###    f,Spectrum_LFP_V_SI=signal.periodogram(LFP_V_SI, 100000,'flattop', scaling='spectrum')
###    f,Spectrum_LFP_V_IB=signal.periodogram(LFP_V_IB, 100000,'flattop', scaling='spectrum')
###    f,Spectrum_LFP_V_RSg=signal.periodogram(LFP_V_RSg, 100000,'flattop', scaling='spectrum')
###    f,Spectrum_LFP_V_FSg=signal.periodogram(LFP_V_FSg, 100000,'flattop', scaling='spectrum')
###    f,Spectrum_LFP_V_SId=signal.periodogram(LFP_V_SId, 100000,'flattop', scaling='spectrum')
#    
#    figure()
#    subplot(211)
#    plot(f,Spectrum_LFP_V_RS)
#    ylabel('Spectrum')
#    xlabel('Frequency (Hz)')
#    xlim(0,100)
#    subplot(212)
#    plot(f,Spectrum_LFP_I_RS)
#    ylabel('Spectrum')
#    xlabel('Frequency (Hz)')
#    xlim(0,100)
    

    
    def flipEnds(mat, end_length):
        beginning = mat[0:end_length, :]
        ending = mat[-end_length-1:-1, :]
        flipped = vstack((flipud(beginning), mat, flipud(ending)))
        return flipped
    
    def pctMean(mat, ax):
        diagMean = diag(nanmean(mat, axis=ax))
        matMean = ones(shape(mat))
        if ax == 0:
            matMean = matmul(ones(shape(mat)), diagMean)
        else:
            matMean = matmul(diagMean, ones(shape(mat)))
        normed = (mat - matMean)/matMean
        return normed
    
#    f, t, Sxx = signal.spectrogram(squeeze(LFPflip), fs, nperseg=25000, noverlap=20000)
#    pctSxx = pctMean(squeeze(Sxx[:, end_length:-end_length]), 0)
#    
#    figure()
#    #f, t, Sxx = signal.spectrogram(LFP_LIP, 100000*Hz,nperseg=30000,noverlap=25000)
#    pcolormesh(t, f, pctSxx)#, cmap='RdBu')#, shading='gouraud')
#    ylabel('Frequency [Hz]')
#    xlabel('Time [sec]')
#    ylim(0, 75)
       
    print(len(LFP_I_RS))
    savetxt('LIP_LFP_V_RS',LFP_V_RS,delimiter=',')
    savetxt('LIP_IB_spikes',IB_spikes,delimiter=',')
    savetxt('LIP_LFP_I_RS',LFP_I_RS,delimiter=',')
    savetxt('LIP_RS_spikes',RS_spikes,delimiter=',')
    f, t, Sxx = signal.spectrogram(LFP_I_RS, 100000*Hz, nperseg=12500, noverlap=10000)
    print(len(f),len(t),Sxx.shape)
    figure()
    pcolormesh(t, f, Sxx)#, cmap=)
    colorbar(format='%.1e')
    ylabel('Frequency (Hz)',fontsize=12)
    xlabel('Time (s)',fontsize=12)
    xticks(fontsize=12)
    yticks(fontsize=12)
    ylim(5,75)
    title('Power ($V^2$)',fontsize=12)

#     ind_tround_min=where(t>=0.375)[0][0]
#     ind_tround_max=where(t<=1.625)[0][-1]
#     print(ind_tround_min,ind_tround_max)
#     print(t[0],t[-1])
#     print(t[ind_tround_min],t[ind_tround_max])
#     ntpertheta=(ind_tround_max+1-ind_tround_min)//5
# #    Sxxfolded = Sxx[:,ind_tround_min:ind_tround_max+1].reshape((len(f), ntpertheta, 5), order='F')
#     Sxxfolded = Sxx[:,ind_tround_min:ind_tround_max].reshape((len(f), ntpertheta, 5), order='F')
#     SxxmeanTheta = nanmean(Sxxfolded, axis = 2)
    
#     figure()
#     pcolormesh(t[:ntpertheta]-t[0], f, hstack((SxxmeanTheta[:,ntpertheta//2:],SxxmeanTheta[:,:ntpertheta//2])))#, cmap=) 
#     colorbar(format='%.1e')
#     ylabel('Frequency (Hz)',fontsize=12)
#     xlabel('Time (s)',fontsize=12)
#     xticks(fontsize=12)
#     yticks(fontsize=12)
#     ylim(0,75)
#     title('Power ($V^2$)',fontsize=12)
    
#     SxxpctTheta = pctMean(hstack((SxxmeanTheta[:,ntpertheta//2:],SxxmeanTheta[:,:ntpertheta//2])), 1)
#     figure()
#     pcolormesh(t[:ntpertheta]-t[0], f, SxxpctTheta)
#     colorbar(format='%.1e')
#     ylabel('Frequency (Hz)',fontsize=12)
#     xlabel('Time (s)',fontsize=12)
#     xticks(fontsize=12)
#     yticks(fontsize=12)
#     ylim(0,75)    

    record_dt=1/100000*second#1/512*second
    t=int(0.3*second/record_dt) #t_debut
    L=int(2*second/record_dt)
    fs = 1/record_dt
    
    end_length = 5000
    LFPflip = flipEnds(LFP_I_RS[:, None], end_length)#transpose(atleast_2d(LFP_V_RS)), end_length)
    
    beta1 = signal.cwt(squeeze(LFPflip), signal.morlet2, fs/(10*Hz))
    beta1 = beta1[:, end_length:-end_length]
    beta1_phase = angle(beta1)
    mrv = mean(multiply(exp(1j*beta1_phase), IB_spikes))
    
    rbins = linspace(0, max(IB_spikes), 30)
    abins = linspace(-pi, pi, 60)
    hist, _, _ = histogram2d(beta1_phase, IB_spikes, bins=(abins,rbins))
    A, R = meshgrid(abins, rbins)
    fig, ax = subplots(subplot_kw = dict(projection = "polar"))
    pc = ax.pcolormesh(A, R, hist.T, cmap="magma_r")
    
    #freq = logspace(0, 2, 50)*Hz
    freq = linspace(1/second, 100*Hz, 100)
    #widths = logspace(log10(3),log10(30),50)*fs/(2*freq*pi)
    #widths = around(linspace(3,7,100)*fs/freq) 
    widths = fs/freq #linspace(3, 30, 100)*fs/(2*freq*pi)
    
    LFPflip = flipEnds(RS_spikes[:, None], end_length)
    
    CWT = signal.cwt(squeeze(LFPflip), signal.morlet2, widths)#, w=6)
    CWT = CWT[:, end_length:-end_length]
    # CWTmean = diag(mean(CWT, axis=1))
    # mean_mat = matmul(CWTmean, ones(shape(CWT)))
    # CWTpct = (CWT - mean_mat)/mean_mat
    CWTpct = pctMean(absolute(CWT), 1)
    
    figure()
    #f, t, Sxx = signal.spectrogram(LFP_LIP, 100000*Hz,nperseg=30000,noverlap=25000)
    pcolormesh(V1.t, freq, absolute(CWT))#, cmap='RdBu')#, shading='gouraud')
    ylabel('Frequency [Hz]')
    xlabel('Time [sec]')
    ylim(10, 75)
    
    figure()
    #f, t, Sxx = signal.spectrogram(LFP_LIP, 100000*Hz,nperseg=30000,noverlap=25000)
    pcolormesh(V1.t, freq, CWTpct)#, cmap='RdBu')#, shading='gouraud')
    ylabel('Frequency [Hz]')
    xlabel('Time [sec]')
    ylim(10, 75)
    
    if theta_phase=='mixed':
        #nanmat = empty((len(freq),5000))
        #nanmat[:] = NaN
        #CWTfilled = hstack((nanmat, CWT))
        figure()
        CWTfolded = CWT.reshape((len(freq), 25000, 8), order='F')
        CWTmeanTheta = nanmean(CWTfolded, axis = 2)
        CWTpctTheta = pctMean(absolute(CWTmeanTheta), 1)
        
        # for i in arange(8): 
        #     subplot(4, 2, i+1)
        #     pcolormesh(V1.t[0:25000], freq, absolute(CWTfolded[:,:,i]))
        
        figure()
        #f, t, Sxx = signal.spectrogram(LFP_LIP, 100000*Hz,nperseg=30000,noverlap=25000)
        pcolormesh(V1.t[0:25000], freq, absolute(CWTmeanTheta))#, cmap='RdBu')#, shading='gouraud')
        ylabel('Frequency [Hz]')
        xlabel('Time [sec]')
        ylim(10, 75)
        colorbar()
        
        figure()
        #f, t, Sxx = signal.spectrogram(LFP_LIP, 100000*Hz,nperseg=30000,noverlap=25000)
        pcolormesh(V1.t[0:25000], freq, CWTpctTheta)#, cmap='RdBu')#, shading='gouraud')
        ylabel('Frequency [Hz]',fontsize=12)
        xlabel('Time [sec]',fontsize=12)
        xticks(fontsize=12)
        yticks(fontsize=12)
        ylim(10, 75)
        colorbar()
#        savetxt('data.csv', CWTpctTheta, delimiter=',')
    
    
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
    
#    figure()
#    plot(R1.t,R1.i+140,'r.',label='RS cells')
#    plot(R2.t,R2.i+120,'b.',label='FS cells')
#    plot(R3.t,R3.i+100,'g.',label='SI cells')
#    plot(R5.t,R5.i+70,'.',label='Granular RS',color='C1')
#    plot(R6.t,R6.i+50,'c.',label='Granular FS')
#    plot(R4.t,R4.i+20,'m.',label='IB cells')
#    plot(R7.t,R7.i,'.',label='Deep SI',color='lime')
#    xlim(0,runtime/second)
#    ylim(0,220)
#    legend(loc='upper left')
#    xlabel('Time (s)')
#    ylabel('Neuron index')
#    
#    figure()
#    plot(R1.t,R1.i+140,'r.',label='RS cells')
#    plot(R2.t,R2.i+120,'b.',label='FS cells')
#    plot(R3.t,R3.i+100,'g.',label='SOM cells')
#    plot([0.2,runtime/second],[95,95],'k--')
#    plot(R5.t,R5.i+70,'r.')
#    plot(R6.t,R6.i+50,'b.')
#    plot([0.2,runtime/second],[45,45],'k--')
#    plot(R4.t,R4.i+20,'m.',label='IB cells')
#    plot(R7.t,R7.i,'g.')
#    xlim(0.2,runtime/second)
#    ylim(0,220)
##    legend(loc='upper left')
#    xlabel('Time (s)')
#    ylabel('Neuron index')
    
    
#    min_t=int(50*ms*100000*Hz)
#    LFP_V_RS=1/N_RS*sum(V1.V,axis=0)[min_t:]
#    LFP_V_FS=1/N_FS*sum(V2.V,axis=0)[min_t:]
#    LFP_V_SI=1/N_SI*sum(V3.V,axis=0)[min_t:]
#    LFP_V_IB=1/N_IB*sum(V4.V,axis=0)[min_t:]
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
#    new_path=path+str(index)
#    os.mkdir(new_path)
#
#    for n in get_fignums():
#        current_fig=figure(n)
#        current_fig.savefig(new_path+'/figure'+str(n)+'.png')
        
#    save_raster('LIP_RS',R1.i,R1.t,new_path)
#    save_raster('LIP_FS',R2.i,R2.t,new_path)
#    save_raster('LIP_SI',R3.i,R3.t,new_path)
#    save_raster('LIP_RS_gran',R5.i,R5.t,new_path)
#    save_raster('LIP_FS_gran',R6.i,R6.t,new_path)
#    save_raster('LIP_IB',R4.i,R4.t,new_path)
#    save_raster('LIP_SI_deep',R7.i,R7.t,new_path)
#    save_raster('LIP_RS_V',V1.V,V1.t,new_path)

    
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
    all_thal=[10* msiemens * cm **-2]
#    all_thal=[0* msiemens * cm **-2]
    all_theta=['mixed']
    #all_theta=['mixed','mixed','mixed','mixed','mixed']
    
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
    
    this_time=datetime.datetime.now()
    path="./results_"+this_time.strftime("%y-%m-%d_%H-%M-%S")
    os.mkdir(path)
        
    all_sim=list(product(all_syn_cond,all_J,all_thal,all_theta))
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
    
    # clear_cache('cython') 
