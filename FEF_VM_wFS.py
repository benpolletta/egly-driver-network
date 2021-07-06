# -*- coding: utf-8 -*-
"""
Created on July 1 2021

@author: benpolletta
"""

from brian2 import *

from scipy import signal
from cells.RS_FEF import *
from cells.FS_FEF import *
from cells.SI_FEF import *
from cells.VIP_FEF import *

from FEF_and_LIP_parallel import save_raster

import os
import sys

def generate_VM_wFS(theta_phase,IappSI,gFSSI,gVIPSI,runtime):
    
    N_VIP,N_RS_gran,N_FS_gran,N_SI_gran=20,20,20,20
    
    if theta_phase=='bad':
        ginp_IB=0* msiemens * cm **-2
        input_beta2_RS=False
        input_beta2_FS_SI=True
        input_thalamus_gran=True
        gFS=0* msiemens * cm **-2
        thal_cond=3* msiemens * cm **-2
#        thal_cond=10* msiemens * cm **-2
        kainate='low'
        
    if theta_phase=='good' or theta_phase=='mixed':
        ginp_IB=10* msiemens * cm **-2#*2754.660086037123/139.46773954954165
#        ginp_IB=2* msiemens * cm **-2#*2754.660086037123/12782.0904181147
#        ginp_IB=0* msiemens * cm **-2
        input_beta2_RS=False
        input_beta2_FS_SI=False
        input_thalamus_gran=True
#        thal_cond=3* msiemens * cm **-2
#        thal_cond=0* msiemens * cm **-2
#        thal_cond=15* msiemens * cm **-2
#        thal_cond=4* msiemens * cm **-2
        thal_cond=5* msiemens * cm **-2
        kainate='low'
        
    
    prefs.codegen.target = 'numpy'
    
    defaultclock.dt = 0.01*ms
    
    #Single column network
    #N_VIP,N_RS_gran,N_FS_gran,N_VIP_gran=20,20,20,20
    
    ##Define neuron groups
    E_gran=NeuronGroup(N_RS_gran,eq_RS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    E_gran.V = '-70*mvolt+10*rand()*mvolt'
    E_gran.h = '0+0.05*rand()'
    E_gran.m = '0+0.05*rand()'
    E_gran.mAR = '0.035+0.025*rand()'
#    E_gran.J='30 * uA * cmeter ** -2'  #article SI=25, code=1
#    E_gran.J='20 * uA * cmeter ** -2'  #article SI=25, code=1
    E_gran.J=str(IappSI)+'* uA * cmeter ** -2'  #article SI=25, code=1
    
    FS_gran=NeuronGroup(N_FS_gran,eq_FS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    FS_gran.V = '-110*mvolt+10*rand()*mvolt'
    FS_gran.h = '0+0.05*rand()'
    FS_gran.m = '0+0.05*rand()'
#    FS_gran.J='0 * uA * cmeter ** -2' #article=code=35
    FS_gran.J='5 * uA * cmeter ** -2' #article=code=35
    
    SI_gran=NeuronGroup(N_SI_gran,eq_SI_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    SI_gran.V = '-110*mvolt+10*rand()*mvolt'
    SI_gran.h = '0+0.05*rand()'
    SI_gran.m = '0+0.05*rand()'
#    SI_gran.J='5 * uA * cmeter ** -2' #article=code=35
    SI_gran.J='0 * uA * cmeter ** -2' 

    VIP=NeuronGroup(N_VIP,eq_VIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    VIP.V = '-63*mvolt'
    VIP.Iapp='0 * uA * cmeter ** -2'
#    VIP.Iapp='8 * uA * cmeter ** -2'
    
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
    
    
    #From RS cells
    S_EgranEgran=generate_syn(E_gran,E_gran,'IsynRS_FEF_VM','',0.6*msiemens * cm **-2,0.125*ms,1*ms,0*mV) #0.4 #0.6
    S_EgranFSgran=generate_syn(E_gran,FS_gran,'IsynRS_FEF_VM','',0.2*msiemens * cm **-2,0.125*ms,1*ms,0*mV) #0.2
    S_EgranSIgran=generate_syn(E_gran,SI_gran,'IsynRS_FEF_VM','',0.5*msiemens * cm **-2,0.125*ms,1*ms,0*mV) #0.6 #0.5
    
    #From FS cells
    S_FSgranFSgran=generate_syn(FS_gran,FS_gran,'IsynFS_FEF_VM','',0.2*msiemens * cm **-2,0.25*ms,5*ms,-80*mV) #0.2
    S_FSgranEgran=generate_syn(FS_gran,E_gran,'IsynFS_FEF_VM','',0.2*msiemens * cm **-2,0.25*ms,5*ms,-80*mV) #0.2
    S_FSgranSIgran=generate_syn(FS_gran,SI_gran,'IsynFS_FEF_VM','',1*msiemens * cm **-2,0.25*ms,5*ms,-80*mV) #0.5

    #From SOM cells, normal timescale
    S_SIgranEgran=generate_syn(SI_gran,E_gran,'IsynSI_FEF_VM','',0.5*msiemens * cm **-2,0.25*ms,20*ms,-80*mV) #0.35 #0.5
    S_SIgranFSgran=generate_syn(SI_gran,FS_gran,'IsynSI_FEF_VM','',0.2*msiemens * cm **-2,0.25*ms,20*ms,-80*mV) #0.2
    S_SIgranSIgran=generate_syn(SI_gran,SI_gran,'IsynSI_FEF_VM','i',0* msiemens * cm **-2,0.25*ms,20*ms,-75*mV) #1 #0.2
    
#    #From SOM cells, gamma timescale
    S_VIPSIgran=generate_syn(VIP,SI_gran,'IsynSI2_FEF_VM','',1*msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
    
    eq_gap='''_post=g_i*(V_post-V_pre) : amp * meter ** -2 (summed)
        g_i : siemens * meter**-2
    '''
    
    gap_SISI=Synapses(SI_gran,SI_gran,model='Igap'+eq_gap,method='exact')
    gap_SISI.connect()
    gap_SISI.g_i=0* msiemens * cm **-2
        
    
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


    if theta_phase=='good':
        VIP.ginp_VIP_good=ginp_IB
        VIP.ginp_VIP_bad=ginp_IB
    elif theta_phase=='mixed':
        VIP.ginp_VIP_good=ginp_IB 
        VIP.ginp_VIP_bad=ginp_IB
#        print(VIP.ginp_VIP_good)
    fIB=13*Hz
#    fIB=30*Hz
    inputs_topdown3=generate_spike_timing(N_VIP,fIB,0*ms,end_time=3000*ms)
    
    Theta_freq=4*Hz
    Theta_pd=1/Theta_freq
    if theta_phase=='mixed':
        t0=0*ms
        t1=Theta_pd/2
        inputs_topdown3=generate_spike_timing(N_VIP,fIB,t0,end_time=t1)
        while t0+Theta_pd<runtime:
            t0,t1=t0+Theta_pd,t1+Theta_pd
            inputs_topdown3=vstack((inputs_topdown3,generate_spike_timing(N_VIP,fIB,t0,end_time=t1)))
    
    
    G_topdown3 = SpikeGeneratorGroup(N_VIP, inputs_topdown3[:,1], inputs_topdown3[:,0]*second)
    topdown_in3=Synapses(G_topdown3,VIP,on_pre='Vinp=Vhigh')
    topdown_in3.connect(j='i')
        
    if input_thalamus_gran:    
        SI_gran.ginp_SI=thal_cond
        E_gran.ginp_RS=thal_cond
        if theta_phase=='good' or theta_phase=='mixed':
            fLIP=50*Hz
#            fLIP=12*Hz #test, if LIP hasn't switched to its good phase activity
        else :
            fLIP=13*Hz
#        print(fLIP)
        bottomup=generate_spike_timing(N_SI_gran,fLIP,0*ms,end_time=2100*ms)
        #theta=4Hz
        if theta_phase=='mixed':
            t0=0*ms
            t1=Theta_pd/2
            fLIP=50*Hz
#            fLIP=12*Hz
            bottomup=generate_spike_timing(N_SI_gran,fLIP,t0,end_time=t1)
            while t0+Theta_pd<runtime:
                t0,t1=t0+Theta_pd/2,t1+Theta_pd/2
                fLIP=50*Hz*int(fLIP==13*Hz)+13*Hz*int(fLIP==50*Hz)
                bottomup=vstack((bottomup,generate_spike_timing(N_SI_gran,fLIP,t0,end_time=t1)))
                      
                 
        Poisson_input = SpikeGeneratorGroup(N_SI_gran, bottomup[:,1], bottomup[:,0]*second)
        bottomup_in=Synapses(Poisson_input,SI_gran,on_pre='Vinp=Vhigh')
        bottomup_in.connect(j='i')
        Poisson_input2 = SpikeGeneratorGroup(N_RS_gran, bottomup[:,1], bottomup[:,0]*second)
        bottomup_in2=Synapses(Poisson_input2,E_gran,on_pre='Vinp=Vhigh')
        bottomup_in2.connect(j='i')
    
    #Define monitors and run network :
    R5=SpikeMonitor(E_gran,record=True)
    R6=SpikeMonitor(SI_gran,record=True)
    R7=SpikeMonitor(VIP,record=True)
    R8=SpikeMonitor(FS_gran,record=True)
    BUinput=SpikeMonitor(Poisson_input,record=True)
    #TDinput=SpikeMonitor(G_topdown3,record=True)
    
    inpmon=StateMonitor(E_gran,'Iinp1',record=True)
    #graninpmon=StateMonitor(FS,'IsynEgran',record=[0])
    #inpIBmon=StateMonitor(IB_bd,'Iapp',record=[0])
    
    V_RS=StateMonitor(E_gran,'V',record=True)
    V_FS=StateMonitor(SI_gran,'V',record=True)
    V_SI=StateMonitor(VIP,'V',record=True)
    
    all_neurons=VIP,E_gran,FS_gran,SI_gran,G_topdown3,Poisson_input,Poisson_input2
    all_synapses=S_EgranEgran,S_EgranFSgran,S_EgranSIgran,S_FSgranEgran,S_FSgranFSgran,S_FSgranSIgran,S_SIgranEgran,S_SIgranFSgran,S_SIgranSIgran,S_VIPSIgran,gap_SISI,topdown_in3,bottomup_in,bottomup_in2
    all_monitors=R5,R6,R7,R8,V_RS,V_FS,V_SI
#    all_monitors=R5,R6,R7,V_RS,V_FS,V_SI,inpmon
    
    return all_neurons,all_synapses,all_monitors



if __name__=='__main__':
    close('all')
    start_scope()    
    
    prefs.codegen.target = 'numpy'
    defaultclock.dt = 0.01*ms
    
    theta_phase='mixed' #'good' or 'bad' or 'mixed'
    runtime=2*second
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp 
    Vhigh=0*mV
    Vlow=-80*mV
    ginp=0* msiemens * cm **-2
    taurinp2=2*ms
    taudinp2=10*ms
    tauinp2=taudinp2   
    
    gFSSI=[]
    gVIPSI=[]
    # J=[]
    for step in range(11):
        gFSSI.append(1+step*3/10)
        gVIPSI.append(1-step/10)
    #    J.append(-10+3*step)
        
    params=transpose([tile(gFSSI, len(gVIPSI)), repeat(gVIPSI, len(gFSSI))])
    # params=transpose([tile(J, len(params)), repeat(params, len(J))])
    
    for pset in range(len(params)):
        
        all_neurons,all_synapses,all_monitors=generate_VM_wFS(theta_phase,0,params[pset][0],params[pset][1],runtime)    
        
        name = 'FEF_VM_wFS_gFS'+str(params[pset][0])+'_gVS'+str(params[pset][1])
        sim_dir = 'sims/'+name
        
        if not os.path.exists(sim_dir):
            os.mkdir(sim_dir)
        
        net=Network()
        net.add(all_neurons)
        net.add(all_synapses)
        net.add(all_monitors)
        
        prefs.codegen.target = 'cython' #cython=faster, numpy = default python
        
        net.run(runtime,report='text',report_period=300*second)
    
        R5,R6,R7,R8,V_RS,V_FS,V_SI=all_monitors
    #    R5,    R6,R7,V_RS,V_FS,V_SI,inpmon=all_monitors
        
        save_raster('RS',R5.i,R5.t,sim_dir)
        save_raster('SI',R6.i,R6.t,sim_dir)
        save_raster('VIP',R7.i,R7.t,sim_dir)
        save_raster('FS',R8.i,R8.t,sim_dir)
        
        figure()
        plot(R7.t,R7.i+0,'k.',label='VIP')
        plot(R5.t,R5.i+20,'r.',label='RS')
        plot(R6.t,R6.i+40,'g.',label='SOM')
        plot(R8.t,R8.i+60,'b.',label='FS')
        xlim(0,runtime/second)
    #    legend(loc='upper left')
        xlabel('Time (s)')
        ylabel('Neuron index')
        ylim(-1,81)
        savefig(sim_dir+'/raster.png')
        
        min_t=int(50*ms*100000*Hz)
        LFP_V_RS=1/20*sum(V_RS.V,axis=0)[min_t:]
        LFP_V_FS=1/20*sum(V_FS.V,axis=0)[min_t:]
        
        f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='spectrum')
        f,Spectrum_LFP_V_FS=signal.periodogram(LFP_V_FS, 100000,'flattop', scaling='spectrum')
        
        figure()
        subplot(221)
        plot((V_RS.t/second)[min_t:],LFP_V_RS)
        ylabel('LFP')
        title('gran RS cell')
        subplot(223)
        plot((V_FS.t/second)[min_t:],LFP_V_FS)
        ylabel('LFP')
        title('gran SOM cell')
        plot()
        
        subplot(222)
        plot(f,Spectrum_LFP_V_RS)
        ylabel('Spectrum')
        yticks([],[])
        xlim(0,100)
        title('gran RS cell')
        subplot(224)
        plot(f,Spectrum_LFP_V_FS)
        ylabel('Spectrum')
        yticks([],[])
        xlim(0,100)
        title('gran SOM cell')
        savefig(sim_dir+'/RSFS_V.png')
        
        figure()
        plot(f,Spectrum_LFP_V_RS)
        ylabel('Spectrum')
        xlabel('Frequency (Hz)')
        xlim(0,50)
        
        f, t, Sxx = signal.spectrogram(LFP_V_RS, 100000*Hz,nperseg=20000,noverlap=15000)
        figure()
        pcolormesh(t, f, Sxx)#, shading='gouraud')
        ylabel('Frequency [Hz]')
        xlabel('Time [sec]')
        ylim(0,50)
    
        f, t, Sxx = signal.spectrogram(LFP_V_RS, 100000*Hz,nperseg=20000,noverlap=15000)
        figure()
        pcolormesh(t, f, Sxx)#, cmap=)
        colorbar(format='%.1e')
        ylabel('Frequency (Hz)')
        xlabel('Time (s)')
        ylim(0,45)
        title('Power ($V^2$)')
        savefig(sim_dir+'/spec.png')
        
        close('all')
        
    clear_cache('cython')