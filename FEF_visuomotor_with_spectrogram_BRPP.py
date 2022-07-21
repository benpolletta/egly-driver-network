# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 14:28:44 2019

@author: aaussel
"""

from brian2 import *

from scipy import signal
from cells.RS_FEF_VM import *
from cells.FS_FEF import *
from cells.SI_FEF_more_h import *
from cells.VIP_FEF import *

def zeros_ones_monitor(spikemon,record_dt,runtime):
    L=int(runtime/record_dt)
    zeros_ones=[0]*L
    for time in spikemon.t:
        zeros_ones[int(time/record_dt)]+=1
    return zeros_ones

def generate_deepSI_and_gran_layers(theta_phase,N_SI,N_RS_gran,N_SI_gran,runtime):
    
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
#        thal_cond=1* msiemens * cm **-2
        thal_cond=5* msiemens * cm **-2
        kainate='low'
        

    
    prefs.codegen.target = 'numpy'
    
    defaultclock.dt = 0.01*ms
    
    #Single column network
    
#    noise = TimedArray(100* uA * cmeter ** -2* rand(200000,20), dt=defaultclock.dt)
    
    ##Define neuron groups
    E_gran=NeuronGroup(N_RS_gran,eq_RS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    E_gran.V = '-70*mvolt+10*rand()*mvolt'
    E_gran.h = '0+0.05*rand()'
    E_gran.m = '0+0.05*rand()'
    E_gran.mAR = '0.035+0.025*rand()'
#    E_gran.J='30 * uA * cmeter ** -2'  #article SI=25, code=1
#    E_gran.J='15 * uA * cmeter ** -2'  #article SI=25, code=1
#    E_gran.J='10 * uA * cmeter ** -2'  #article SI=25, code=1
    #0   
#    E_gran.J='75 * uA * cmeter ** -2'
    E_gran.J_fixed='45 * uA * cmeter ** -2'
#    E_gran.J_fixed='80* uA * cmeter ** -2'
#    E_gran.J='30 * uA * cmeter ** -2+300*uamp*cmeter**-2*rand()'
    
    SI_gran=NeuronGroup(N_SI_gran,eq_SI_FEF_more_h,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    SI_gran.V = '-110*mvolt+10*rand()*mvolt'
    SI_gran.h = '0+0.05*rand()'
    SI_gran.m = '0+0.05*rand()'
#    SI_gran.J='15 * uA * cmeter ** -2' #article=code=35
#    SI_gran.J='0 * uA * cmeter ** -2' #article=code=35
     #-30
    SI_gran.J='40 * uA * cmeter ** -2'
#    SI_gran.J='50 * uA * cmeter ** -2'
#    SI_gran.J='75 * uA * cmeter ** -2'
     
#    SI_deep=NeuronGroup(N_SI,eq_SIdeep,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
#    SI_deep.V = '-100*mvolt+10*rand()*mvolt'
#    SI_deep.h = '0+0.05*rand()'
#    SI_deep.m = '0+0.05*rand()'
#    SI_deep.mAR = '0.02+0.04*rand()'
#    SI_deep.J='35* uA * cmeter ** -2' #article SI=50, code=35, Mark = 45

    SI_deep=NeuronGroup(N_SI,eq_VIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    SI_deep.V = '-63*mvolt'
    SI_deep.Iapp='0 * uA * cmeter ** -2'
#    SI_deep.Iapp='8 * uA * cmeter ** -2'
    
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
    #S_EgranEgran=generate_syn(E_gran,E_gran,'IsynEgran','',0.4*usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)
    #S_EgranEgran=generate_syn(E_gran,E_gran,'IsynEgran','',1/160*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_EgranEgran=generate_syn(E_gran,E_gran,'IsynRS_FEF_VM','',0.6*msiemens * cm **-2,0.125*ms,1*ms,0*mV) #0.4
#    S_EgranEgran=generate_syn(E_gran,E_gran,'IsynRS_FEF_VM','',1*msiemens * cm **-2,0.125*ms,1*ms,0*mV) #0.4
    #S_EgranFSgran=generate_syn(E_gran,SI_gran,'IsynEgran','',0.2*usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)
    S_EgranFSgran=generate_syn(E_gran,SI_gran,'IsynRS_FEF_VM','',0.5*msiemens * cm **-2,0.125*ms,1*ms,0*mV) #0.6
    #S_EgranRS=generate_syn(E_gran,RS,'IsynEgran','',0.2*usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)

    #From SOM cells, normal timescale
    #S_FSgranEgran=generate_syn(SI_gran,E_gran,'IsynFSgran','',1* usiemens * cm **-2*FLee,0.25*ms,5*ms,-80*mV)
    #S_FSgranEgran=generate_syn(SI_gran,E_gran,'IsynSI_FEF_VM','',0.6* msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
#    S_FSgranEgran=generate_syn(SI_gran,E_gran,'IsynSI_FEF_VM','',0.5*msiemens * cm **-2,0.25*ms,20*ms,-80*mV) #0.35
    S_FSgranEgran=generate_syn(SI_gran,E_gran,'IsynSI_FEF_VM','',0.8*msiemens * cm **-2,0.25*ms,20*ms,-80*mV) #0.35
    
    
    #S_FSgranFSgran=generate_syn(SI_gran,SI_gran,'IsynFSgran','',0.1* usiemens * cm **-2*FLee,0.25*ms,5*ms,-75*mV)
    S_FSgranFSgran=generate_syn(SI_gran,SI_gran,'IsynSI_FEF_VM','',0.2* msiemens * cm **-2,0.25*ms,20*ms,-75*mV) #1
    
#    #From SOM cells, gamma timescale
#    S_FSgranEgran=generate_syn(SI_gran,E_gran,'IsynSI_FEF_VM','',0.5*msiemens * cm **-2,0.25*ms,5*ms,-80*mV) #0.35
#    S_FSgranFSgran=generate_syn(SI_gran,SI_gran,'IsynSI_FEF_VM','',0.2* msiemens * cm **-2,0.25*ms,5*ms,-75*mV) #1
    
    
    #From VIP cells    
    #S_SIdeepFSgran=generate_syn(SI_deep,SI_gran,'IsynSIdeep','',0.4* usiemens * cm **-2*FLee,0.25*ms,20*ms,-80*mV)
#    S_SIdeepFSgran=generate_syn(SI_deep,SI_gran,'IsynSI2_FEF_VM','',1*msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
    S_SIdeepFSgran=generate_syn(SI_deep,SI_gran,'IsynSI2_FEF_VM','',0*msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
        
    
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
        SI_deep.ginp_VIP_good=ginp_IB
        SI_deep.ginp_VIP_bad=ginp_IB
    elif theta_phase=='mixed':
        SI_deep.ginp_VIP_good=ginp_IB 
        SI_deep.ginp_VIP_bad=ginp_IB
#        print(SI_deep.ginp_VIP_good)
    fIB=13*Hz
#    fIB=50*Hz
    inputs_topdown3=generate_spike_timing(N_SI,fIB,0*ms,end_time=3000*ms)
    
    #Theta=4Hz
    if theta_phase=='mixed':
        t0=0*ms
        t1=125*ms
        inputs_topdown3=generate_spike_timing(N_SI,fIB,t0,end_time=t1)
        while t0+250*ms<runtime:
            t0,t1=t0+250*ms,t1+250*ms
            inputs_topdown3=vstack((inputs_topdown3,generate_spike_timing(N_SI,fIB,t0,end_time=t1)))
            
#    #Theta=8Hz
#    if theta_phase=='mixed':
#        t0=0*ms
#        t1=62.5*ms
#        inputs_topdown3=generate_spike_timing(N_SI,fIB,t0,end_time=t1)
#        while t0+125*ms<runtime:
#            t0,t1=t0+125*ms,t1+125*ms
#            inputs_topdown3=vstack((inputs_topdown3,generate_spike_timing(N_SI,fIB,t0,end_time=t1)))

#    #Theta=2Hz
#    if theta_phase=='mixed':
#        t0=0*ms
#        t1=250*ms
#        inputs_topdown3=generate_spike_timing(N_SI,fIB,t0,end_time=t1)
#        while t0+500*ms<runtime:
#            t0,t1=t0+500*ms,t1+500*ms
#            inputs_topdown3=vstack((inputs_topdown3,generate_spike_timing(N_SI,fIB,t0,end_time=t1)))
                        
    
    G_topdown3 = SpikeGeneratorGroup(N_SI, inputs_topdown3[:,1], inputs_topdown3[:,0]*second)
#    topdown_in3=Synapses(G_topdown3,SI_deep,on_pre='Vinp=Vhigh')
#    topdown_in3.connect(j='i')
    
    topdown_in3=Synapses(G_topdown3,E_gran,on_pre='Vinp2=Vhigh')
    topdown_in3.connect(j='i')
    
    topdown_in4=Synapses(G_topdown3,SI_deep,on_pre='Vinp=Vhigh')
    topdown_in4.connect(j='i')
        
#    if input_beta2_RS:    
#        RS.ginp_RS=4* msiemens * cm **-2
#        inputs_topdown2=generate_spike_timing(N_RS,25*Hz,0*ms,end_time=2100*ms)
#        G_topdown2 = SpikeGeneratorGroup(N_RS, inputs_topdown2[:,1], inputs_topdown2[:,0]*second)
#        topdown_in2=Synapses(G_topdown2,RS,on_pre='Vinp=Vhigh')
#        topdown_in2.connect(j='i')
        
    #if input_beta2_FS_SI:
    #    FS.ginp_FS=gFS
    #    inputs_lateral=generate_spike_timing(N_FS,25*Hz,0*ms,end_time=2100*ms)
    #    G_lateral = SpikeGeneratorGroup(N_FS, inputs_lateral[:,1], inputs_lateral[:,0]*second)
    #    lateral_in=Synapses(G_lateral,FS,on_pre='Vinp=Vhigh')
    #    lateral_in.connect(j='i')
    #    
    #    inputs_lateral2=generate_spike_timing(N_SI,25*Hz,0*ms,end_time=2100*ms)
    #    G_lateral2 = SpikeGeneratorGroup(N_SI, inputs_lateral2[:,1], inputs_lateral2[:,0]*second)
    #    lateral_in2=Synapses(G_lateral2,SI,on_pre='Vinp=Vhigh')
    #    lateral_in2.connect(j='i')
        
    if input_thalamus_gran:    
        SI_gran.ginp_SI=thal_cond
        E_gran.ginp_RS=thal_cond
        E_gran.ginp_RS2=5* msiemens * cm **-2
#        E_gran.ginp_RS2=1* msiemens * cm **-2
    #    SI_gran.ginp_FS=thal_cond
#        Poisson_input = PoissonGroup(N_SI_gran,100*Hz)
#        bottomup_in = Synapses(Poisson_input,SI_gran, on_pre='Vinp=Vhigh')
#        bottomup_in.connect(j='i')
#        
#        Poisson_input2 = PoissonGroup(N_RS_gran,100*Hz)
#        bottomup_in2 = Synapses(Poisson_input2,E_gran, on_pre='Vinp=Vhigh')
#        bottomup_in2.connect(j='i')
    #    print(bottomup_in,bottomup_in2)
        if theta_phase=='good' or theta_phase=='mixed':
            fLIP=50*Hz
#            fLIP=13*Hz #test, if LIP hasn't switched to its good phase activity
        else :
            fLIP=13*Hz
#        print(fLIP)
        bottomup=generate_spike_timing(N_SI_gran,fLIP,0*ms,end_time=2100*ms)
        #theta=4Hz
        if theta_phase=='mixed':
            t0=0*ms
            t1=125*ms
            fLIP=50*Hz
#            fLIP=12*Hz
            bottomup=generate_spike_timing(N_SI_gran,fLIP,t0,end_time=t1)
#            while t0+250*ms<runtime:
#                t0,t1=t0+250*ms,t1+250*ms
#                fLIP=50*Hz*int(fLIP==13*Hz)+13*Hz*int(fLIP==50*Hz)
#                bottomup=vstack((bottomup,generate_spike_timing(N_SI_gran,fLIP,t0,end_time=t1)))
            while t0+250*ms<runtime:
                t0,t1=t0+125*ms,t1+125*ms
                fLIP=50*Hz*int(fLIP==13*Hz)+13*Hz*int(fLIP==50*Hz)
#                fLIP=12*Hz*int(fLIP==13*Hz)+13*Hz*int(fLIP==12*Hz)
#                print(fLIP)
                bottomup=vstack((bottomup,generate_spike_timing(N_SI_gran,fLIP,t0,end_time=t1)))
            
#        #theta=8Hz
#        if theta_phase=='mixed':
#            t0=0*ms
#            t1=62.5*ms
#            fLIP=50*Hz
#            bottomup=generate_spike_timing(N_SI_gran,fLIP,t0,end_time=t1)
#            while t0+125*ms<runtime:
#                t0,t1=t0+62.5*ms,t1+62.5*ms
#                fLIP=50*Hz*int(fLIP==13*Hz)+13*Hz*int(fLIP==50*Hz)
#                bottomup=vstack((bottomup,generate_spike_timing(N_SI_gran,fLIP,t0,end_time=t1)))
                      
#        #theta=2Hz
#        if theta_phase=='mixed':
#            t0=0*ms
#            t1=250*ms
#            fLIP=50*Hz
#            bottomup=generate_spike_timing(N_SI_gran,fLIP,t0,end_time=t1)
#            while t0+500*ms<runtime:
#                t0,t1=t0+250*ms,t1+250*ms
#                fLIP=50*Hz*int(fLIP==13*Hz)+13*Hz*int(fLIP==50*Hz)
#                bottomup=vstack((bottomup,generate_spike_timing(N_SI_gran,fLIP,t0,end_time=t1)))
                      
                 
        Poisson_input = SpikeGeneratorGroup(N_SI_gran, bottomup[:,1], bottomup[:,0]*second)
        bottomup_in=Synapses(Poisson_input,SI_gran,on_pre='Vinp=Vhigh')
        bottomup_in.connect(j='i')
        Poisson_input2 = SpikeGeneratorGroup(N_RS_gran, bottomup[:,1], bottomup[:,0]*second)
        bottomup_in2=Synapses(Poisson_input2,E_gran,on_pre='Vinp=Vhigh')
        bottomup_in2.connect(j='i')
    
    #Define monitors and run network :
    R5=SpikeMonitor(E_gran,record=True)
    R6=SpikeMonitor(SI_gran,record=True)
    R7=SpikeMonitor(SI_deep,record=True)
    
#    inpmon=StateMonitor(E_gran,'Iinp2',record=True)
    inpmon=StateMonitor(E_gran,'J',record=True)
    #graninpmon=StateMonitor(FS,'IsynEgran',record=[0])
    #inpIBmon=StateMonitor(IB_bd,'Iapp',record=[0])
    
    V_RS=StateMonitor(E_gran,'V',record=True)
    V_FS=StateMonitor(SI_gran,'V',record=True)
    V_SI=StateMonitor(SI_deep,'V',record=True)
    
    I_RS=StateMonitor(E_gran,'Isyn',record=True)
    I_FS=StateMonitor(SI_gran,'Isyn',record=True)
    I_SI=StateMonitor(SI_deep,'Isyn',record=True)
    
    all_neurons=SI_deep,E_gran,SI_gran,G_topdown3,Poisson_input,Poisson_input2
#    all_synapses=S_EgranEgran,S_EgranFSgran,S_FSgranEgran,S_FSgranFSgran,S_SIdeepFSgran,topdown_in3,bottomup_in,bottomup_in2
    #no SOM to SOM connections
#    all_synapses=S_EgranEgran,S_EgranFSgran,S_FSgranEgran,S_SIdeepFSgran,topdown_in3,bottomup_in,bottomup_in2
    all_synapses=S_EgranEgran,S_EgranFSgran,S_FSgranEgran,topdown_in3,bottomup_in,bottomup_in2,topdown_in4
#    all_monitors=R5,R6,R7,V_RS,V_FS,V_SI
    all_monitors=R5,R6,R7,V_RS,V_FS,V_SI,inpmon,I_RS,I_FS,I_SI
    
    return all_neurons,all_synapses,all_monitors



if __name__=='__main__':
    close('all')
    start_scope()    
    
    prefs.codegen.target = 'numpy'
    defaultclock.dt = 0.01*ms
    
    runtime=2*second
    
    noise = TimedArray(175* uA * cmeter ** -2* randn(200000,20), dt=defaultclock.dt)
#    noise = TimedArray(175* uA * cmeter ** -2* randn(500000,20), dt=defaultclock.dt)
    
    FLee=(0.05*mS/cm**2)/(0.4*uS/cm**2)*0.5
#    theta_phase='bad' #'good' or 'bad' or 'mixed'
    theta_phase='mixed' #'good' or 'bad' or 'mixed'


    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
#    taurinp=2*ms
#    taudinp=10*ms
#    tauinp=taudinp   
    Vhigh=0*mV
    Vlow=-80*mV
    ginp=0* msiemens * cm **-2
    
    
    N_SI,N_RS_gran,N_SI_gran=20,20,20
    all_neurons,all_synapses,all_monitors=generate_deepSI_and_gran_layers(theta_phase,N_SI,N_RS_gran,N_SI_gran,runtime)    
    
    net=Network()
    net.add(all_neurons)
    net.add(all_synapses)
    net.add(all_monitors)
    
#    taurinp=2*ms
#    taudinp=10*ms
#    tauinp=taudinp
    
#    taurinp2=0.1*ms
#    taudinp2=0.5*ms
    taurinp2=2*ms
    taudinp2=10*ms
    tauinp2=taudinp2   
    
    taurinp3=2*ms
    taudinp3=40*ms
#    taurinp3=2*ms
#    taudinp3=42.5*ms
#    taudinp3=10*ms
    tauinp3=taudinp3   
    
    noise_good=0* uA * cmeter ** -2
    noise_level=-30* uA * cmeter ** -2
    if theta_phase=='mixed':
        t0,t1=125*ms,250*ms
        i0,i1=int(t0//defaultclock.dt)+1,int(t1//defaultclock.dt)+1
        noise_array=ones((200000,20))* noise_good
#        noise_array=ones((500000,20))* noise_good
        noise_array[i0:i1,:]=noise_level* rand(12500,20)
        while t0+250*ms<runtime:
            t0,t1=t0+250*ms,t1+250*ms
            i0,i1=int(t0//defaultclock.dt)+1,int(t1//defaultclock.dt)+1
            noise_array[i0:i1,:]=noise_level* rand(12500,20)
#    print(noise_array)
    noise=TimedArray(noise_array,dt=defaultclock.dt)
    
    prefs.codegen.target = 'cython' #cython=faster, numpy = default python
    
    net.run(runtime,report='text',report_period=300*second)

#    R5,R6,R7,V_RS,V_FS,V_SI=all_monitors
    R5,R6,R7,V_RS,V_FS,V_SI,inpmon,I_RS,I_FS,I_SI=all_monitors
    
#    figure()
#    plot(R7.t,R7.i+0,'b.',label='deep SI cells')
#    plot(R5.t,R5.i+20,'r.',label='gran RS')
#    plot(R6.t,R6.i+40,'k.',label='gran SI')
#    xlim(0,runtime/second)
#    legend(loc='upper left')
    
#    figure()
##    plot(inpmon.t,inpmon.Iinp2[0])
#    plot(inpmon.t,inpmon.J[0])
#    plot(inpmon.t,inpmon.J[1])
#    xlabel('Time (s)')
#    ylabel('Iinp2')
#    tight_layout()
    
    figure()
    plot(R7.t,R7.i+0,'k.',label='input')
    plot(R5.t,R5.i+20,'r.',label='RS')
    plot(R6.t,R6.i+40,'g.',label='SOM')
    xlim(0,runtime/second)
#    legend(loc='upper left')
    xlabel('Time (s)')
    ylabel('Neuron index')
    ylim(-1,61)
    
    figure()
    plot(R5.t,R5.i+0,'r.',label='RS')
    plot(R6.t,R6.i+20,'g.',label='SOM')
    xlim(0,runtime/second)
    xlabel('Time (s)',fontsize=12)
    ylabel('Neuron index',fontsize=12)
    xticks(fontsize=12)
    yticks(fontsize=12)
    legend(loc='upper left',fontsize=12)
    ylim(-1,41)
    
#    figure()
#    plot(V_RS.t,V_RS.V[0])
    
    min_t=int(50*ms*100000*Hz)
    LFP_V_RS=1/20*sum(V_RS.V,axis=0)#[min_t:]
    savetxt('FEFvm_LFP_V_RS',LFP_V_RS,delimiter=",")
    LFP_I_RS=1/20*sum(I_RS.Isyn,axis=0)#[min_t:]
    savetxt('FEFvm_LFP_I_RS',LFP_I_RS,delimiter=",")
    LFP_V_FS=1/20*sum(V_FS.V,axis=0)#[min_t:]
    
    f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_FS=signal.periodogram(LFP_V_FS, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_I_RS=signal.periodogram(LFP_I_RS, 100000,'flattop', scaling='spectrum')
    
    figure()
    subplot(221)
    plot((V_RS.t/second),LFP_V_RS)#[min_t:],LFP_V_RS)
    ylabel('LFP')
    title('gran RS cell')
    subplot(223)
    plot((V_FS.t/second),LFP_I_RS)#[min_t:],LFP_I_RS)#LFP_V_FS)
    ylabel('LFP')
    title('gran FS cell')
    
    subplot(222)
    plot(f,Spectrum_LFP_V_RS)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,100)
    title('gran RS cell')
    subplot(224)
    plot(f,Spectrum_LFP_I_RS)#LFP_V_FS)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,100)
    title('gran FS cell')
    
    figure()
    plot(f,Spectrum_LFP_V_RS)
    ylabel('Spectrum')
    xlabel('Frequency (Hz)')
    xlim(0,50)
    
    f, t, Sxx = signal.spectrogram(LFP_I_RS, 100000*Hz,nperseg=20000,noverlap=15000)
    figure()
    pcolormesh(t, f, Sxx)#, shading='gouraud')
    ylabel('Frequency [Hz]')
    xlabel('Time [sec]')
    ylim(0,50)

    f, t, Sxx = signal.spectrogram(LFP_I_RS, 100000*Hz,nperseg=20000,noverlap=15000)
    figure()
    pcolormesh(t, f, Sxx)#, cmap=)
    colorbar(format='%.1e')
    ylabel('Frequency (Hz)')
    xlabel('Time (s)')
    ylim(0,45)
    title('Power ($V^2$)')
    
    #get spectrum from RS spikes
    N_RS_spikes=array(zeros_ones_monitor(R5,defaultclock.dt,runtime))
    savetxt('FEFvm_N_RS_spikes',N_RS_spikes,delimiter=",")
#    f, t, Sxx = signal.spectrogram(array(N_RS_spikes), 100000*Hz,nperseg=20000,noverlap=15000)
#    f, t, Sxx = signal.spectrogram(array(N_RS_spikes), 100000*Hz,nperseg=12500,noverlap=10000)
#    f, t, Sxx = signal.spectrogram(array(N_RS_spikes), 100000*Hz,nperseg=25000,noverlap=24000)
#    f, t, Sxx = signal.spectrogram(array(N_RS_spikes), 100000*Hz,nperseg=50000,noverlap=49500)
    f, t, Sxx = signal.spectrogram(LFP_I_RS, 100000*Hz,nperseg=50000,noverlap=49500)
    figure()
    pcolormesh(t, f, Sxx)#, cmap=)
    colorbar(format='%.1e')
    ylabel('Frequency (Hz)',fontsize=12)
    xlabel('Time (s)',fontsize=12)
    xticks(fontsize=12)
    yticks(fontsize=12)
    ylim(0,45)
    title('Power ($V^2$)',fontsize=12)
    
    spike_shape = atleast_2d(shape(LFP_I_RS))
    print("Size of LFP_I_RS = "+str(spike_shape[0]))
    
    
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
    
# #    ind_tround_min=where(t>=0.125)[0][0]
# #    ind_tround_max=where(t<=1.875)[0][-1]
# #    ntpertheta=(ind_tround_max+1-ind_tround_min)//7
#     ind_tround_min=where(t>=0.25)[0][0]
#     ind_tround_max=where(t<=1.75)[0][-1]
#     ntpertheta=(ind_tround_max+1-ind_tround_min)//6
# #    Sxxfolded = Sxx[:,ind_tround_min:ind_tround_max+1].reshape((len(f), ntpertheta, 7), order='F')
#     Sxxfolded = Sxx[:,ind_tround_min:ind_tround_max+1].reshape((len(f), ntpertheta, 6), order='F')
#     SxxmeanTheta = nanmean(Sxxfolded, axis = 2)
# #    figure()
# #    pcolormesh(t[:ntpertheta], f, SxxmeanTheta)#, cmap=) 
# ##    pcolormesh(t[:ntpertheta], f, Sxxfolded[:,:,0])#, cmap=) 
# #    colorbar(format='%.1e')
# #    ylabel('Frequency (Hz)',fontsize=12)
# #    xlabel('Time (s)',fontsize=12)
# #    xticks(fontsize=12)
# #    yticks(fontsize=12)
# #    ylim(0,45)
# #    title('Power ($V^2$)',fontsize=12)
    
#     figure()
#     pcolormesh(t[:ntpertheta]-t[ind_tround_min], f, hstack((SxxmeanTheta[:,ntpertheta//2:],SxxmeanTheta[:,:ntpertheta//2])))#, cmap=) 
# #    pcolormesh(t[:ntpertheta], f, Sxxfolded[:,:,0])#, cmap=) 
#     colorbar(format='%.1e')
#     ylabel('Frequency (Hz)',fontsize=12)
#     xlabel('Time (s)',fontsize=12)
#     xticks(fontsize=12)
#     yticks(fontsize=12)
#     ylim(0,45)
#     title('Power ($V^2$)',fontsize=12)
    
#     SxxpctTheta = pctMean(hstack((SxxmeanTheta[:,ntpertheta//2:],SxxmeanTheta[:,:ntpertheta//2])), 1)
#     figure()
#     pcolormesh(t[:ntpertheta]-t[ind_tround_min], f, SxxpctTheta)
#     colorbar(format='%.1e')
#     ylabel('Frequency (Hz)',fontsize=12)
#     xlabel('Time (s)',fontsize=12)
#     xticks(fontsize=12)
#     yticks(fontsize=12)
#     ylim(0,45)    
    
#     record_dt=1/100000*second#1/512*second
#     t=int(0.3*second/record_dt) #t_debut
#     L=int(2*second/record_dt)
#     fs = 1/record_dt
#     #LFP_V_RS=1/20*sum(V_RS.V,axis=0)
#     LFP_V_RS=array([0]+N_RS_spikes)
    
    end_length = 5000
    #LFPflip = flipEnds(N_RS_spikes[:,None], end_length)
    LFPflip = flipEnds(LFP_I_RS[:, None], end_length)#transpose(atleast_2d(LFP_V_RS)), end_length)
    #LFPflip = flipEnds(LFP_V_RS[:, None], end_length)#transpose(atleast_2d(LFP_V_RS)), end_length)

    
#    freq = linspace(1/second, 100*Hz, 100)
#    widths = 6*fs/freq#linspace(3, 30, 100)*fs/(2*freq*pi)
    fs=100000*Hz
    freq = linspace(1/second, 50*Hz, 50)
    #widths = around(linspace(3,5,50)*fs/freq) 
    widths = fs/freq #linspace(3, 30, 100)*fs/(2*freq*pi)
    
    CWT = signal.cwt(squeeze(LFPflip), signal.morlet2, widths)#, w=6)
    CWT = CWT[:, end_length:-end_length]
    figure()
    pcolormesh(V_RS.t, freq, absolute(CWT))
    ylabel('Frequency [Hz]',fontsize=12)
    xlabel('Time [sec]',fontsize=12)
    xticks(fontsize=12)
    yticks(fontsize=12)
    ylim(0, 50)
    colorbar()

    
    if theta_phase=='mixed':
        CWTfolded = CWT.reshape((len(freq), 25000, 8), order='F')
        CWTmeanTheta = nanmean(CWTfolded, axis = 2)
        CWTpctTheta = pctMean(absolute(CWTmeanTheta), 1)
        
#        
        figure()
        #f, t, Sxx = signal.spectrogram(LFP_LIP, 100000*Hz,nperseg=30000,noverlap=25000)
        pcolormesh(V_RS.t[0:25000], freq, absolute(CWTmeanTheta))#, cmap='RdBu')#, shading='gouraud')
        ylabel('Frequency [Hz]',fontsize=12)
        xlabel('Time [sec]',fontsize=12)
        xticks(fontsize=12)
        yticks(fontsize=12)
        ylim(0, 50)
        colorbar()

#        
        figure()
        #f, t, Sxx = signal.spectrogram(LFP_LIP, 100000*Hz,nperseg=30000,noverlap=25000)
        pcolormesh(V_RS.t[0:25000], freq, CWTpctTheta)#, cmap='RdBu')#, shading='gouraud')
        ylabel('Frequency [Hz]',fontsize=12)
        xlabel('Time [sec]',fontsize=12)
        xticks(fontsize=12)
        yticks(fontsize=12)
        ylim(0, 50)
        colorbar()
#        savetxt('data.csv', CWTpctTheta, delimiter=',')

    
    #clear_cache('cython')