# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 14:28:44 2019

@author: aaussel
"""

from brian2 import *

from scipy import signal
from cells.RS_FEF import *
from cells.FS_FEF import *
from cells.SI_FEF import *
from cells.VIP_FEF import *

def generate_deepSI_and_gran_layers(theta_phase,N_VIP,N_RS_gran,N_FS_gran,N_SI_gran,runtime):
    
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
    
    ##Define neuron groups
    E_gran=NeuronGroup(N_RS_gran,eq_RS_FEF,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    E_gran.V = '-70*mvolt+10*rand()*mvolt'
    E_gran.h = '0+0.05*rand()'
    E_gran.m = '0+0.05*rand()'
    E_gran.mAR = '0.035+0.025*rand()'
#    E_gran.J='30 * uA * cmeter ** -2'  #article SI=25, code=1
#    E_gran.J='20 * uA * cmeter ** -2'  #article SI=25, code=1
    E_gran.J='10 * uA * cmeter ** -2'  #article SI=25, code=1
    #0   
#    E_gran.J='20 * uA * cmeter ** -2'
    
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
    SI_gran.J='5 * uA * cmeter ** -2' #article=code=35
     #-30
#    SI_gran.J='-25 * uA * cmeter ** -2'
     
#    SI_deep=NeuronGroup(N_SI,eq_SIdeep,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
#    SI_deep.V = '-100*mvolt+10*rand()*mvolt'
#    SI_deep.h = '0+0.05*rand()'
#    SI_deep.m = '0+0.05*rand()'
#    SI_deep.mAR = '0.02+0.04*rand()'
#    SI_deep.J='35* uA * cmeter ** -2' #article SI=50, code=35, Mark = 45

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
    #S_EgranEgran=generate_syn(E_gran,E_gran,'IsynEgran','',0.4*usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)
    #S_EgranEgran=generate_syn(E_gran,E_gran,'IsynEgran','',1/160*msiemens * cm **-2,0.125*ms,1*ms,0*mV)
    S_EgranEgran=generate_syn(E_gran,E_gran,'IsynRS_FEF_VM','',0.6*msiemens * cm **-2,0.125*ms,1*ms,0*mV) #0.4 #0.6
    S_EgranFSgran=generate_syn(E_gran,FS_gran,'IsynRS_FEF_VM','',0.2*msiemens * cm **-2,0.125*ms,1*ms,0*mV) #0.2
    #S_EgranSIgran=generate_syn(E_gran,SI_gran,'IsynEgran','',0.2*usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)
    S_EgranSIgran=generate_syn(E_gran,SI_gran,'IsynRS_FEF_VM','',0.5*msiemens * cm **-2,0.125*ms,1*ms,0*mV) #0.6 #0.5
    #S_EgranRS=generate_syn(E_gran,RS,'IsynEgran','',0.2*usiemens * cm **-2*FLee,0.125*ms,1*ms,0*mV)

    #From FS cells
    S_FSgranFSgran=generate_syn(FS_gran,FS_gran,'IsynFS_FEF_VM','',0.2*msiemens * cm **-2,0.25*ms,5*ms,-80*mV) #0.2
    S_FSgranEgran=generate_syn(FS_gran,E_gran,'IsynFS_FEF_VM','',0.2*msiemens * cm **-2,0.25*ms,5*ms,-80*mV) #0.2
    S_FSgranSIgran=generate_syn(FS_gran,SI_gran,'IsynFS_FEF_VM','',1*msiemens * cm **-2,0.25*ms,5*ms,-80*mV) #0.5

    #From SOM cells, normal timescale
    #S_FSgranEgran=generate_syn(SI_gran,E_gran,'IsynFSgran','',1* usiemens * cm **-2*FLee,0.25*ms,5*ms,-80*mV)
    #S_FSgranEgran=generate_syn(SI_gran,E_gran,'IsynSI_FEF_VM','',0.6* msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
    S_SIgranEgran=generate_syn(SI_gran,E_gran,'IsynSI_FEF_VM','',0.5*msiemens * cm **-2,0.25*ms,20*ms,-80*mV) #0.35 #0.5
    S_SIgranFSgran=generate_syn(SI_gran,FS_gran,'IsynSI_FEF_VM','',0.2*msiemens * cm **-2,0.25*ms,20*ms,-80*mV) #0.2
    #S_FSgranFSgran=generate_syn(SI_gran,SI_gran,'IsynFSgran','',0.1* usiemens * cm **-2*FLee,0.25*ms,5*ms,-75*mV)
    S_SIgranSIgran=generate_syn(SI_gran,SI_gran,'IsynSI_FEF_VM','i',1* msiemens * cm **-2,0.25*ms,20*ms,-75*mV) #1 #0.2
    
#    #From SOM cells, gamma timescale
#    S_FSgranEgran=generate_syn(SI_gran,E_gran,'IsynSI_FEF_VM','',0.5*msiemens * cm **-2,0.25*ms,5*ms,-80*mV) #0.35
#    S_FSgranFSgran=generate_syn(SI_gran,SI_gran,'IsynSI_FEF_VM','',0.2* msiemens * cm **-2,0.25*ms,5*ms,-75*mV) #1
    
    
    #From VIP cells    
    #S_SIdeepFSgran=generate_syn(VIP,SI_gran,'IsynSIdeep','',0.4* usiemens * cm **-2*FLee,0.25*ms,20*ms,-80*mV)
#    S_SIdeepFSgran=generate_syn(VIP,SI_gran,'IsynSI2_FEF_VM','',1*msiemens * cm **-2,0.25*ms,20*ms,-80*mV)
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
    
    #Theta=4Hz
    if theta_phase=='mixed':
        t0=0*ms
        t1=125*ms
        inputs_topdown3=generate_spike_timing(N_VIP,fIB,t0,end_time=t1)
        while t0+250*ms<runtime:
            t0,t1=t0+250*ms,t1+250*ms
            inputs_topdown3=vstack((inputs_topdown3,generate_spike_timing(N_VIP,fIB,t0,end_time=t1))) # t0+1/(4*fIB),end_time=t1)))
            
#    #Theta=8Hz
#    if theta_phase=='mixed':
#        t0=0*ms
#        t1=62.5*ms
#        inputs_topdown3=generate_spike_timing(N_VIP,fIB,t0,end_time=t1)
#        while t0+125*ms<runtime:
#            t0,t1=t0+125*ms,t1+125*ms
#            inputs_topdown3=vstack((inputs_topdown3,generate_spike_timing(N_VIP,fIB,t0,end_time=t1)))

#    #Theta=2Hz
#    if theta_phase=='mixed':
#        t0=0*ms
#        t1=250*ms
#        inputs_topdown3=generate_spike_timing(N_VIP,fIB,t0,end_time=t1)
#        while t0+500*ms<runtime:
#            t0,t1=t0+500*ms,t1+500*ms
#            inputs_topdown3=vstack((inputs_topdown3,generate_spike_timing(N_VIP,fIB,t0,end_time=t1)))
                        
    
    G_topdown3 = SpikeGeneratorGroup(N_VIP, inputs_topdown3[:,1], inputs_topdown3[:,0]*second)
    topdown_in3=Synapses(G_topdown3,VIP,on_pre='Vinp=Vhigh')
    topdown_in3.connect(j='i')
        
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
    #    inputs_lateral2=generate_spike_timing(N_VIP,25*Hz,0*ms,end_time=2100*ms)
    #    G_lateral2 = SpikeGeneratorGroup(N_VIP, inputs_lateral2[:,1], inputs_lateral2[:,0]*second)
    #    lateral_in2=Synapses(G_lateral2,SI,on_pre='Vinp=Vhigh')
    #    lateral_in2.connect(j='i')
        
    if input_thalamus_gran:    
        SI_gran.ginp_SI=thal_cond
        E_gran.ginp_RS=thal_cond
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
#            fLIP=12*Hz #test, if LIP hasn't switched to its good phase activity
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
                bottomup=vstack((bottomup,generate_spike_timing(N_SI_gran,fLIP,t0,end_time=t1))) # t0+1/(4*fLIP),end_time=t1)))
            
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
#         Poisson_input3 = SpikeGeneratorGroup(N_FS_gran, bottomup[:,1], bottomup[:,0]*second)
#         bottomup_in3=Synapses(Poisson_input3,FS_gran,on_pre='Vinp=Vhigh')
#         bottomup_in3.connect(j='i')
    
    #Define monitors and run network :
    R5=SpikeMonitor(E_gran,record=True)
    R6=SpikeMonitor(SI_gran,record=True)
    R7=SpikeMonitor(VIP,record=True)
    R8=SpikeMonitor(FS_gran,record=True)
    BUinp=SpikeMonitor(Poisson_input,record=True)
    TDinp=SpikeMonitor(G_topdown3,record=True)
    
    inpmon=StateMonitor(E_gran,'Iinp1',record=True)
    #graninpmon=StateMonitor(FS,'IsynEgran',record=[0])
    #inpIBmon=StateMonitor(IB_bd,'Iapp',record=[0])
    
    V_RS=StateMonitor(E_gran,'V',record=True)
    V_SI=StateMonitor(SI_gran,'V',record=True)
    V_VIP=StateMonitor(VIP,'V',record=True)
    
    all_neurons=VIP,E_gran,FS_gran,SI_gran,G_topdown3,Poisson_input,Poisson_input2#,Poisson_input3
    all_synapses=S_EgranEgran,S_EgranFSgran,S_EgranSIgran,S_FSgranEgran,S_FSgranFSgran,S_FSgranSIgran,S_SIgranEgran,S_SIgranFSgran,S_SIgranSIgran,S_VIPSIgran,gap_SISI,topdown_in3,bottomup_in,bottomup_in2#,bottomup_in3
    all_monitors=R5,R6,R7,R8,V_RS,V_SI,V_VIP,BUinp,TDinp
#    all_monitors=R5,R6,R7,V_RS,V_FS,V_SI,inpmon
    
    return all_neurons,all_synapses,all_monitors



if __name__=='__main__':
    close('all')
    start_scope()    
    
    prefs.codegen.target = 'numpy'
    defaultclock.dt = 0.01*ms
    
    FLee=(0.05*mS/cm**2)/(0.4*uS/cm**2)*0.5
    theta_phase='mixed' #'good' or 'bad' or 'mixed'
    runtime=2*second
    
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
    
    N_VIP,N_RS_gran,N_FS_gran,N_SI_gran=20,20,20,20
    all_neurons,all_synapses,all_monitors=generate_deepSI_and_gran_layers(theta_phase,N_VIP,N_RS_gran,N_FS_gran,N_SI_gran,runtime)    
    
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
#    taurinp2=10*ms
#    taudinp2=50*ms
    tauinp2=taudinp2   
    
    prefs.codegen.target = 'cython' #cython=faster, numpy = default python
    
    net.run(runtime,report='text',report_period=300*second)

    R5,R6,R7,R8,V_RS,V_SI,V_VIP,BUinp,TDinp=all_monitors
#    R5,R6,R7,V_RS,V_FS,V_SI,inpmon=all_monitors
    
#    figure()
#    plot(R7.t,R7.i+0,'b.',label='deep SI cells')
#    plot(R5.t,R5.i+20,'r.',label='gran RS')
#    plot(R6.t,R6.i+40,'k.',label='gran SI')
#    xlim(0,runtime/second)
#    legend(loc='upper left')
    
#    figure()
#    plot(inpmon.t,inpmon.Iinp1[0])
#    xlabel('Time (s)')
#    ylabel('Iinp1')
#    tight_layout()
    
    figure()
    plot(BUinp.t,BUinp.i+0,'k.',label='BU')
    plot(TDinp.t,TDinp.i+20,'r.',label='TD')
    xlim(0,runtime/second)
#    legend(loc='upper left')
    xlabel('Time (s)')
    ylim(-1,41)
    
    figure()
    plot(R7.t,R7.i+0,'k.',label='VIP')
    plot(R5.t,R5.i+20,'r.',label='RS')
    plot(R6.t,R6.i+40,'g.',label='SOM')
    plot(R8.t,R8.i+60,'b.',label='FS')
    xlim(0,runtime/second)
#    legend(loc='upper left')
    xlabel('Time (s)')
#    ylabel('Neuron index')
    ylim(-1,81)
    
#    figure()
#    plot(V_RS.t,V_RS.V[0])
    
    min_t=int(50*ms*100000*Hz)
    LFP_V_RS=1/20*sum(V_RS.V,axis=0)[min_t:]
    LFP_V_SI=1/20*sum(V_SI.V,axis=0)[min_t:]
    
    f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='spectrum')
    f,Spectrum_LFP_V_SI=signal.periodogram(LFP_V_SI, 100000,'flattop', scaling='spectrum')
    
    figure()
    subplot(221)
    plot((V_RS.t/second)[min_t:],LFP_V_RS)
    ylabel('LFP')
    title('gran RS cell')
    subplot(223)
    plot((V_SI.t/second)[min_t:],LFP_V_SI)
    ylabel('LFP')
    title('gran SOM cell')
    
    subplot(222)
    plot(f,Spectrum_LFP_V_RS)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,100)
    title('gran RS cell')
    subplot(224)
    plot(f,Spectrum_LFP_V_SI)
    ylabel('Spectrum')
    yticks([],[])
    xlim(0,100)
    title('gran SOM cell')
    
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
    
    clear_cache('cython')