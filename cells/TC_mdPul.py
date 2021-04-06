# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 08:38:28 2019

@author: Amelie Aussel
"""

from brian2 import *

defaultclock.dt = 0.01*ms
prefs.codegen.target = 'numpy'

eq_TC_mdPul='''
    dV/dt=(-INa-IK-IL-IKL-IH-ITLT-J-Isyn-Iapp-Itheta)/Cm_TC : volt
    
    J : amp * meter ** -2
    
    Itheta = amptheta * int(sin(2*pi*time_total*ftheta)<0) : amp * meter ** -2
    amptheta : amp * meter ** -2
    ftheta : Hz    
    dtime_total/dt=1 : second    
       
    Isyn=IsynHTC+IsynTC+IsynI+IsynREGABAA+IsynREGABAB+IsynREA+IsynREB+Isyn_FEF+Isyn_LIP : amp * meter ** -2
    IsynHTC : amp * meter ** -2
    IsynTC : amp * meter ** -2
    IsynI : amp * meter ** -2
    IsynREGABAA : amp * meter ** -2
    IsynREGABAB : amp * meter ** -2
    IsynREA : amp * meter ** -2
    IsynREB : amp * meter ** -2
    Isyn_FEF : amp * meter ** -2
    Isyn_LIP : amp * meter ** -2
    
    Vt=V+25*mV : volt
    INa=gNa_TC*mNa**3*hNa*(V-ENa_TC)*int(mNa>0)*int(hNa>0) : amp * meter ** -2
        dmNa/dt = (mNainf-mNa)/taumNa : 1
        dhNa/dt = (hNainf-hNa)/tauhNa : 1
        mNainf = alphamNa/(alphamNa+betamNa) : 1
        taumNa=1*msecond/(alphamNa+betamNa) : second      
        hNainf = alphahNa/(alphahNa+betahNa) : 1
        tauhNa=1*msecond/(alphahNa+betahNa) : second
        alphamNa=(0.32*(13*mV-Vt))/(exp((13*mV-Vt)/4/mV)-1)/mV : 1
        betamNa = (0.28*(Vt-40*mV))/(exp((Vt-40*mV)/5/mV)-1)/mV : 1
        alphahNa=(0.128*exp((17*mV-Vt)/18/mV)) : 1
        betahNa = 4/(1+exp((40*mV-Vt)/5/mV)) : 1  
        
    IK=gK_TC*mK**4*(V-EK_TC) : amp * meter ** -2
        dmK/dt=(mKinf-mK)/taumK : 1
        mKinf = alphamK/(alphamK+betamK) : 1
        taumK=1/(alphamK+betamK) : second
        alphamK = (0.032*(15*mV-Vt)/mV)/(exp((15*mV-Vt)/5/mV)-1)/ms : hertz
        betamK=0.5 * exp((10*mV-Vt)/40/mV)/ms : hertz
        
    IL=gL_TC*(V-EL_TC) : amp * meter ** -2
    IKL=gKL_TC*(V-EKL_TC) : amp * meter ** -2
    
    ITLT = gTLT_TC * mTLT**2 * hTLT * (V-ETLT) : amp * meter ** -2
        dmTLT/dt = (mTLTinf-mTLT)/taumTLT : 1
        dhTLT/dt = (hTLTinf-hTLT)/tauhTLT : 1 
        mTLTinf = 1/(1+exp(-(V+59*mV)/6.2/mV)) : 1
        taumTLT= 0.1*ms: second      
        hTLTinf = 1/(1+exp((V+83*mV)/4/mV)) : 1
        tauhTLT= 0.1*msecond*(30.8+ (211.4 + exp((V+115.2*mV)/5/mV))/(1+exp((V+86*mV)/3.2/mV)))/3.7372 : second
        dCa_i/dt = (-10000000*ITLT)/(2*96485.3* coulomb * mole ** -1 * meter)*int(ITLT<0*amp * meter ** -2)- (Ca_i-0.00024*mmolar)/5/ms : mole * meter**-3
        ETLT=120*mV : volt
    
    IH = gH_TC * r * (V-EH_TC) : amp * meter ** -2
        dr/dt=(rinf-r)/tausr : 1
        tausr = 20*ms + 1000*ms / (exp((V+71.5*mV)/14.2/mV)+exp(-(V+89*mV)/11.6/mV)) : second
        rinf = 1/(1+exp((V+75*mV)/5.5/mV)) : 1

    Iapp=sinp*ginp_TC*(V-Vrev_inp) : amp * meter ** -2
        dsinp/dt=-sinp/taudinp + (1-sinp)/taurinp*0.5*(1+tanh(Vinp/10/mV)) : 1
        dVinp/dt=1/tauinp*(Vlow-Vinp) : volt
        ginp_TC : siemens * metre **-2
        
    Iapp2=sinp2*ginp_TC2*(V-Vrev_inp) : amp * meter ** -2
        dsinp2/dt=-sinp2/taudinp2 + (1-sinp2)/taurinp2*0.5*(1+tanh(Vinp/10/mV)) : 1
        dVinp2/dt=1/tauinp2*(Vlow2-Vinp2) : volt
        ginp_TC2 : siemens * metre **-2   
    '''

#ETLT=R*T/(z*F)*log(2*mmolar/Ca_i)/log(2) : volt
    
Cm_TC = 2.5* ufarad * cm ** -2
gNa_TC=90e-3 * siemens * cm **-2
ENa_TC=50*mV
gK_TC=10e-3 * siemens * cm **-2 #10
EK_TC=-100*mV
gL_TC=0.01e-3 * siemens * cm **-2   #0.01
#gL_TC=0.015e-3 * siemens * cm **-2 
EL_TC=-70*mV
#gKL_TC=0.006e-3 * siemens * cm **-2   #0.006 0.0028
gKL_TC=0.0028e-3 * siemens * cm **-2 
EKL_TC=-100*mV
#gTLT_TC= 10e-3 * siemens * cm **-2
gTLT_TC= 30e-3 * siemens * cm **-2
R = 8.314 * joule * kelvin**-1 * mole **-1
T = (273.15+37.2)*kelvin
z = 2
F = 96485.3* coulomb * mole ** -1
gH_TC = 0.02e-3 * siemens * cm **-2
EH_TC = -43 * mV
a=1

N_TC=1

if __name__=='__main__':
    close('all')
    start_scope()
    TC=NeuronGroup(N_TC,eq_TC_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    TC.V = '-80*mvolt'
    TC.Ca_i = '0.00024*mmolar'
#    TC.o1 = '0.5'
#    TC.c1 = '0.5'
#    TC.p0 = '0.5'
    TC.J = '0 * nA * cmeter ** -2'
    TC.mTLT = '1'
    
    
    V1=StateMonitor(TC,'V',record=[0])
#    V2=StateMonitor(TC,'IH',record=[0])
    V3=StateMonitor(TC,'Ca_i',record=[0])
    #R1=SpikeMonitor(SI)
    V6=StateMonitor(TC,'INa',record=[0])
    V7=StateMonitor(TC,'IK',record=[0])
    V8=StateMonitor(TC,'IL',record=[0])
    V9=StateMonitor(TC,'IKL',record=[0])
#    V10=StateMonitor(TC,'cSK',record=[0])
    V11=StateMonitor(TC,'ITLT',record=[0])
    V12=StateMonitor(TC,'hTLT',record=[0])
    V13=StateMonitor(TC,'mTLT',record=[0])
    V14=StateMonitor(TC,'ETLT',record=[0])
    
    V15=StateMonitor(TC,'hNa',record=[0])
    V16=StateMonitor(TC,'mNa',record=[0])
    
    
    H1=StateMonitor(TC,'IH',record=[0])
#    H2=StateMonitor(TC,'o1',record=[0])
#    H3=StateMonitor(TC,'p0',record=[0])
#    H4=StateMonitor(TC,'c1',record=[0])
#    H5=StateMonitor(TC,'tausH',record=[0])  
#    H6=StateMonitor(TC,'hinfH',record=[0])  
    
    run(1*second)
    
    figure()
    plot(V1.t/second,V1.V[0]/volt)
    xlabel('Time (s)')
    ylabel('Membrane potential (V)')
    title('TC cell')
    
#    figure()
#    plot(V10.t/second,V10.cSK[0]/volt)
    
    figure()
    plot(V3.t/second,V3.Ca_i[0])
    xlabel('Time (s)')
    title('Ca concentration')
    
    figure()
    subplot(131)
    plot(V11.t/second,V11.ITLT[0])
    xlabel('Time (s)')
    title('I_TLT')
    subplot(132)
    plot(V12.t/second,V12.hTLT[0],label='h')
    plot(V13.t/second,V13.mTLT[0],label='m')
    xlabel('Time (s)')
    legend()
    subplot(133)
    plot(V14.t/second,V14.ETLT[0])
    xlabel('Time (s)')
    title('E_TLT')
    
    figure()
    subplot(131)
    plot(V6.t/second,V6.INa[0])
    xlabel('Time (s)')
    title('I_Na')
    subplot(132)
    plot(V15.t/second,V15.hNa[0])
    xlabel('Time (s)')
    title('h_Na')
    subplot(133)
    plot(V16.t/second,V16.mNa[0])
    xlabel('Time (s)')
    title('mNa')
   
    figure()
    plot(V6.t/second,V6.INa[0],label='INa')
    plot(V7.t/second,V7.IK[0],label='IK')
    plot(V8.t/second,V8.IL[0],label='IL')
    plot(V9.t/second,V9.IKL[0],label='IKL')
    plot(V11.t/second,V11.ITLT[0],label='ITLT')
    plot(H1.t/second,H1.IH[0],label='IH')
#    plot(V10.t/second,V10.Iapp[0],label='Iapp')
#    plot(V11.t/second,V11.IappGABA[0],label='IappGABA')
#    plot(V11.t/second,V6.INa[0]+V7.IK[0]+V8.IL[0]+V9.IKL[0]+V4.ITLT[0]+V2.IH[0]+V10.Iapp[0]+V11.IappGABA[0],label='sum of all I')
    legend()
    
#    figure()
#    plot(H1.t/second,H1.IH[0],label='IH')
#    plot(H1.t/second,H2.o1[0],label='o1')
#    plot(H1.t/second,H3.p0[0],label='p0')
#    plot(H1.t/second,H4.c1[0],label='c1')
#    plot(H1.t/second,H5.tausH[0],label='tauh')
#    plot(H1.t/second,H6.hinfH[0],label='hinf')
#    legend()
#    

    
    #figure()
    #plot(V2.t/second,-V2.Iapp[0])
    #xlabel('Time (s)')
    #ylabel('Input current')


#    IH = gH_TC * (o1 + a*(1-c1-o1)) * (V-EH_TC) : amp * meter ** -2
#        do1/dt = 0.0001/ms * (1-c1-o1) - 0.001/ms * ((1-p0)/0.01) : 1
#        dp0/dt = 0.0004/ms * (1-p0) - 0.004/ms * (Ca_i/(0.0002*mmolar))**2 : 1
#        dc1/dt = betaH * o1 - alphaH * c1 : 1
#        betaH = (1-hinfH)/tausH : hertz
#        alphaH = hinfH / tausH : hertz
#        tausH = 20 * ms + 1000* ms/(exp((V+71.5*mV)/14.2/mV)+exp(-(V+89*mV)/11.6/mV)) : second
#        hinfH = 1 / (1 + exp((V+75*mV)/5.5/mV)) : 1
    
#*int(ITLT<0*amp * meter ** -2)
    
#mTLTinf = 1/(1+exp(-(V+59*mV)/6.2/mV)) : 1
#taumTLT= 1*ms: second      
#hTLTinf = 1/(1+exp((V+83*mV)/4/mV)) : 1
#tauhTLT=  (30.8+ (211.4 + exp((V+115.2*mV)/5/mV))/(1+exp((V+86*mV)/3.2/mV)))/3.7372*ms : second
       
    
#    ISK = gSK_HTC * cSK * (V-EK_TC) : amp * meter ** -2
#    dcSK/dt = (cSKinf-cSK)/taucSK : 1
#    cSKinf = 0.81 / (1 + exp(-(log(Ca_i/mmolar)+0.3)/0.46)): 1
#    taucSK = 6.1*msecond : second