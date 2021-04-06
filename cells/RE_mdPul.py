# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 15:19:08 2019

@author: Amelie Aussel
"""
from brian2 import *
defaultclock.dt = 0.01*ms

prefs.codegen.target = 'numpy'

eq_RE_mdPul='''
    dV/dt=(-INa-IK-IL-IKL-ITRE-Isyn-J-Iapp-Itheta)/Cm_RE : volt
    
    J : amp * meter ** -2
    
    Isyn=IsynHTC+IsynTC+IsynI+IsynRE+IsynREA+IsynREB+Isyn_FEF+Isyn_LIP : amp * meter ** -2
    IsynHTC : amp * meter ** -2
    IsynTC : amp * meter ** -2
    IsynI : amp * meter ** -2
    IsynRE : amp * meter ** -2
    IsynREA : amp * meter ** -2
    IsynREB : amp * meter ** -2
    Isyn_FEF : amp * meter ** -2
    Isyn_LIP : amp * meter ** -2
    
    Itheta = amptheta * int(sin(2*pi*time_total*ftheta+offsettheta)<0) : amp * meter ** -2
    amptheta : amp * meter ** -2
    ftheta : Hz    
    dtime_total/dt=1 : second
    offsettheta : 1
    
    Vt=V+55*mV : volt
    IK=gK_RE*mK**4*(V-EK_RE) : amp * meter ** -2
        dmK/dt=(mKinf-mK)/taumK : 1
        mKinf = alphamK/(alphamK+betamK) : 1
        taumK=1/(alphamK+betamK) : second
        alphamK = (0.032*(15*mV-Vt)/mV)/(exp((15*mV-Vt)/5/mV)-1)/ms : hertz
        betamK=0.5 * exp((10*mV-Vt)/40/mV)/ms : hertz
        
    INa=gNa_RE*mNa**3*hNa*(V-ENa_RE) : amp * meter ** -2
        dmNa/dt = (mNainf-mNa)/taumNa : 1
        dhNa/dt = (hNainf-hNa)/tauhNa : 1
        mNainf = alphamNa/(alphamNa+betamNa) : 1
        taumNa=1*ms/(alphamNa+betamNa) : second      
        hNainf = alphahNa/(alphahNa+betahNa) : 1
        tauhNa=1*ms/(alphahNa+betahNa) : second
        alphamNa=(0.32*(13*mV-Vt))/(exp((13*mV-Vt)/4/mV)-1)/mV : 1
        betamNa = (0.28*(Vt-40*mV))/(exp((Vt-40*mV)/5/mV)-1)/mV : 1
        alphahNa=(0.128*exp((17*mV-Vt)/18/mV)) : 1
        betahNa = 4/(1+exp((40*mV-Vt)/5/mV)) : 1    
        
    IL=gL_RE*(V-EL_RE) : amp * meter ** -2
    IKL=gKL_RE*(V-EKL_RE) : amp * meter ** -2
    
    ITRE = gTRE * mTRE**2 * hTRE * (V-ETRE) : amp * meter ** -2
        dmTRE/dt = (mTREinf-mTRE)/taumTRE : 1
        dhTRE/dt = (hTREinf-hTRE)/tauhTRE : 1 
        mTREinf = 1/(1+exp(-(V+52*mV)/7.4/mV)) : 1
        taumTRE= .999 * ms + .333 / (exp((V+27*mV)/2/mV)+exp(-(V+102*mV)/15/mV)) * ms : second      
        hTREinf = 1/(1+exp((V+80*mV)/5/mV)) : 1
        tauhTRE= 0.1*(28.307 * ms + .333 / (exp((V+48*mV)/4/mV)+exp(-(V+407*mV)/50/mV)) * ms) : second
        dCa_i/dt = (-10*ITRE)/(2*96485.3* coulomb * mole ** -1 * meter)*int((-10*ITRE)/(2*96485.3* coulomb * mole ** -1 * meter)>0*mmolar/ms)- (Ca_i-0.00024*mmolar)/5/ms : mole * meter**-3
        ETRE=R*T/(z*F)*log(2*mmolar/Ca_i)/log(2) : volt

    Iapp=sinp*ginp_RE*(V-Vrev_inp) : amp * meter ** -2
        dsinp/dt=-sinp/taudinp + (1-sinp)/taurinp*0.5*(1+tanh(Vinp/10/mV)) : 1
        dVinp/dt=1/tauinp*(Vlow-Vinp) : volt
        ginp_RE : siemens * metre **-2    
'''


##Constants :
Cm_RE = 2.5* ufarad * cm ** -2
gK_RE=10e-3 * siemens * cm **-2
EK_RE=-100*mV
gNa_RE=100e-3 * siemens * cm **-2
ENa_RE=50*mV
gL_RE=0.01e-3 * siemens * cm **-2   ### correct value ?
EL_RE=-73*mV
gKL_RE=0.0028e-3 * siemens * cm **-2   ### correct value ?
EKL_RE=-100*mV
gTRE= 2.3e-3 * siemens * cm **-2
R = 8.314 * joule * kelvin**-1 * mole **-1
T = (273.15+37.2)*kelvin
z = 2
F = 96485.3* coulomb * mole ** -1


if __name__=='__main__' :
    start_scope()
    RE=NeuronGroup(1,eq_RE_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    RE.V = '-60*mvolt'
    RE.J = '-100 * nA * cmeter ** -2'
    RE.Ca_i = '1e-7 * mole * metre**-3'
    
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp_IB=0* msiemens * cm **-2
    ginp=0* msiemens * cm **-2
    
    V1=StateMonitor(RE,'V',record=[0])
    V2=StateMonitor(RE,'Ca_i',record=[0])
    V3=StateMonitor(RE,'Iapp',record=[0])
#    V5=StateMonitor(RE,'last_t_EPSP',record=[0])
    #R1=SpikeMonitor(SI)
    
    run(1*second)
    
    figure()
    plot(V1.t/second,V1.V[0]/volt)
    xlabel('Time (s)')
    ylabel('Membrane potential (V)')
    title('RE cell')
    
#    figure()
#    plot(V2.t/second,V2.Ca_i[0])
#    xlabel('Time (s)')
#    ylabel('Calcium concentration')