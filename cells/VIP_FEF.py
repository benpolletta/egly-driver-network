#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 13:17:48 2020

@author: amelie
"""

from brian2 import *

defaultclock.dt = 0.01*ms

eq_VIP='''
dV/dt = (1/Cm)*(-INa-IK-ID-IL-Isyn+Irand+Iapp-Iapp2) : volt
    
INa=gna*h*(minf*minf*minf)*(V-Ena) : amp * meter ** -2
    minf = 1/(1+exp(-(V+24*mV)/11.5/mV)) : 1
    dh/dt= (hinf - h)/tauh : 1
    hinf = 1/(1+exp((V+58.3*mV)/6.7/mV)) : 1
    tauh =  0.5*msecond + 14*msecond/(1+exp((V+60*mV)/12/mV)) : second
    
IK=gk*(n*n)*(V-Ek) : amp * meter ** -2
    dn/dt = (ninf - n)/taun : 1
    ninf = 1/(1+ exp(-(V+12.4*mV)/6.8/mV)) : 1
    taun = 1*msecond* (0.087 + 11.4/(1+exp((V+14.6*mV)/8.6/mV))) * (0.087 + 11.4/(1+exp(-(V-1.3*mV)/18.7/mV))) : second

ID=gd*a*a*a*b*(V-Ek) : amp * meter ** -2
    da/dt=(ainfD - a)/2/msecond : 1
    ainfD =  1/(1 + exp(-(V+50*mV)/20/mV)) : 1
    db/dt=(binfD - b)/150/msecond : 1
    binfD =  1/(1 + exp((V+70*mV)/6/mV)) : 1

IL=gl*(V-El) : amp * meter ** -2
Irand=0*4*sqrt(0.05)*rand()*mamp * cmeter ** -2 : amp * meter ** -2 (constant over dt)
Iapp : amp * meter ** -2 
Isyn=IsynRS_FEF_VM+IsynSI_FEF_VM+IsynSI2_FEF_VM+IsynRS_FEF_V+IsynFS_FEF_V+Isyn_mdPul : amp * meter ** -2
Iapp2=sinp*ginpVIP*(V-Vrev_inp) : amp * meter ** -2
    dsinp/dt=-sinp/taudinp2 + (1-sinp)/taurinp2*0.5*(1+tanh(Vinp/10/mV)) : 1
    dVinp/dt=1/tauinp2*(Vlow-Vinp) : volt
    ginpVIP = ginp_VIP_good* int(sin(2*pi*t*4*Hz)>0) + ginp_VIP_bad* int(sin(2*pi*t*4*Hz)<=0) : siemens * meter **-2
    ginp_VIP_good : siemens * meter **-2
    ginp_VIP_bad : siemens * meter **-2
IsynRS_FEF_VM : amp * meter ** -2
IsynSI_FEF_VM : amp * meter ** -2
IsynSI2_FEF_VM : amp * meter ** -2
IsynRS_FEF_V : amp * meter ** -2
IsynFS_FEF_V : amp * meter ** -2
Isyn_mdPul : amp * meter ** -2
'''


##Constants :
Cm = 2 * ufarad * cm ** -2
gna = 112.5 * msiemens * cm **-2
gk = 225 * msiemens * cm **-2
gd = 4 * msiemens * cm **-2 #4
gl = 0.25 * msiemens * cm **-2
Ena = 50 *mV
Ek  = -90 *mV
El  = -70 *mV


if __name__=='__main__' :
    start_scope()
    
    Vrev_inp=0*mV
#    taurinp=0.1*ms
#    taudinp=0.5*ms
    taurinp=1*ms
    taudinp=5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    
    VIP=NeuronGroup(1,eq_VIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    VIP.V = '-63*mvolt'
#    VIP.Iapp='5 * uA * cmeter ** -2'
    VIP.Iapp='0 * uA * cmeter ** -2'
    VIP.Vinp=Vlow
    
    G_topdown3 = SpikeGeneratorGroup(1, array([0,0]), array([0.25,0.75])*second)
    topdown_in3=Synapses(G_topdown3,VIP,on_pre='Vinp=Vhigh')
    topdown_in3.connect(j='i')
    VIP.ginp_VIP_bad=10* msiemens * cm **-2
    VIP.ginp_VIP_good=10* msiemens * cm **-2
    
    V1=StateMonitor(VIP,'V',record=[0])
    
#    I1=StateMonitor(FS,'IL',record=[0])
#    I2=StateMonitor(FS,'INa',record=[0])
#    I3=StateMonitor(FS,'IK',record=[0])
    
    run(1*second)
    
    figure()
    plot(V1.t/second,V1.V[0]/volt)
    xlabel('Time (s)')
    ylabel('Membrane potential (V)')
    title('VIP cell')
    
#    figure()
#    plot(I1.t/second,I1.IL[0],label='L')
#    plot(I1.t/second,I2.INa[0],label='Na')
#    plot(I1.t/second,I3.IK[0],label='K')
#    plot(I1.t/second,I4.IAR[0],label='AR')
#    title('Synaptic currents')
#    legend()