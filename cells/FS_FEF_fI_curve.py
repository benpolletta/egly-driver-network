# -*- coding: utf-8 -*-


from brian2 import *

defaultclock.dt = 0.01*ms

eq_FS_FEF='''
dV/dt=1/C_FS*(-J-IL-INa-IK) : volt
J : amp * meter ** -2

IL=gL_FS*(V-VL_FS) : amp * meter ** -2
INa=gNa_FS*m0**3*h*(V-VNa_FS) : amp * meter ** -2
    m0=1/(1+exp((-V-38*mV)/10/mV)) : 1
    dh/dt=1/tauh*(hinf-h) : 1
    hinf=1/(1+exp((V+58.3*mV)/6.7/mV)) : 1
    tauh=0.225*ms+1.125*ms/(1+exp((V+37*mV)/15/mV)) : second
IK=gK_FS*m**4*(V-VK_FS) : amp * meter ** -2
    dm/dt=1/taum*(minf-m) : 1
    minf=1/(1+exp((-V-27*mV)/11.5/mV)) : 1
    taum=0.25*ms+4.35*ms*exp(-abs(V+10*mV)/10/mV) : second
    
'''


##Constants :
C_FS = 0.9* ufarad * cm ** -2

gL_FS=1 * msiemens * cm **-2
VL_FS=-65*mV

gNa_FS=200 * msiemens * cm **-2
VNa_FS=50*mV

gK_FS=20 * msiemens * cm **-2
VK_FS=-100*mV


if __name__=='__main__' :
    start_scope()
    
    FS=NeuronGroup(100,eq_FS_FEF,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    FS.V = '-70*mvolt+10*rand()*mvolt'
    FS.h = '0+0.05*rand()'
    FS.m = '0+0.05*rand()'
#    FS.J='5 * uA * cmeter ** -2'
    FS.J='5* uA * cmeter ** -2-i/100 *6 * uA * cmeter ** -2'
    
    V1=StateMonitor(FS,'V',record=[0])
    
    R1=SpikeMonitor(FS,record=True)
    
#    I1=StateMonitor(FS,'IL',record=[0])
#    I2=StateMonitor(FS,'INa',record=[0])
#    I3=StateMonitor(FS,'IK',record=[0])
    
    run(10*second)
    
    figure()
    plot(V1.t/second,V1.V[0]/volt)
    xlabel('Time (s)')
    ylabel('Membrane potential (V)')
    title('FS cell')
    
    figure()
    plot((-FS.J/ (uA * cmeter ** -2)),R1.count/10)
    xlabel(r'I ($\mu A/cm^2$)', fontsize=14)
    ylabel('f (Hz)', fontsize=14)
    tick_params(axis='both', which='major', labelsize=12)
    
#    figure()
#    plot(I1.t/second,I1.IL[0],label='L')
#    plot(I1.t/second,I2.INa[0],label='Na')
#    plot(I1.t/second,I3.IK[0],label='K')
#    plot(I1.t/second,I4.IAR[0],label='AR')
#    title('Synaptic currents')
#    legend()