# -*- coding: utf-8 -*-


from brian2 import *

defaultclock.dt = 0.01*ms

eq_IB_soma='''
dV/dt=1/C_IB_soma*(-J-Isyn-Igap-Iran-Iapp-IL-INa-IK) : volt
J : amp * meter ** -2
Isyn : amp * meter ** -2
Iran: amp * meter ** -2 
Igap = Igap_soma+Igap_axon+Igap_ad+Igap_bd : amp * meter ** -2
Igap_soma : amp * meter ** -2
Igap_axon : amp * meter ** -2
Igap_ad : amp * meter ** -2
Igap_bd : amp * meter ** -2
IL=gL_IB_soma*(V-VL_IB_soma) : amp * meter ** -2
INa=gNa_IB_soma*m0**3*h*(V-VNa_IB_soma) : amp * meter ** -2
    m0=1/(1+exp((-V-34.5*mV)/10/mV)) : 1
    dh/dt=1/tauh*(hinf-h) : 1
    hinf=1/(1+exp((V+59.4*mV)/10.7/mV)) : 1
    tauh=0.15*ms+1.15*ms/(1+exp((V+33.5*mV)/15/mV)) : second
IK=gK_IB_soma*m**4*(V-VK_IB_soma) : amp * meter ** -2
    dm/dt=1/taum*(minf-m) : 1
    minf=1/(1+exp((-V-29.5*mV)/10/mV)) : 1
    taum=0.25*ms+4.35*ms*exp(-abs(V+10*mV)/10/mV) : second
    
Iapp=sinp*ginp*(V-Vrev_inp) : amp * meter ** -2
    dsinp/dt=-sinp/taudinp + (1-sinp)/taurinp*0.5*(1+tanh(Vinp/10/mV)) : 1
    dVinp/dt=1/tauinp*(Vlow-Vinp) : volt
    ginp : siemens * meter **-2
'''


##Constants :
C_IB_soma = 0.9* ufarad * cm ** -2
gL_IB_soma=1 * msiemens * cm **-2
VL_IB_soma=-70*mV
gNa_IB_soma=50 * msiemens * cm **-2
VNa_IB_soma=50*mV
gK_IB_soma=10 * msiemens * cm **-2
VK_IB_soma=-95*mV

if __name__=='__main__' :
    start_scope()
    IB_soma=NeuronGroup(1,eq_IB_soma,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    IB_soma.V = '-100*mvolt+10*rand()*mvolt'
    IB_soma.h = '0+0.05*rand()'
    IB_soma.m = '0+0.05*rand()'
    IB_soma.J='-13 * uA * cmeter ** -2'
    
    V1=StateMonitor(IB_soma,'V',record=[0])
    
#    I1=StateMonitor(IB_soma,'IL',record=[0])
#    I2=StateMonitor(IB_soma,'INa',record=[0])
#    I3=StateMonitor(IB_soma,'IK',record=[0])
#    I4=StateMonitor(IB_soma,'IAR',record=[0])
    
    run(1*second)
    
    figure()
    plot(V1.t/second,V1.V[0]/volt)
    xlabel('Time (s)')
    ylabel('Membrane potential (V)')
    title('IB_soma cell')
    
#    figure()
#    plot(I1.t/second,I1.IL[0],label='L')
#    plot(I1.t/second,I2.INa[0],label='Na')
#    plot(I1.t/second,I3.IK[0],label='K')
#    plot(I1.t/second,I4.IAR[0],label='AR')
#    title('Synaptic currents')
#    legend()