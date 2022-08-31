# -*- coding: utf-8 -*-


from brian2 import *

defaultclock.dt = 0.01*ms

eq_IB_axon='''
dV/dt=1/C_IB_axon*(-J-Isyn-Igap-Iran-Iapp-IL-INa-IK-IKM) : volt
J : amp * meter ** -2
Isyn : amp * meter ** -2
Igap = Igap_soma+Igap_axon+Igap_ad+Igap_bd : amp * meter ** -2
Igap_soma : amp * meter ** -2
Igap_axon : amp * meter ** -2
Igap_ad : amp * meter ** -2
Igap_bd : amp * meter ** -2
IL=gL_IB_axon*(V-VL_IB_axon) : amp * meter ** -2
INa=gNa_IB_axon*m0**3*h*(V-VNa_IB_axon) : amp * meter ** -2
    m0=1/(1+exp((-V-34.5*mV)/10/mV)) : 1
    dh/dt=1/tauh*(hinf-h) : 1
    hinf=1/(1+exp((V+59.4*mV)/10.7/mV)) : 1
    tauh=0.15*ms+1.15*ms/(1+exp((V+33.5*mV)/15/mV)) : second
IK=gK_IB_axon*m**4*(V-VK_IB_axon) : amp * meter ** -2
    dm/dt=1/taum*(minf-m) : 1
    minf=1/(1+exp((-V-29.5*mV)/10/mV)) : 1
    taum=0.25*ms+4.35*ms*exp(-abs(V+10*mV)/10/mV) : second
IKM=gKM_IB_axon*mKM*(V-VKM_IB_axon) : amp * meter ** -2
    dmKM/dt=alphaKM*(1-mKM)-betaKM*mKM : 1
    alphaKM= 0.02/(1+exp((-V-20*mV)/5/mV))/ms : hertz
    betaKM= 0.01*exp((-V-43*mV)/18/mV)/ms: hertz
    
Iran=sig_ranIB_axon*randn(): amp * meter ** -2 (constant over dt)

Iapp=sinp*ginp*(V-Vrev_inp) : amp * meter ** -2
    dsinp/dt=-sinp/taudinp + (1-sinp)/taurinp*0.5*(1+tanh(Vinp/10/mV)) : 1
    dVinp/dt=1/tauinp*(Vlow-Vinp) : volt
    ginp : siemens * meter **-2
'''


##Constants :
C_IB_axon = 0.9* ufarad * cm ** -2
gL_IB_axon=0.25 * msiemens * cm **-2
VL_IB_axon=-70*mV
gNa_IB_axon=100 * msiemens * cm **-2
VNa_IB_axon=50*mV
gK_IB_axon=5 * msiemens * cm **-2
VK_IB_axon=-95*mV
gKM_IB_axon=1.5 * msiemens * cm **-2
VKM_IB_axon=-95*mV

sig_ranIB_axon=0.025* mamp * cm **-2
sig_ranIB_axon=0.025* mamp * cm **-2*0.5

if __name__=='__main__' :
    start_scope()
    IB_axon=NeuronGroup(1,eq_IB_axon,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    IB_axon.V = '-100*mvolt+10*rand()*mvolt'
    IB_axon.h = '0+0.05*rand()'
    IB_axon.m = '0+0.05*rand()'
    IB_axon.mKM = '0+0.05*rand()'
    IB_axon.J='-1 * uA * cmeter ** -2'
    
    V1=StateMonitor(IB_axon,'V',record=[0])
    
#    I1=StateMonitor(IB_axon,'IL',record=[0])
#    I2=StateMonitor(IB_axon,'INa',record=[0])
#    I3=StateMonitor(IB_axon,'IK',record=[0])
#    I4=StateMonitor(IB_axon,'IAR',record=[0])
    
    run(1*second)
    
    figure()
    plot(V1.t/second,V1.V[0]/volt)
    xlabel('Time (s)')
    ylabel('Membrane potential (V)')
    title('IB_axon cell')
    
#    figure()
#    plot(I1.t/second,I1.IL[0],label='L')
#    plot(I1.t/second,I2.INa[0],label='Na')
#    plot(I1.t/second,I3.IK[0],label='K')
#    plot(I1.t/second,I4.IAR[0],label='AR')
#    title('Synaptic currents')
#    legend()