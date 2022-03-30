# -*- coding: utf-8 -*-


from brian2 import *

defaultclock.dt = 0.01*ms

eq_SI_FEF='''
dV/dt=1/C_SI*(-J-Isyn-Igap-IL-INa-IK-IAR) : volt
J : amp * meter ** -2
Isyn=IsynRS_FEF_VM+IsynSI_FEF_VM+IsynSI2_FEF_VM+IsynRS_FEF_V+IsynFS_FEF_V+Isyn_LIP+Isyn_mdPul : amp * meter ** -2
IsynRS_FEF_VM : amp * meter ** -2
IsynSI_FEF_VM : amp * meter ** -2
IsynSI2_FEF_VM : amp * meter ** -2
IsynRS_FEF_V : amp * meter ** -2
IsynFS_FEF_V : amp * meter ** -2
Isyn_LIP : amp * meter ** -2
Isyn_mdPul : amp * meter ** -2

Igap : amp * meter ** -2
IL=gL_SI*(V-VL_SI) : amp * meter ** -2
INa=gNa_SI*m0**3*h*(V-VNa_SI) : amp * meter ** -2
    m0=1/(1+exp((-V-38*mV)/10/mV)) : 1
    dh/dt=1/tauh*(hinf-h) : 1
    hinf=1/(1+exp((V+58.3*mV)/6.7/mV)) : 1
    tauh=0.225*ms+1.125*ms/(1+exp((V+37*mV)/15/mV)) : second
IK=gK_SI*m**4*(V-VK_SI) : amp * meter ** -2
    dm/dt=1/taum*(minf-m) : 1
    minf=1/(1+exp((-V-27*mV)/11.5/mV)) : 1
    taum=0.25*ms+4.35*ms*exp(-abs(V+10*mV)/10/mV) : second
IAR=gAR_SI*mAR*(V-VAR_SI) : amp * meter ** -2
    dmAR/dt=1/taumAR*(mARinf-mAR) : 1
    mARinf=1/(1+exp((V+75*mV)/5.5/mV)) : 1
    taumAR=1*ms/(exp((-14.6*mV-0.086*V)/mV)+exp((-1.87*mV+0.07*V)/mV)) : second
'''


##Constants :
C_SI = 0.9* ufarad * cm ** -2
gL_SI=6 * msiemens * cm **-2
VL_SI=-65*mV

gNa_SI=200 * msiemens * cm **-2
VNa_SI=50*mV

gK_SI=10 * msiemens * cm **-2
VK_SI=-100*mV

gAR_SI=50 * msiemens * cm **-2
VAR_SI=-35*mV

sig_ranSI=0.05* mamp * cm **-2
sig_ranSI=0.05* mamp * cm **-2*0.5
sig_ranSI=0.05* mamp * cm **-2*0

if __name__=='__main__' :
    start_scope()
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp_IB=0* msiemens * cm **-2
    ginp=0* msiemens * cm **-2
    
    SI=NeuronGroup(100,eq_SI_FEF,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    SI.V = '-70*mvolt'
    SI.h = '0'
    SI.m = '0'
    SI.mAR = '0.02'
    SI.J='-i/100 *20 * uA * cmeter ** -2'
    
    V1=StateMonitor(SI,'V',record=[0])
    R1=SpikeMonitor(SI,record=True)
    
#    I1=StateMonitor(SI,'IL',record=[0])
#    I2=StateMonitor(SI,'INa',record=[0])
#    I3=StateMonitor(SI,'IK',record=[0])
    
    run(10*second)
    
    figure()
    plot(V1.t/second,V1.V[0]/volt)
    xlabel('Time (s)')
    ylabel('Membrane potential (V)')
    title('SI cell')
    
    figure()
    plot((-SI.J/ (uA * cmeter ** -2)),R1.count/10)
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