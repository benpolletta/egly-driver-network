# -*- coding: utf-8 -*-


from brian2 import *

defaultclock.dt = 0.01*ms

eq_IB_bd='''
dV/dt=1/C_IB_bd*(-J-Isyn-Igap-Iran-Iapp-IL-INa-IK-IAR-IKM-ICaH) : volt
J : amp * meter ** -2
Isyn=IsynRS_LIP_sup+IsynFS_LIP_sup+IsynSI_LIP_sup+IsynRS_LIP_gran+IsynFS_LIP_gran+IsynIB_LIP+IsynSI_LIP_deep+IsynSI2_LIP_deep+Isyn_FEF+Isyn_mdPul : amp * meter ** -2
IsynRS_LIP_sup : amp * meter ** -2
IsynFS_LIP_sup : amp * meter ** -2
IsynSI_LIP_sup : amp * meter ** -2
IsynRS_LIP_gran : amp * meter ** -2
IsynFS_LIP_gran : amp * meter ** -2
IsynIB_LIP : amp * meter ** -2
IsynSI_LIP_deep : amp * meter ** -2
IsynSI2_LIP_deep : amp * meter ** -2
Isyn_FEF : amp * meter ** -2
Isyn_mdPul : amp * meter ** -2
Igap = Igap_soma+Igap_axon+Igap_ad+Igap_bd : amp * meter ** -2
Igap_soma : amp * meter ** -2
Igap_axon : amp * meter ** -2
Igap_ad : amp * meter ** -2
Igap_bd : amp * meter ** -2
IL=gL_IB_bd*(V-VL_IB_bd) : amp * meter ** -2
INa=gNa_IB_bd*m0**3*h*(V-VNa_IB_bd) : amp * meter ** -2
    m0=1/(1+exp((-V-34.5*mV)/10/mV)) : 1
    dh/dt=1/tauh*(hinf-h) : 1
    hinf=1/(1+exp((V+59.4*mV)/10.7/mV)) : 1
    tauh=0.15*ms+1.15*ms/(1+exp((V+33.5*mV)/15/mV)) : second
IK=gK_IB_bd*m**4*(V-VK_IB_bd) : amp * meter ** -2
    dm/dt=1/taum*(minf-m) : 1
    minf=1/(1+exp((-V-29.5*mV)/10/mV)) : 1
    taum=0.25*ms+4.35*ms*exp(-abs(V+10*mV)/10/mV) : second
IAR=gAR_IB_bd*mAR*(V-VAR_IB_bd) : amp * meter ** -2
    dmAR/dt=1/taumAR*(mARinf-mAR) : 1
    mARinf=1/(1+exp((V+75*mV)/5.5/mV)) : 1
    taumAR=1*ms/(exp((-14.6*mV-0.086*V)/mV)+exp((-1.87*mV+0.07*V)/mV)) : second
IKM=gKM_IB_bd*mKM*(V-VKM_IB_bd) : amp * meter ** -2
    dmKM/dt=alphaKM*(1-mKM)-betaKM*mKM : 1
    alphaKM= 0.02/(1+exp((-V-20*mV)/5/mV))/ms : hertz
    betaKM= 0.01*exp((-V-43*mV)/18/mV)/ms: hertz
ICaH=gCaH_IB_bd*mKM**2*(V-VCaH_IB_bd) : amp * meter ** -2
    dmCaH/dt=alphaCaH*(1-mCaH)-betaCaH*mCaH : 1
    alphaCaH= 1.6/(1+exp(-0.072*(-V-5*mV)/mV))/ms : hertz
    betaCaH= 0.02/ms*(V+8.9*mV)/mV/(exp((V+8.9*mV)/5/mV)-1): hertz
    
Iran=sig_ranIB_bd*randn(): amp * meter ** -2 (constant over dt)

Iapp=sinp*ginp_IB*(V-Vrev_inp) : amp * meter ** -2
    dsinp/dt=-sinp/taudinp2 + (1-sinp)/taurinp2*0.5*(1+tanh(Vinp/10/mV)) : 1
    dVinp/dt=1/tauinp2*(Vlow-Vinp) : volt
    ginp_IB : siemens * meter **-2
'''


##Constants :
C_IB_bd = 0.9* ufarad * cm ** -2
gL_IB_bd=2 * msiemens * cm **-2
VL_IB_bd=-70*mV
gNa_IB_bd=125 * msiemens * cm **-2
VNa_IB_bd=50*mV
gK_IB_bd=10 * msiemens * cm **-2
VK_IB_bd=-95*mV
gAR_IB_bd=115 * msiemens * cm **-2
VAR_IB_bd=-25*mV
gKM_IB_bd=0.75 * msiemens * cm **-2
VKM_IB_bd=-95*mV
gCaH_IB_bd=6.5 * msiemens * cm **-2
VCaH_IB_bd=125*mV

sig_ranIB_bd=0.005* mamp * cm **-2
sig_ranIB_bd=0.005* mamp * cm **-2*0.5

taurinp2=0.1*ms
taudinp2=0.5*ms
tauinp2=taudinp2


if __name__=='__main__' :
    start_scope()
    IB_bd=NeuronGroup(1,eq_IB_bd,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    IB_bd.V = '-100*mvolt+10*rand()*mvolt'
    IB_bd.h = '0+0.05*rand()'
    IB_bd.m = '0+0.05*rand()'
    IB_bd.mAR = '0+0.001*rand()'
    IB_bd.mKM = '0+0.05*rand()'
    IB_bd.mCaH = '0+0.01*rand()'
    IB_bd.J='-13 * uA * cmeter ** -2'
    
    V1=StateMonitor(IB_bd,'V',record=[0])
    
#    I1=StateMonitor(IB_bd,'IL',record=[0])
#    I2=StateMonitor(IB_bd,'INa',record=[0])
#    I3=StateMonitor(IB_bd,'IK',record=[0])
#    I4=StateMonitor(IB_bd,'IAR',record=[0])
    
    run(1*second)
    
    figure()
    plot(V1.t/second,V1.V[0]/volt)
    xlabel('Time (s)')
    ylabel('Membrane potential (V)')
    title('IB_bd cell')
    
#    figure()
#    plot(I1.t/second,I1.IL[0],label='L')
#    plot(I1.t/second,I2.INa[0],label='Na')
#    plot(I1.t/second,I3.IK[0],label='K')
#    plot(I1.t/second,I4.IAR[0],label='AR')
#    title('Synaptic currents')
#    legend()