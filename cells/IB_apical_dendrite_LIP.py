# -*- coding: utf-8 -*-


from brian2 import *

defaultclock.dt = 0.01*ms

eq_IB_ad='''
dV/dt=1/C_IB_ad*(-J-Isyn-Igap-Iran-Iapp-IL-INa-IK-IAR-IKM-ICaH) : volt
J : amp * meter ** -2
Isyn=IsynRS_LIP_sup_AMPA+IsynRS_LIP_sup_NMDA+IsynFS_LIP_sup+IsynSI_LIP_sup+IsynRS_LIP_gran+IsynFS_LIP_gran+IsynIB_LIP+IsynSI_LIP_deep+Isyn_FEF+Isyn_mdPul : amp * meter ** -2
IsynRS_LIP_sup_AMPA : amp * meter ** -2
IsynRS_LIP_sup_NMDA : amp * meter ** -2
IsynFS_LIP_sup : amp * meter ** -2
IsynSI_LIP_sup : amp * meter ** -2
IsynRS_LIP_gran : amp * meter ** -2
IsynFS_LIP_gran : amp * meter ** -2
IsynIB_LIP : amp * meter ** -2
IsynSI_LIP_deep : amp * meter ** -2
Isyn_FEF : amp * meter ** -2
Isyn_mdPul : amp * meter ** -2
Igap = Igap_soma+Igap_axon+Igap_ad+Igap_bd : amp * meter ** -2
Igap_soma : amp * meter ** -2
Igap_axon : amp * meter ** -2
Igap_ad : amp * meter ** -2
Igap_bd : amp * meter ** -2
IL=gL_IB_ad*(V-VL_IB_ad) : amp * meter ** -2
INa=gNa_IB_ad*m0**3*h*(V-VNa_IB_ad) : amp * meter ** -2
    m0=1/(1+exp((-V-34.5*mV)/10/mV)) : 1
    dh/dt=1/tauh*(hinf-h) : 1
    hinf=1/(1+exp((V+59.4*mV)/10.7/mV)) : 1
    tauh=0.15*ms+1.15*ms/(1+exp((V+33.5*mV)/15/mV)) : second
IK=gK_IB_ad*m**4*(V-VK_IB_ad) : amp * meter ** -2
    dm/dt=1/taum*(minf-m) : 1
    minf=1/(1+exp((-V-29.5*mV)/10/mV)) : 1
    taum=0.25*ms+4.35*ms*exp(-abs(V+10*mV)/10/mV) : second
IAR=gAR_IB_ad*mAR*(V-VAR_IB_ad) : amp * meter ** -2
    dmAR/dt=1/taumAR*(mARinf-mAR) : 1
    mARinf=1/(1+exp((V+75*mV)/5.5/mV)) : 1
    taumAR=1*ms/(exp((-14.6*mV-0.086*V)/mV)+exp((-1.87*mV+0.07*V)/mV)) : second
IKM=gKM_IB_ad*mKM*(V-VKM_IB_ad) : amp * meter ** -2
    dmKM/dt=alphaKM*(1-mKM)-betaKM*mKM : 1
    alphaKM= 0.02/(1+exp((-V-20*mV)/5/mV))/ms : hertz
    betaKM= 0.01*exp((-V-43*mV)/18/mV)/ms: hertz
ICaH=gCaH_IB_ad*mKM**2*(V-VCaH_IB_ad) : amp * meter ** -2
    dmCaH/dt=alphaCaH*(1-mCaH)-betaCaH*mCaH : 1
    alphaCaH= 1.6/(1+exp(-0.072*(-V-5*mV)/mV))/ms : hertz
    betaCaH= 0.02/ms*(V+8.9*mV)/mV/(exp((V+8.9*mV)/5/mV)-1): hertz
    
Iran=sig_ranIB_ad*randn(): amp * meter ** -2 (constant over dt)

Iapp=sinp*ginp*(V-Vrev_inp) : amp * meter ** -2
    dsinp/dt=-sinp/taudinp + (1-sinp)/taurinp*0.5*(1+tanh(Vinp/10/mV)) : 1
    dVinp/dt=1/tauinp*(Vlow-Vinp) : volt
    ginp : siemens * meter **-2
'''


##Constants :
C_IB_ad = 0.9* ufarad * cm ** -2
gL_IB_ad=2 * msiemens * cm **-2
VL_IB_ad=-70*mV
gNa_IB_ad=125 * msiemens * cm **-2
VNa_IB_ad=50*mV
gK_IB_ad=10 * msiemens * cm **-2
VK_IB_ad=-95*mV
gAR_IB_ad=155 * msiemens * cm **-2
VAR_IB_ad=-25*mV
gKM_IB_ad=0.75 * msiemens * cm **-2
VKM_IB_ad=-95*mV
gCaH_IB_ad=6.5 * msiemens * cm **-2
VCaH_IB_ad=125*mV

sig_ranIB_ad=0.005* mamp * cm **-2
sig_ranIB_ad=0.005* mamp * cm **-2*0.5

if __name__=='__main__' :
    start_scope()
    IB_ad=NeuronGroup(1,eq_IB_ad,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    IB_ad.V = '-100*mvolt+10*rand()*mvolt'
    IB_ad.h = '0+0.05*rand()'
    IB_ad.m = '0+0.05*rand()'
    IB_ad.mAR = '0+0.001*rand()'
    IB_ad.mKM = '0+0.05*rand()'
    IB_ad.mCaH = '0+0.01*rand()'
    IB_ad.J='-13 * uA * cmeter ** -2'
    
    V1=StateMonitor(IB_ad,'V',record=[0])
    
#    I1=StateMonitor(IB_ad,'IL',record=[0])
#    I2=StateMonitor(IB_ad,'INa',record=[0])
#    I3=StateMonitor(IB_ad,'IK',record=[0])
#    I4=StateMonitor(IB_ad,'IAR',record=[0])
    
    run(1*second)
    
    figure()
    plot(V1.t/second,V1.V[0]/volt)
    xlabel('Time (s)')
    ylabel('Membrane potential (V)')
    title('IB_ad cell')
    
#    figure()
#    plot(I1.t/second,I1.IL[0],label='L')
#    plot(I1.t/second,I2.INa[0],label='Na')
#    plot(I1.t/second,I3.IK[0],label='K')
#    plot(I1.t/second,I4.IAR[0],label='AR')
#    title('Synaptic currents')
#    legend()