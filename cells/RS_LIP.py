# -*- coding: utf-8 -*-


from brian2 import *

defaultclock.dt = 0.01*ms

eq_RS_LIP='''
dV/dt=1/C_RS*(-J-Isyn-Igap-Iran-Iapp-IL-INa-IK-IAR) : volt
J : amp * meter ** -2
Isyn=IsynRS_LIP_sup+IsynFS_LIP_sup+IsynSI_LIP_sup+IsynRS_LIP_gran+IsynFS_LIP_gran+IsynIB_LIP+IsynSI_LIP_deep+Isyn_FEF+Isyn_mdPul : amp * meter ** -2
IsynRS_LIP_sup : amp * meter ** -2
IsynFS_LIP_sup : amp * meter ** -2
IsynSI_LIP_sup : amp * meter ** -2
IsynRS_LIP_gran : amp * meter ** -2
IsynFS_LIP_gran : amp * meter ** -2
IsynIB_LIP : amp * meter ** -2
IsynSI_LIP_deep : amp * meter ** -2
Isyn_FEF : amp * meter ** -2
Isyn_mdPul : amp * meter ** -2
Igap : amp * meter ** -2
Iapp = Iinp1 + Iinp2 + Iinp3 : amp * meter ** -2
IL=gL_RS*(V-VL_RS) : amp * meter ** -2
INa=gNa_RS*m0**3*h*(V-VNa_RS) : amp * meter ** -2
    m0=1/(1+exp((-V-34.5*mV)/10/mV)) : 1
    dh/dt=1/tauh*(hinf-h) : 1
    hinf=1/(1+exp((V+59.4*mV)/10.7/mV)) : 1
    tauh=0.15*ms+1.15*ms/(1+exp((V+33.5*mV)/15/mV)) : second
IK=gK_RS*m**4*(V-VK_RS) : amp * meter ** -2
    dm/dt=1/taum*(minf-m) : 1
    minf=1/(1+exp((-V-29.5*mV)/10/mV)) : 1
    taum=0.25*ms+4.35*ms*exp(-abs(V+10*mV)/10/mV) : second
IAR=gAR_RS*mAR*(V-VAR_RS) : amp * meter ** -2
    dmAR/dt=1/taumAR*(mARinf-mAR) : 1
    mARinf=1/(1+exp((V+87.5*mV)/5.5/mV)) : 1
    taumAR=1*ms/(exp((-14.6*mV-0.086*V)/mV)+exp((-1.87*mV+0.07*V)/mV)) : second

Iran=sig_ranRS*randn()+s_ran*g_ranRS*(V-0*mV) : amp * meter ** -2 (constant over dt)
    ds_ran/dt=-s_ran/4/ms : 1
    
Iinp1=sinp*ginp_RS*(V-Vrev_inp) : amp * meter ** -2
    dsinp/dt=-sinp/taudinp + (1-sinp)/taurinp*0.5*(1+tanh(Vinp/10/mV)) : 1
    dVinp/dt=1/tauinp*(Vlow-Vinp) : volt
    ginp_RS = ginp_RS_good* int(sin(2*pi*t*4*Hz)>0) + ginp_RS_bad* int(sin(2*pi*t*4*Hz)<=0) : siemens * meter **-2
    ginp_RS_good : siemens * meter **-2
    ginp_RS_bad : siemens * meter **-2    
Iinp2: amp * meter ** -2
Iinp3: amp * meter ** -2

Iapp2=sinp2*ginp_RS2*(V-Vrev_inp) : amp * meter ** -2
    dsinp2/dt=-sinp2/taudinp + (1-sinp2)/taurinp*0.5*(1+tanh(Vinp2/10/mV)) : 1
    dVinp2/dt=1/tauinp*(Vlow-Vinp2) : volt
    ginp_RS2 : siemens * meter **-2
'''


##Constants :
C_RS = 0.9* ufarad * cm ** -2
gL_RS=1 * msiemens * cm **-2
VL_RS=-70*mV
gNa_RS=200 * msiemens * cm **-2
VNa_RS=50*mV
gK_RS=20 * msiemens * cm **-2
VK_RS=-95*mV
#gAR_RS=40 * msiemens * cm **-2 #25 in Mark model, but other channel properties have been changed as well (forward,backward rates)
gAR_RS=25 * msiemens * cm **-2
VAR_RS=-35*mV

sig_ranRS=0.15* mamp * cm **-2
sig_ranRS=0.15* mamp * cm **-2*0.5
g_ranRS=0.03* msiemens * cm **-2

#sig_ranRS=0* mamp * cm **-2
#g_ranRS=0* msiemens * cm **-2


if __name__=='__main__' :
    start_scope()
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp=0* msiemens * cm **-2
        
    RS=NeuronGroup(1,eq_RS_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
#    RS.V = '-70*mvolt+10*rand()*mvolt'
#    RS.h = '0+0.05*rand()'
#    RS.m = '0+0.05*rand()'
#    RS.mAR = '0.035+0.025*rand()'
#    RS.J='1 * uA * cmeter ** -2'
#    RS.V = '-60*mvolt'
#    RS.h = '0.56'
#    RS.m = '0.038'
#    RS.mAR = '0.01'
#    RS.J='1 * uA * cmeter ** -2'

    
#    Poisson_input = PoissonGroup(1,0.1/ms)
#    in_syn = Synapses(Poisson_input, RS, on_pre='s_ran+=0.0001') #defaultclock.dt
#    in_syn.connect(j='i')
    
    
    V1=StateMonitor(RS,'V',record=[0])
    
#    I1=StateMonitor(RS,'IL',record=[0])
#    I2=StateMonitor(RS,'INa',record=[0])
#    I3=StateMonitor(RS,'IK',record=[0])
#    I4=StateMonitor(RS,'IAR',record=[0])
    I5=StateMonitor(RS,'J',record=[0])
    
    M1=StateMonitor(RS,'m0',record=[0])
    M2=StateMonitor(RS,'h',record=[0])
    M3=StateMonitor(RS,'m',record=[0])
    M4=StateMonitor(RS,'mAR',record=[0])
    
    run(1*second)
    
    figure()
    plot(V1.t/second,V1.V[0]/volt)
    xlabel('Time (s)')
    ylabel('Membrane potential (V)')
    title('RS cell')
    
#    figure()
#    plot(I1.t/second,I1.IL[0],label='L')
#    plot(I1.t/second,I2.INa[0],label='Na')
#    plot(I1.t/second,I3.IK[0],label='K')
#    plot(I1.t/second,I4.IAR[0],label='AR')
#    title('Synaptic currents')
#    legend()
    
#    figure()
#    plot(I5.t/second,I5.J[0])
    
    figure()
    plot(M1.t/second,M1.m0[0],label='m0')
    plot(M1.t/second,M2.h[0],label='h')
    plot(M1.t/second,M3.m[0],label='m')
    plot(M1.t/second,M4.mAR[0],label='mAR')
    title('Gating variables')
    legend()