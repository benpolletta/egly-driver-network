# -*- coding: utf-8 -*-


from brian2 import *

defaultclock.dt = 0.01*ms

eq_RS_FEF='''
dV/dt=1/C_RS*(-J-IL-INa-IK-IAR) : volt
J : amp * meter ** -2
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

sig_ranRS=0.15* mamp * cm **-2*0
g_ranRS=0.03* msiemens * cm **-2*0

#gL_RS=0.9 * msiemens * cm **-2
#gNa_RS=200 * msiemens * cm **-2
#gK_RS=20 * msiemens * cm **-2
#sig_ranRS=0* mamp * cm **-2
#g_ranRS=0* msiemens * cm **-2


if __name__=='__main__' :
    start_scope()
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    taurinp3=0.1*ms
    taudinp3=0.5*ms
    tauinp3=taudinp3
    Vhigh=0*mV
    Vlow=-80*mV
    ginp=0* msiemens * cm **-2
        
    RS=NeuronGroup(100,eq_RS_FEF,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
#    RS.V = '-70*mvolt+10*rand()*mvolt'
#    RS.h = '0+0.05*rand()'
#    RS.m = '0+0.05*rand()'
#    RS.mAR = '0.035+0.025*rand()'
#    RS.J='1 * uA * cmeter ** -2'
    RS.V = '-60*mvolt'
    RS.h = '0.56'
    RS.m = '0.038'
    RS.mAR = '0.01'
    RS.J='-i/100 *20 * uA * cmeter ** -2'

    
#    Poisson_input = PoissonGroup(1,0.1/ms)
#    in_syn = Synapses(Poisson_input, RS, on_pre='s_ran+=0.0001') #defaultclock.dt
#    in_syn.connect(j='i')
    
    
    V1=StateMonitor(RS,'V',record=[0])
    R1=SpikeMonitor(RS,record=True)
    
#    I1=StateMonitor(RS,'IL',record=[0])
#    I2=StateMonitor(RS,'INa',record=[0])
#    I3=StateMonitor(RS,'IK',record=[0])
#    I4=StateMonitor(RS,'IAR',record=[0])
    I5=StateMonitor(RS,'J',record=[0])
    
    M1=StateMonitor(RS,'m0',record=[0])
    M2=StateMonitor(RS,'h',record=[0])
    M3=StateMonitor(RS,'m',record=[0])
    M4=StateMonitor(RS,'mAR',record=[0])
    
    run(10*second)
    
    figure()
    plot(V1.t/second,V1.V[0]/volt)
    xlabel('Time (s)')
    ylabel('Membrane potential (V)')
    title('RS cell')
    tick_params(axis='both', which='major', labelsize=12)
    
    figure()
    plot((-RS.J/ (uA * cmeter ** -2)),R1.count/10)
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
    
#    figure()
#    plot(I5.t/second,I5.J[0])
    
#    figure()
#    plot(M1.t/second,M1.m0[0],label='m0')
#    plot(M1.t/second,M2.h[0],label='h')
#    plot(M1.t/second,M3.m[0],label='m')
#    plot(M1.t/second,M4.mAR[0],label='mAR')
#    title('Gating variables')
#    legend()