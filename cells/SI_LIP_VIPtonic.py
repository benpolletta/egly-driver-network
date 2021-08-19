# -*- coding: utf-8 -*-


from brian2 import *

defaultclock.dt = 0.01*ms

eq_SI_LIP='''
dV/dt=1/C_SI*(-J-Isyn-Igap-Iran-Iapp-Iapp1-IL-INa-IK-IAR) : volt
J : amp * meter ** -2
Isyn=IsynRS_LIP_sup+IsynFS_LIP_sup+IsynSI_LIP_sup+IsynRS_LIP_gran+IsynFS_LIP_gran+IsynIB_LIP+IsynSI_LIP_deep+IsynVIP+Isyn_FEF+Isyn_mdPul : amp * meter ** -2
IsynRS_LIP_sup : amp * meter ** -2
IsynFS_LIP_sup : amp * meter ** -2
IsynSI_LIP_sup : amp * meter ** -2
IsynRS_LIP_gran : amp * meter ** -2
IsynFS_LIP_gran : amp * meter ** -2
IsynIB_LIP : amp * meter ** -2
IsynSI_LIP_deep : amp * meter ** -2
IsynVIP : amp * meter ** -2
Isyn_FEF : amp * meter ** -2
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
    
Iran=sig_ranSI*randn(): amp * meter ** -2 (constant over dt)

Iapp=sinp_SI(t)*ginp_SI*(V-VK_SI) : amp * meter ** -2
     ginp_SI : siemens * meter **-2

Iapp1=sinp*ginp_SI1*(V-VK_SI) : amp * meter ** -2
    dsinp/dt=-sinp/taudinp2 + (1-sinp)/taurinp2*0.5*(1+tanh(Vinp/10/mV)) : 1
    dVinp/dt=1/tauinp2*(Vlow-Vinp) : volt
    ginp_SI1 : siemens * meter **-2
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

taurinp2=0.2*ms
taudinp2=5*ms
tauinp2=taudinp2

def generate_spike_timing(N,f,start_time,end_time):
    list_time_and_i=[]
    for i in range(N):
        list_time=[(start_time,i)]
        next_spike=list_time[-1][0]+(1+0.001*rand())/f
        while next_spike<end_time:
            list_time.append((next_spike,i))
            next_spike=list_time[-1][0]+(1+0.001*rand())/f
        list_time_and_i+=list_time
    return array(list_time_and_i)

if __name__=='__main__' :
    start_scope()
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp_IB=0* msiemens * cm **-2
    ginp_SI=20* msiemens * cm **-2
        
    VIP_timing=generate_spike_timing(1,50*Hz,500*ms,600*ms)
    VIPinput = SpikeGeneratorGroup(1, VIP_timing[:,1], VIP_timing[:,0]*second)
    
    #sinp_SI=TimedArray([0,1,0], dt=667*ms)
    
    SI=NeuronGroup(1,eq_SI_LIP,threshold='V>-20*mvolt',refractory=3*ms,method='rk4')
    SI.V = '-100*mvolt+10*rand()*mvolt'
    SI.h = '0+0.05*rand()'
    SI.m = '0+0.05*rand()'
    SI.mAR = '0.02+0.04*rand()'
    SI.J='-20 * uA * cmeter ** -2'
    SI.ginp_SI='20 * msiemens * cm ** -2'
    
    S_sinp_SI = Synapses(VIPinput,SI,on_pre='Vinp=Vhigh')
    S_sinp_SI.connect(j='i')
    
    V1=StateMonitor(SI,'V',record=[0])
    Inp=SpikeMonitor(VIPinput,record=True)
    #Inp=StateMonitor(SI,'Iapp',record=[0])
    
    net = Network(collect())
    net.add((SI,VIPinput,S_sinp_SI,V1,Inp))
    
#    I1=StateMonitor(SI,'IL',record=[0])
#    I2=StateMonitor(SI,'INa',record=[0])
#    I3=StateMonitor(SI,'IK',record=[0])
    
    net.run(1*second,report='text',report_period=300*second)
    
    figure()
    subplot(211)
    plot(V1.t/second,V1.V[0]/volt)
    xlabel('Time (s)')
    ylabel('Membrane potential (V)')
    title('SI cell')
    subplot(212)
    plot(Inp.t,Inp.i,'k.')
    #plot(Inp.t/second,Inp.Iapp[0]/volt)
    xlabel('Time (s)')
    ylabel('Input Current')
    xlim([0,1])
    
    
#    figure()
#    plot(I1.t/second,I1.IL[0],label='L')
#    plot(I1.t/second,I2.INa[0],label='Na')
#    plot(I1.t/second,I3.IK[0],label='K')
#    plot(I1.t/second,I4.IAR[0],label='AR')
#    title('Synaptic currents')
#    legend()