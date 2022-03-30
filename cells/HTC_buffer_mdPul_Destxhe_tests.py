# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 10:46:47 2019

@author: Amelie Aussel
"""

from brian2 import *
defaultclock.dt = 0.01*ms
prefs.codegen.target = 'numpy'

eq_HTC_buffer_mdPul='''
    dV/dt=(-INa-IK-IL-IKL-IH-ITLT-ITHT-IAHP-IGJ-Iapp-Iapp2-J-Isyn-Itheta)/Cm_HTC : volt
    
    J : amp * meter ** -2

    Itheta = amptheta * int(sin(2*pi*time_total*ftheta+offsettheta)<0) : amp * meter ** -2
    amptheta : amp * meter ** -2
    ftheta : Hz    
    dtime_total/dt=1 : second
    offsettheta : 1
    
    Iapp = randn() * gapp : amp * meter ** -2 (constant over dt)

    Isyn=IsynHTC+IsynTC+IsynI+IsynREintGABAA+IsynREintGABAB+IsynREA+IsynREB+IsynREC+Isyn_FEF+Isyn_LIP : amp * meter ** -2
    IsynHTC : amp * meter ** -2
    IsynTC : amp * meter ** -2
    IsynI : amp * meter ** -2
    IsynREintGABAA : amp * meter ** -2
    IsynREintGABAB : amp * meter ** -2
    IsynREA : amp * meter ** -2
    IsynREB : amp * meter ** -2
    IsynREC : amp * meter ** -2
    Isyn_FEF : amp * meter ** -2
    Isyn_LIP : amp * meter ** -2
    
    Vt=V+25*mV : volt
    INa=gNa_HTC*mNa**3*hNa*(V-ENa_HTC)*int(mNa>0)*int(hNa>0) : amp * meter ** -2
        dmNa/dt = (mNainf-mNa)/taumNa : 1
        dhNa/dt = (hNainf-hNa)/tauhNa : 1
        mNainf = alphamNa/(alphamNa+betamNa) : 1
        taumNa=1*ms/(alphamNa+betamNa) : second      
        hNainf = alphahNa/(alphahNa+betahNa) : 1
        tauhNa=1*ms/(alphahNa+betahNa) : second
        alphamNa=(0.32*(13*mV-Vt))/(exp((13*mV-Vt)/4/mV)-1)/mV : 1
        betamNa = (0.28*(Vt-40*mV))/(exp((Vt-40*mV)/5/mV)-1)/mV : 1
        alphahNa=(0.128*exp((17*mV-Vt)/18/mV)) : 1
        betahNa = 4/(1+exp((40*mV-Vt)/5/mV)) : 1 
        
    IK=gK_HTC*mK**4*(V-EK_HTC) : amp * meter ** -2
        dmK/dt=(mKinf-mK)/taumK : 1
        mKinf = alphamK/(alphamK+betamK) : 1
        taumK=1*ms/(alphamK+betamK) : second
        alphamK = (0.032*(15*mV-Vt)/mV)/(exp((15*mV-Vt)/5/mV)-1) : 1
        betamK=0.5 * exp((10*mV-Vt)/40/mV) : 1
        
    IL=gL_HTC*(V-EL_HTC) : amp * meter ** -2
    IKL=gKL_HTC*(V-EKL_HTC) : amp * meter ** -2
    
    ITLT = gTLT_HTC * mTLT**2 * hTLT * (V-ETLT) : amp * meter ** -2
        dmTLT/dt = (mTLTinf-mTLT)/taumTLT : 1
        dhTLT/dt = (hTLTinf-hTLT)/tauhTLT : 1 
        mTLTinf = 1/(1+exp(-(V+59*mV)/6.2/mV)) : 1
        taumTLT= 1*ms: second    
        hTLTinf = 1/(1+exp((V+83*mV)/4/mV)) : 1
        tauhTLT=  (30.8+ (211.4 + exp((V+115.2*mV)/5/mV))/(1+exp((V+86*mV)/3.2/mV)))/3.7372*ms : second
        dCa_iTLT/dt = ((-10000000*ITLT)/(2*96485.3* coulomb * mole ** -1 * meter)- (Ca_iTLT-0.00024*mmolar)/5/ms ) : mole * meter**-3
        ETLT=R*T/(z*F)*log(2*mmolar/Ca_iTLT)/log(2) : volt
    
    
    ITHT = gTHT_HTC * mTHT**2 * hTHT * (V-ETHT) : amp * meter ** -2
        dmTHT/dt = (mTHTinf-mTHT)/taumTHT : 1
        dhTHT/dt = (hTHTinf-hTHT)/tauhTHT : 1 
        mTHTinf = 1/(1+exp(-(V+40.1*mV)/3.5/mV)) : 1
        taumTHT= 1*ms: second    
        hTHTinf = 1/(1+exp((V+62.2*mV)/5.5/mV)) : 1
        tauhTHT=  0.1483*ms * exp(-0.09398 * V/mV) + 5.284*ms * exp(0.008855 * V/mV) : second
        dCa_iTHT/dt = 0.4*(-24*ITHT)*umolar/msecond/(mamp * cm **-2)- Ca_iTHT/100/ms: mole * meter**-3
        ETHT=120*mvolt : volt

        
    
    IH = gH_HTC * r * (V-EH_HTC) : amp * meter ** -2
        dr/dt=(rinf-r)/tausr : 1
        tausr = 1*(20*ms + 1000*ms / (exp((V+71.5*mV)/14.2/mV)+exp(-(V+89*mV)/11.6/mV))) : second
        rinf = 1/(1+exp((V+75*mV)/5.5/mV)) : 1
        
    IAHP = gAHP_HTC * mAHP**2  * (V-EAHP) : amp * meter ** -2
        dmAHP/dt = (mAHPinf-mAHP)/taumAHP : 1
        aKCa= 0.1*(Ca_iTHT/(mmolar))*int(0.1*(Ca_iTHT/(mmolar))<=1)+1*int(0.1*(Ca_iTHT/(mmolar))>1): 1
        mAHPinf= aKCa / (aKCa + 0.002) : 1
        taumAHP = 1*ms/(aKCa+ 0.002) : second
    
    IGJ : amp * meter ** -2
    
    Iapp2=sinp*ginp_HTC*(V-Vrev_inp) : amp * meter ** -2
        dsinp/dt=-sinp/taudinp + (1-sinp)/taurinp*0.5*(1+tanh(Vinp/10/mV)) : 1
        dVinp/dt=1/tauinp*(Vlow-Vinp) : volt
        ginp_HTC : siemens * metre **-2   
        
    buffer_pointer : integer (shared)
    delay_steps : integer
    voltage_delayed:volt
    '''
    
Cm_HTC = 2.5* ufarad * cm ** -2
gNa_HTC=90e-3 * siemens * cm **-2
ENa_HTC=50*mV
gK_HTC=10e-3 * siemens * cm **-2
#gK_HTC=5e-3 * siemens * cm **-2
EK_HTC=-100*mV
gL_HTC=0.010e-3 * siemens * cm **-2  
#gL_HTC=0.001e-3 * siemens * cm **-2  
EL_HTC=-70*mV
gKL_HTC=0.0069e-3 * siemens * cm **-2   
#gKL_HTC=0.001e-3 * siemens * cm **-2  
#gKL_HTC=0e-3 * siemens * cm **-2  
EKL_HTC=-100*mV

#gTLT_HTC= 2.1e-3 * siemens * cm **-2
#gTLT_HTC= 15e-3 * siemens * cm **-2
#gTLT_HTC= 5e-3 * siemens * cm **-2
gTLT_HTC= 15e-3 * siemens * cm **-2

#gTHT_HTC= 6e-3 * siemens * cm **-2
gTHT_HTC= 3e-3 * siemens * cm **-2
#gTHT_HTC= 1e-3 * siemens * cm **-2

#gAHP_HTC= 15e-3 * siemens * cm **-2
#gAHP_HTC= 3e-3 * siemens * cm **-2
gAHP_HTC= 2.9e-3 * siemens * cm **-2
#EAHP=-95*mV
EAHP=-80*mV

#gH_HTC = 0.36e-3 * siemens * cm **-2
gH_HTC = 0.5e-3 * siemens * cm **-2
#gH_HTC = 0.1e-3 * siemens * cm **-2

EH_HTC = -40 * mV

R = 8.314 * joule * kelvin**-1 * mole **-1
T = (273.15+37.2)*kelvin
z = 2
F = 96485.3* coulomb * mole ** -1

a=1

N_HTC=20


#gapp=0.1*mamp * cmeter ** -2 
#gapp=0.05*mamp * cmeter ** -2 
#gapp=0.4*mamp * cmeter ** -2 
gapp=0*mamp * cmeter ** -2 


if __name__=='__main__':
    Vrev_inp=0*mV
    taurinp=0.1*ms
    taudinp=0.5*ms
    tauinp=taudinp
    Vhigh=0*mV
    Vlow=-80*mV
    ginp_IB=0* msiemens * cm **-2
    ginp=0* msiemens * cm **-2
    
    close('all')
    start_scope()
    HTC=NeuronGroup(N_HTC,eq_HTC_buffer_mdPul,threshold='V>0*mvolt',refractory=3*ms,method='rk4')
    HTC.V = '-25*mvolt'
    HTC.Ca_iTLT = '1e-7 * mole * metre**-3'
    HTC.Ca_iTHT = '1e-7 * mole * metre**-3'
    HTC.Ca_iTHT = '0.01*mmolar'
#    HTC.J='-150 * nA * cmeter ** -2'
#    HTC.J='-1000 * nA * cmeter ** -2'
    HTC.J='0 * nA * cmeter ** -2'
    
    HTC.mAHP='0.3'
    
    HTC.delay_steps = [3999]  # delay in time steps per neuron
    buffer_size = 4000  # 1+Maximum delay (in time steps)
    
    HTC.variables.add_array('voltage_buffer', dimensions=volt.dim, size=(buffer_size, len(HTC)))
    
    update_code = '''buffer_pointer = (buffer_pointer + 1) % buffer_size
                     voltage_delayed = update_voltage_buffer(V, voltage_buffer, buffer_pointer, delay_steps, buffer_size)'''
       
    buffer_updater = HTC.run_regularly(update_code, codeobj_class=NumpyCodeObject)
        
    @check_units(V=volt, voltage_buffer=volt, buffer_pointer=1, delay_steps=1, buffer_size=1, result=volt)
    def update_voltage_buffer(V, voltage_buffer, buffer_pointer, delay_steps, buffer_size):
        # Write current rate into the buffer
        voltage_buffer[buffer_pointer, :] = V
        # Get delayed rates 
        rows = (buffer_pointer - delay_steps) % buffer_size    
        return voltage_buffer[rows, arange(len(rows))]
        
    GJ_HTC=Synapses(HTC,HTC,'''
                 w : siemens * meter **-2 # gap junction conductance
                 IGJ_post = w * (V_post - V_pre) : amp * meter ** -2 (summed)
                 ''')
    GJ_HTC.connect()
    GJ_HTC.w = 0.02e-3 * siemens * cm **-2 #for 20 HTC cells
    GJ_HTC.w = 0.08e-3 * siemens * cm **-2 #for 20 HTC cells
#    GJ_HTC.w = 0.3e-3 * siemens * cm **-2 #for 2 HTC cells
#    GJ_HTC.w = 0e-3 * siemens * cm **-2 #for 1 cell
    
    V1=StateMonitor(HTC,'V',record=[0])
    V2=StateMonitor(HTC,'IAHP',record=[0])
    V3=StateMonitor(HTC,'ITHT',record=[0])
    V4=StateMonitor(HTC,'Ca_iTHT',record=[0])
    V5=StateMonitor(HTC,'mAHP',record=[0])
    
    V6=StateMonitor(HTC,'INa',record=[0])
    V7=StateMonitor(HTC,'IK',record=[0])
    V8=StateMonitor(HTC,'IL',record=[0])
    V9=StateMonitor(HTC,'IKL',record=[0])
    V10=StateMonitor(HTC,'IH',record=[0])
    V11=StateMonitor(HTC,'ITLT',record=[0])
    
    V12=StateMonitor(HTC,'taumAHP',record=[0])
    V13=StateMonitor(HTC,'mAHPinf',record=[0])
    
    V14=StateMonitor(HTC,'mTHT',record=[0])
    V15=StateMonitor(HTC,'hTHT',record=[0])
    
    R1=SpikeMonitor(HTC)
    
    prefs.codegen.target = 'cython'
    run(1*second)
#    run(10*second)
    
    figure()
    plot(V1.t/second,V1.V[0]/volt)
#    plot(V2.t/second,V2.J[0])
    xlabel('Time (s)')
    ylabel('Membrane potential (V)')
    
    figure()
    plot(V2.t/second,V2.IAHP[0]/volt)
#    plot(V2.t/second,V2.J[0])
    xlabel('Time (s)')
    ylabel('IAHP')
    
    figure()
    plot(V3.t/second,V3.ITHT[0])
    xlabel('Time (s)')
    ylabel('ITHT')
    
    figure()
    plot(V4.t/second,V4.Ca_iTHT[0]/mmolar)
    xlabel('Time (s)')
    ylabel('Ca (mmolar)')
    
#    figure()
#    plot(V5.t/second,V5.mAHP[0])
#    xlabel('Time (s)')
#    ylabel('mAHP')
#    
#    figure()
#    plot(V12.t/second,V12.taumAHP[0]/second)
#    xlabel('Time (s)')
#    ylabel('taumAHP(s)')
#        
#    figure()
#    plot(V13.t/second,V13.mAHPinf[0]/second)
#    xlabel('Time (s)')
#    ylabel('mAHPinf')
    
    figure()
    plot(V14.t/second,V14.mTHT[0], label='mTHT')
    plot(V15.t/second,V15.hTHT[0], label='hTHT')
    plot(V15.t/second,V14.mTHT[0]*V14.mTHT[0]*V15.hTHT[0], label='m2*h')
    xlabel('Time (s)')
    legend()
    
#    figure()
#    plot(V15.t/second,V15.hTHT[0])
#    xlabel('Time (s)')
#    ylabel('hTHT')
    
    figure()
    plot(R1.t,R1.i,'r.')
    
    figure()
    plot(V2.t/second,V2.IAHP[0], label='IAHP')
    plot(V3.t/second,V3.ITHT[0],label='ITHT')
    plot(V6.t/second,V6.INa[0],label='INa')
    plot(V7.t/second,V7.IK[0],label='IK')
    plot(V8.t/second,V8.IL[0],label='IL')
    plot(V9.t/second,V9.IKL[0],label='IKL')
    plot(V10.t/second,V10.IH[0],label='IH')
    plot(V11.t/second,V11.ITLT[0],label='ITLT')
    xlabel('Time (s)')
    legend()
    
    clear_cache('cython') 
    
#dCa_iTLT/dt = ((-10*ITLT)/(2*96485.3* coulomb * mole ** -1 * meter)- (Ca_iTLT-0.00024*mmolar)/5/ms ) : mole * meter**-3
    
    
#    ITLT = gTLT_HTC * mTLT**2 * hTLT * (V-ETLT) : amp * meter ** -2
#        dmTLT/dt = (mTLTinf-mTLT)/taumTLT : 1
#        dhTLT/dt = (hTLTinf-hTLT)/tauhTLT : 1 
#        mTLTinf = 1/(1+exp(-(V+59*mV)/6.2/mV)) : 1
#        taumTLT= 1*ms: second    
#        hTLTinf = 1/(1+exp((V+83*mV)/4/mV)) : 1
#        tauhTLT=  (30.8+ (211.4 + exp((V+115.2*mV)/5/mV))/(1+exp((V+86*mV)/3.2/mV)))/3.7372*ms : second
#        dCa_iTLT/dt = ((-10*ITLT)/(2*96485.3* coulomb * mole ** -1 * meter)- (Ca_iTLT-0.00024*mmolar)/5/ms ) : mole * meter**-3
#        ETLT=R*T/(z*F)*log(2*mmolar/Ca_iTLT)/log(2) : volt
    
#    ITHT = gTHT_HTC * mTHT**2 * hTHT * (V-ETHT) : amp * meter ** -2
#        dmTHT/dt = (mTHTinf-mTHT)/taumTHT : 1
#        dhTHT/dt = (hTHTinf-hTHT)/tauhTHT : 1 
#        mTHTinf = 1/(1+exp(-(V+40.1*mV)/3.5/mV)) : 1
#        taumTHT= 1*ms: second    
#        hTHTinf = 1/(1+exp((V+62.2*mV)/5.5/mV)) : 1
#        tauhTHT=  0.1483*ms * exp(-0.09398 * V/mV) + 5.284*ms * exp(0.008855 * V/mV) : second
#        dCa_iTHT/dt = ((-10*ITHT)/(2*96485.3* coulomb * mole ** -1 * meter)- (Ca_iTHT-0.00024*mmolar)/5/ms ) : mole * meter**-3
#        ETHT=R*T/(z*F)*log(2*mmolar/Ca_iTHT)/log(2) : volt    
#        
    
    
#SK current : Quantifying the Neural Elements Activated and Inhibited by Globus Pallidus Deep Brain Stimulation
#Matthew D. Johnson and Cameron C. McIntyre  2008
    
    #    J = Jmax * int(t>100*msecond) * int(t<800*msecond) : amp * meter ** -2

#    ISK = gSK_HTC * cSK * (V-EK_HTC) : amp * meter ** -2
#    dcSK/dt = (cSKinf-cSK)/taucSK : 1
#    cSKinf = 0.81 / (1 + exp(-(log(Ca_iTLT/mmolar)+0.3)/0.46)): 1
#    taucSK = 6.1*msecond : second
    
    
#    Iapp = randn() *gapp : amp * meter ** -2 (constant over dt)


#    ITLT = gTLT_HTC * mTLT**2 * hTLT * (V-ETLT) : amp * meter ** -2
#        dmTLT/dt = (mTLTinf-mTLT)/taumTLT : 1
#        dhTLT/dt = (hTLTinf-hTLT)/tauhTLT : 1 
#        mTLTinf = 1/(1+exp(-(V+59*mV)/6.2/mV)) : 1
#        taumTLT= 1*ms: second    
#        hTLTinf = 1/(1+exp((V+83*mV)/4/mV)) : 1
#        tauhTLT=  (30.8+ (211.4 + exp((V+113.2*mV)/5/mV))/(1+exp((V+84*mV)/3.2/mV)))/3.7372*ms : second
#        dCa_iTLT/dt = ((-10000000*ITLT)/(2*96485.3* coulomb * mole ** -1 * meter)*int(ITLT<0*amp * meter ** -2)- (Ca_iTLT-0.00024*mmolar)/5/ms ) : mole * meter**-3
#        ETLT=R*T/(z*F)*log(2*mmolar/Ca_iTLT) : volt
#        
#    ITHT = gTHT_HTC * mTHT**2 * hTHT * (V-ETHT) : amp * meter ** -2
#        dmTHT/dt = (mTHTinf-mTHT)/taumTHT : 1
#        dhTHT/dt = (hTHTinf-hTHT)/tauhTHT : 1 
#        mTHTinf = 1/(1+exp(-(V+40.1*mV)/3.5/mV)) : 1
#        taumTHT= 1*ms: second    
#        hTHTinf = 1/(1+exp((V+62.2*mV)/5.5/mV)) : 1
#        tauhTHT=  0.1483*ms * exp(-0.09398 * V/mV) + 5.284*ms * exp(0.008855 * V/mV) : second
#        dCa_iTHT/dt = ((-10000000*ITHT)/(2*96485.3* coulomb * mole ** -1 * meter)*int(ITHT<0*amp * meter ** -2)- (Ca_iTHT-0.00024*mmolar)/5/ms ) : mole * meter**-3
#        ETHT=R*T/(z*F)*log(2*mmolar/Ca_iTHT) : volt
    
#dCa_iTHT/dt = ((-10000000*ITHT)/(2*96485.3* coulomb * mole ** -1 * meter)- (Ca_iTHT-0.00024*mmolar)/5/ms ) : mole * meter**-3
#mAHPinf= 48*(Ca_iTHT/(mmolar))**2 / (48*(Ca_iTHT/(mmolar))**2 + 0.09) : 1
    
    
#        ITHT = gTHT_HTC * mTHT**2 * hTHT * (V-ETHT) : amp * meter ** -2
#        dmTHT/dt = (mTHTinf-mTHT)/taumTHT : 1
#        dhTHT/dt = (hTHTinf-hTHT)/tauhTHT : 1 
#        mTHTinf = 1/(1+exp(-(V+40.1*mV)/3.5/mV)) : 1
#        taumTHT= 1*ms: second    
#        hTHTinf = 1/(1+exp((V+62.2*mV)/5.5/mV)) : 1
#        tauhTHT=  0.1483*ms * exp(-0.09398 * V/mV) + 5.284*ms * exp(0.008855 * V/mV) : second
#        
    
    
#        ITHT = gTHT_HTC * mTHT**2 * (V-ETHT) : amp * meter ** -2
#        dmTHT/dt = (mTHTinf-mTHT)/taumTHT : 1
#        mTHTinf = alphamTHT/(alphamTHT+betamTHT) : 1
#        taumTHT= 1/(alphamTHT+betamTHT): second    
#        alphamTHT= 1.6/ms/(1+exp(-0.072*(V-65*mV)/mV)) : Hz
#        betamTHT= 0.0002/msecond * (V-51.1*mV)/mV / (exp((V-51.1*mV)/5/mV)-1): Hz
#        dCa_iTHT/dt = (-24*ITHT)*mmolar/msecond/(mamp * cm **-2)- Ca_iTHT/100/ms: mole * meter**-3
#        ETHT=120*mV : volt   
    
    
#    V2=V+25*mV : volt
#    ITHT = gTHT_HTC * mTHT**2 * (V-ETHT) : amp * meter ** -2
#        dmTHT/dt = (mTHTinf-mTHT)/taumTHT : 1
#        mTHTinf = alphamTHT/(alphamTHT+betamTHT) : 1
#        taumTHT= 1*msecond/(alphamTHT+betamTHT): second    
#        alphamTHT= 1.6/(1+exp(-0.072*(V2-65*mV)/mV)) : 1
#        betamTHT= 0.02 * (V2-51.1*mV)/mvolt / (exp((V2-51.1*mV)/5/mV)-1): 1
#        dCa_iTHT/dt = (-24*ITHT)*mmolar/msecond/(mamp * cm **-2)- Ca_iTHT/100/ms: mole * meter**-3
#        ETHT=120*mvolt : volt
    
#    *mmolar *umolar ?