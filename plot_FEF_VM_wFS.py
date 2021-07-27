# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 14:40:22 2020

@author: ameli
"""

from brian2 import *
import matplotlib.gridspec as gridspec
import scipy.signal as signal

def read_raster_times_and_indexes(file_t,file_i):
    all_times=[]
    raster_t=open(file_t,'r')
    for line in raster_t:
        time_str=line.split(',')[:-1]
    for line in time_str:
        if line[-2]=='m':
            time=float(line[:-2])*msecond
        elif line[-2]=='u':
            time=float(line[:-2])*usecond
        else :
            time=float(line[:-1])*second
        all_times.append(time)
    raster_t.close()
        
    all_i=[]
    raster_i=open(file_i,'r')
    for line in raster_i:
        all_i=line.split(',')[:-1]
    all_i=[int(i) for i in all_i]
    raster_i.close()
    return array(all_times),array(all_i) 

def spike_train_from_raster(ras_t,runtime):
    record_dt=1/512*second
    L=int(runtime*second/record_dt)
    tvec=arange(0,L*second+record_dt,record_dt)
    spiketrain=zeros(shape(tvec))
    for i in range(len(ras_t)):
        tdist=abs(tvec-ras_t[i]*second)
        spiketrain[tdist.argmin()]+=1
    return spiketrain, tvec

def get_one_raster(group_name):
    t_name=group_name+"_t.txt"
    i_name=group_name+"_i.txt"
    ras_t,ras_i=read_raster_times_and_indexes(t_name,i_name)
    ras=(ras_t,ras_i)
    return ras

def get_raster_Jgg(J,gFS,gVS):
    base="sims/FEF_VM_wFS_J"+str(J)+"_gFS"+str(gFS)+"_gVS"+str(gVS)+"/raster_"
    all_names=["RS","FS","SI","VIP"]
    all_N=[20,20,20,20]
    
    all_raster=[0]*len(all_names)
    
    for i in range(len(all_names)):
        r=get_one_raster(base+all_names[i])
        all_raster[i]=r
    return all_raster
    
J=5

runtime=2*second

gFSSI=[]
gVIPSI=[]
# J=[]
for step in range(11):
    gFSSI.append(1+step*3/10)
    gVIPSI.append(round((1-step/10)*10)/10)
#    J.append(-10+3*step)
    
params=transpose([tile(gFSSI, len(gVIPSI)), repeat(gVIPSI, len(gFSSI))])
# params=transpose([tile(J, len(params)), repeat(params, len(J))])

record_dt=1/512*second
t=int(0.3*second/record_dt) #t_debut
L=int(2*second/record_dt)

fig = figure(figsize=(12,12))
    
outer = gridspec.GridSpec(11, 11)

for pset in range(len(params)):
    
    all_raster=get_raster_Jgg(J,params[pset][0],params[pset][1])
    
    spiketrain, tvec=spike_train_from_raster(all_raster[0][0], 2)
    
    record_dt=1/100000*second#1/512*second
    t=int(0.3*second/record_dt) #t_debut
    L=int(2*second/record_dt)
    fs = 1/record_dt
    
    def flipEnds(mat, end_length):
        beginning = mat[0:end_length, :]
        ending = mat[-end_length-1:-1, :]
        flipped = vstack((flipud(beginning), mat, flipud(ending)))
        return flipped
    
    end_length = 5000
    LFPflip = flipEnds(LFP_I_RS[:, None], end_length)#transpose(atleast_2d(LFP_V_RS)), end_length)
    
    def pctMean(mat, ax):
        diagMean = diag(nanmean(mat, axis=ax))
        matMean = ones(shape(mat))
        if ax == 0:
            matMean = matmul(ones(shape(mat)), diagMean)
        else:
            matMean = matmul(diagMean, ones(shape(mat)))
        normed = (mat - matMean)/matMean
        return normed
    
    f, t, Sxx = signal.spectrogram(transpose(LFPflip), fs, nperseg=25000, noverlap=20000)
    pctSxx = pctMean(Sxx, 1)
    
    inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[floor(pset/11),rem(pset,11)])

    fig.add_subplot(outer[floor(pset/11),rem(pset,11)])
    ax=Subplot(fig, inner[0])
    plot(all_raster[0][0],all_raster[0][1],'k.',label="RS")
            
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    
    xlim(0,runtime/second)
#    legend(loc='upper left')
#    xlabel('Time (s)')
#    ylabel('Neuron index')
    ylim(-1,21)
    
    tight_layout()
    
    figure()
    #f, t, Sxx = signal.spectrogram(LFP_LIP, 100000*Hz,nperseg=30000,noverlap=25000)
    pcolormesh(t, f, pctSxx)#, cmap='RdBu')#, shading='gouraud')
    ylabel('Frequency [Hz]')
    xlabel('Time [sec]')
    ylim(0, 75)
    
    #freq = logspace(0, 2, 50)*Hz
    freq = linspace(1/second, 100*Hz, 100)
    #widths = logspace(log10(3),log10(30),50)*fs/(2*freq*pi)
    widths = linspace(3, 30, 100)*fs/(2*freq*pi)
    
    CWT = signal.cwt(squeeze(LFPflip), signal.morlet2, widths)
    CWT = CWT[:, end_length:-end_length]
    CWTpct = pctMean(CWT, 1)
            
    CWT = signal.cwt(spiketrain, signal.morlet2, widths, w=6)
    
    ax=Subplot(fig, inner[1])
    pcolormesh(tvec, freq, CWT)#, shading='gouraud')
    ylabel('Frequency [Hz]')
    xlabel('Time [sec]')
    ylim(0,50)

savefig("sims/FEF_VM_wFS_J"+str(J)+'_gFS'+str(gFSSI[gf])+'_raster.png')


# all_names=["RS","FS","SI","VIP"]
# all_colors=['r.','b.','g.','k.']
# all_offsets=[60, 40, 20, 0]
# all_N=[20,20,20,20]

# for gf in range(len(gFSSI)):
    
#     figure(figsize=(12,12))
    
#     for gv in range(len(gVIPSI)):
        
#         subplot(3,4,gv+1)
        
#         all_raster=get_raster_Jgg(J,gFSSI[gf],gVIPSI[gv])    
        
#         for r in range(len(all_names)):
#             plot(all_raster[r][0],all_raster[r][1]+all_offsets[r],all_colors[r],label=all_names[r])
            
#         this_ax=gca()
#         this_ax.get_xaxis().set_visible(False)
#         this_ax.get_yaxis().set_visible(False)
        
#         xlim(0,runtime/second)
#     #    legend(loc='upper left')
#     #    xlabel('Time (s)')
#     #    ylabel('Neuron index')
#         ylim(-1,81)
        
#         tight_layout()
    
#     savefig("sims/FEF_VM_wFS_J"+str(J)+'_gFS'+str(gFSSI[gf])+'_raster.png')
    
# for gf in range(len(gFSSI)):
    
#     figure(figsize=(12,12))
    
#     for gv in range(len(gVIPSI)):
        
#         subplot(3,4,gv+1)
        
#         all_raster=get_raster_Jgg(J,gFSSI[gf],gVIPSI[gv])
#         f,Spectrum_LFP_V_RS=signal.periodogram(LFP_V_RS, 100000,'flattop', scaling='spectrum')    
        
#         for r in range(len(all_names)):
#             plot(all_raster[r][0],all_raster[r][1]+all_offsets[r],all_colors[r],label=all_names[r])
            
#         this_ax=gca()
#         this_ax.get_xaxis().set_visible(False)
#         this_ax.get_yaxis().set_visible(False)
        
#         xlim(0,runtime/second)
#     #    legend(loc='upper left')
#     #    xlabel('Time (s)')
#     #    ylabel('Neuron index')
#         ylim(-1,81)
        
#         tight_layout()
    
#     savefig("sims/FEF_VM_wFS_J"+str(J)+'_gFS'+str(gFSSI[gf])+'_raster.png')
    
