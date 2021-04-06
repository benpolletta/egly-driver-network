# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 14:40:22 2020

@author: ameli
"""

from brian2 import *
import matplotlib.animation as animation

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

def get_one_raster_and_pos(group_name,center,N):
    t_name=group_name+"_t.txt"
    i_name=group_name+"_i.txt"
    ras_t,ras_i=read_raster_times_and_indexes(t_name,i_name)
    ras=(ras_t,ras_i)
    R=rand(N)
    theta=rand(N)*2*pi
    pos=array([center[0]+R*cos(theta),center[1]+R*sin(theta)])
    return (ras,pos)

def get_raster_sim_number(nsim):
    base="results_FEF_LIP2/results_"+str(nsim)+"/raster_"
    all_names=["LIP RS","LIP FS","LIP SI","LIP IB","LIP RS gran","LIP FS gran", "LIP SI deep"]
    all_centers=[(2,12),(4.5,15),(7,12),(3,2),(3,7),(6,7),(6,2)]
    all_N=[100,20,20,20,20,20,20]
    all_names+=["FEF RS vm","FEF SI2 vm","FEF SI1 vm","FEF RS v","FEF FS v","FEF VIP v","FEF SI v","FEF RS m"]
    all_centers+=[(15,5),(18,5),(18,2),(13,12),(16,12),(13,9),(16,9),(21,9)]
    all_N+=[20,20,20,20,20,20,20,20]
    
    all_raster=[0]*len(all_names)
    all_pos=[0]*len(all_names)
    
    for i in range(len(all_names)):
        r,p=get_one_raster_and_pos(base+all_names[i],all_centers[i],all_N[i])
        all_raster[i]=r
        all_pos[i]=p
    return all_raster,all_pos


def init_fig(all_raster,all_pos): 
    global all_plots,text_time
    all_centers=[(2,12),(4.5,15),(7,12),(3,2),(3,7),(6,7),(6,2)]
    all_centers+=[(15,5),(18,5),(18,2),(13,12),(16,12),(13,9),(16,9),(21,9)]
    all_colors=['r','b','lime','purple','r','b','lime']
    all_colors+=['r','b','k','lime','r','lime','lime','r']
    theta = linspace(0, 2*pi, 100)
    for i in range(len(all_centers)):
        elem=all_centers[i]
        x1 = elem[0]+1.2*cos(theta)
        x2 = elem[1]+1.2*sin(theta)
        plot(x1, x2,all_colors[i],linewidth=2)
    
    all_plots=[]
    for i in range(len(all_raster)):
#        print('init '+str(i))
        p=scatter(all_pos[i][0,:], all_pos[i][1,:], c=[0 for k in range(len(all_pos[i][0,:]))], marker ='o',cmap='RdGy',vmin=-4, vmax=2)
        all_plots.append(p)
    test_time=text(10,18,'time=0',fontsize=10)
    text(16.5,17,'FEF',fontsize=10)
    text(4.5,17,'LIP',fontsize=10)
    
    plot([10,10],[2,15],'k--')
    xticks([],[])
    yticks([],[])

def update_fig(i):
    global t
    global all_plots,text_time
    all_N=[100,20,20,20,20,20,20]
    all_N+=[20,20,20,20,20,20,20,20]

    text_time.set_text('time = '+str((int((t*1./512-0.3)*1000)))+' ms')
    for j in range(len(all_raster)):
        ras_t,ras_i=all_raster[j]
#        if len(where(ras_t>t*1./512)[0])>0:
        active_tind= array([value for value in where(ras_t>t*1./512)[0] if value in where(ras_t<((t+1)*1./512))[0]]) 
#            active_tind=where(ras_t<((t+1)*1./512))[0][-1]-where(ras_t>t*1./512)[0][0]
#        print('test3a')
#        else :
#            active_tind=array([])
#            print('test3b')
#        print(active_tind)
        if len(active_tind)==0:
            active_i=array([])
        else:
            active_i=ras_i[active_tind]

        color_array=array([0 for k in range(all_N[j])])
        if len(active_i)>0:
            color_array[active_i]=-3

#        all_plots[j]=scatter(all_pos[j][0,:], all_pos[j][1,:], c=color_array, marker ='o')
        all_plots[j].set_array(color_array)
#        print('test4')
    t+=1
    print(t)
    return

record_dt=1/512*second
t=int(0.3*second/record_dt) #t_debut
L=int(2*second/record_dt)
all_raster,all_pos=get_raster_sim_number(250)
seed(1)

fig = figure()
init_fig(all_raster,all_pos)
print('test')
ani = animation.FuncAnimation(fig, update_fig,frames=int(1*second/record_dt),repeat=False, interval=100)
print('test')
ani.save('try_animation.gif',writer='imagemagick', fps=512)
