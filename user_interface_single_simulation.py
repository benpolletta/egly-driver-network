#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:27:34 2019

@author: aussel
"""

import matplotlib
import matplotlib.pyplot as plt
import time
#matplotlib.use('agg')
plt.switch_backend('agg')

from brian2 import *

import os
import datetime

#from model_files.FEF_and_LIP_no_SI_rec_noVIP import *
from model_files.FEF_and_LIP_single_simulation import *

from tkinter import *
from tkinter import ttk
from PIL import Image, ImageTk

import ntpath

BG='white'
myfont=('Helvetica',12)


interface=Tk()
interface.minsize(1000, 900)
interface.option_add("*Font", "courier")
interface.option_add("*Background", "white")

modeled_zone=StringVar(interface,'FEF and LIP full network')
theta_phase=StringVar(interface,'mixed')
target_presence=StringVar(interface,'True')
target_presentation_time=DoubleVar(interface,0.6)
runtime=DoubleVar(interface,2)
g_LIP_FEF_v=DoubleVar(interface,0.15)
t_FS=DoubleVar(interface,0.005)
t_SOM=DoubleVar(interface,0.020)


aborted=True


def start():
    global aborted,modeled_zone,theta_phase,target_presence,target_presentation_time,target_presentation_duration,runtime,g_LIP_FEF_v,t_FS,t_SOM,path
    
    modeled_zone=modeled_zone.get()
    theta_phase=theta_phase.get()
    target_presence=target_presence.get()
    target_presentation_time=target_presentation_time.get()*second
    runtime=runtime.get()*second
    g_LIP_FEF_v=g_LIP_FEF_v.get()*msiemens * cm **-2
    t_FS=t_FS.get()*second
    t_SOM=t_SOM.get()*second
    
    
    interface.destroy()
    start_scope()
    
    path=''
    
    
    if os.name == 'nt':
        path=os.path.join(ntpath.dirname(os.path.abspath(__file__)),"results_"+str(datetime.datetime.now()).replace(':','-'))
    else :
        path="./results_"+str(datetime.datetime.now())


    os.mkdir(path)
    
    param_file=open(path+'/parameters.txt','w')
    param_file.write('Modeled region: '+str(modeled_zone)+'\n')
    param_file.write('Theta phase: '+str(theta_phase)+'\n')
    param_file.write('Target presence: '+str(target_presence)+'\n')
    param_file.write('(if Target presence) Target presentation time: '+str(target_presentation_time)+'\n')
    param_file.write('Simulation duration: '+str(runtime)+'\n')
    param_file.write('LIP to FEFv connection strength: '+str(g_LIP_FEF_v)+'\n')
    param_file.write('FS inhibition time constant: '+str(t_FS)+'\n')
    param_file.write('SOM inhibition time constant: '+str(t_SOM))
    param_file.close()    
    
    simu=(target_presentation_time,0,t_SOM,t_FS,theta_phase,g_LIP_FEF_v,target_presence,runtime)
    
    if modeled_zone=='FEF and LIP full network':
        FEF_and_LIP(simu,path,plot_raster=True)
    elif modeled_zone=='LIP alone':
        run_one_LIP_simulation(simu,path,plot_raster=True)
    elif modeled_zone=='FEF visuo_motor module alone':
        sim_FEF_vm_alone(simu,path,plot_raster=True)
    elif modeled_zone=='FEF visual module alone':
        sim_FEF_v_alone(simu,path,plot_raster=True)

    aborted=False
    return


s = ttk.Style()
s.configure('TNotebook.Tab', font=('URW Gothic L','9','bold') )
s.configure('TNotebook', font=('URW Gothic L','9','bold') )
interface.option_add("*Label.Font", "times 8")
interface.option_add("*Font", "times 8") #"Verdana 10 bold"

tab_parent=ttk.Notebook(interface)
tab1=Frame(tab_parent,bg=BG)
tab2=Frame(tab_parent,bg=BG)
tab3=Frame(tab_parent,bg=BG)
tab4=Frame(tab_parent,bg=BG)
tab5=Frame(tab_parent,bg=BG)
tab6=Frame(tab_parent,bg=BG)

tab_parent.add(tab1,text='Zones to model')
tab_parent.add(tab2,text='Network parameters')

tab_parent.pack(expand=1,fill='both')

bquit=Button(interface,text='Quit',command=interface.destroy)
bquit.place(x=10,y=860) 
bstart=Button(interface,text='End setup and start',command=start)
bstart.place(x=80,y=860) 

#bquit=Button(interface,text='Quit',command=interface.destroy)
#bquit.place(x=10,y=560) 
#bstart=Button(interface,text='End setup and start',command=start)
#bstart.place(x=80,y=560) 



### Tab 1 : Zones to model:

question_type=Label(tab1,text='Choose the parts of the network to model:')
question_type.place(x=10,y=10) 

image_modele = Image.open("model_files/model_diagram.png")
photo_modele = ImageTk.PhotoImage(image_modele)

label_modele = Label(tab1,image=photo_modele)
label_modele.image = photo_modele # keep a reference!
label_modele.place(x=50,y=120)

b0 = Radiobutton(tab1, variable=modeled_zone, text='FEF and LIP full network', value='FEF and LIP full network')
b0.place(x=210, y=40) 
b1 = Radiobutton(tab1, variable=modeled_zone, text='LIP alone', value='LIP alone')
b1.place(x=10, y=60)
b2 = Radiobutton(tab1, variable=modeled_zone, text='FEF visuo_motor module alone', value='FEF visuo_motor module alone')
b2.place(x=210, y=60)
b3 = Radiobutton(tab1, variable=modeled_zone, text='FEF visual module alone', value='FEF visual module alone')
b3.place(x=510, y=60)


### Tab 2 : Network parameters:

question_runtime=Label(tab2,text='Choose simulation duration (in s):')
question_runtime.place(x=10,y=10) 
entrytime=Entry(tab2,textvariable=runtime)
entrytime.place(x=250,y=10)

question_input=Label(tab2,text='Choose mdPul input:')
question_input.place(x=10,y=70) 
b0 = Radiobutton(tab2, variable=theta_phase, text='Good theta phase only', value='good')
b0.place(x=200, y=70) 
b1 = Radiobutton(tab2, variable=theta_phase, text='Poor theta phase only', value='bad')
b1.place(x=400, y=70)
b2 = Radiobutton(tab2, variable=theta_phase, text='Good and poor theta phase alternation', value='mixed')
b2.place(x=600, y=70)

question_targetYN=Label(tab2,text='Choose whether to present a target to the network:')
question_targetYN.place(x=10,y=130) 
b0 = Radiobutton(tab2, variable=target_presence, text='Yes', value='True')
b0.place(x=350, y=130) 
b1 = Radiobutton(tab2, variable=target_presence, text='No', value='False')
b1.place(x=450, y=130)

question_targetparam=Label(tab2,text='If so, choose when to present target (in s):')
question_targetparam.place(x=10,y=160) 
entrytime=Entry(tab2,textvariable=target_presentation_time)
entrytime.place(x=300,y=160)

question_gLIPFEF=Label(tab2,text='If the full LIP and FEF network is modeled, choose LIP to FEF visual module synaptic conductance:')
question_gLIPFEF.place(x=10,y=220) 
entrytime=Entry(tab2,textvariable=g_LIP_FEF_v)
entrytime.place(x=650,y=220)

question_tFS=Label(tab2,text='Choose FS cells inhibition decay time constant (in s):')
question_tFS.place(x=10,y=280) 
entrytime=Entry(tab2,textvariable=t_FS)
entrytime.place(x=400,y=280)

question_tSOM=Label(tab2,text='Choose SOM cells inhibition decay time constant (in s):')
question_tSOM.place(x=10,y=310) 
entrytime=Entry(tab2,textvariable=t_SOM)
entrytime.place(x=400,y=310)



interface.mainloop()

#print(aborted)
if not aborted : 
    print('Saving figures')
    os.mkdir(path+'/figures')
    for i in get_fignums():
        current_fig=figure(i)
        current_fig.savefig(path+'/figures/figure'+str(i)+'.png')

    clear_cache('cython')
