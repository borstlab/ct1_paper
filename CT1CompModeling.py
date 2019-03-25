# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 14:24:40 2015

@author: aborst
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy import sparse
from scipy.sparse.linalg import spsolve

# 2.5 ms computing time per tstep

init=1

mychoice=0

Rm=8000.0 # Ohm cm^2
Ra=0400.0 # Ohm cm
Cm=0.6    # microFarad/(cm**2)
tau=Rm*Cm # microsecond
deltat=0.001

maxtime=50
injtime=30      
injcomp=0
curramp=10*(10**(-12))  # 10 picoAmpere

swcfile='CT1_swc_deposit/BigCT1.swc'
 
# --------- load BigCT1 swc file ----------------------

if init==1:
    
    mycell=np.loadtxt(swcfile)
    mycell[:,2:6]=0.01*mycell[:,2:6]
    interim=mycell.shape
    nofcomps=interim[0]
    
# --------calculate Conductance Matrices ---

def calc_Conductance_M():
    
    # --- define matrix M --------
            
    M=np.zeros((nofcomps,nofcomps))
    
    # -------------------------------------
    
    compdiam=mycell[:,5]*2.0    # in micrometer
    complength=np.zeros(nofcomps)
    
    # complength defined backwards
    
    for i in range(1,nofcomps,1):
            
            aind=mycell[i,0]-1
            bind=mycell[i,6]-1
            axyz=mycell[aind,2:5]
            bxyz=mycell[bind,2:5]
            
            complength[i]=np.sqrt(np.sum((axyz-bxyz)**2)) # in micrometer
            
            meandiam=(compdiam[aind]+compdiam[bind])*0.5
            area=meandiam**2.0/4.0*np.pi
            M[bind,aind]=-area/complength[aind]/Ra*10**(-4)
            M[aind,bind]=M[bind,aind]
            
    complength[0]=complength[1]
    
    gleak=(compdiam*np.pi*complength)/(Rm*10**8)
    memcap=(compdiam*np.pi*complength)*Cm*(10**-6)/(10**8)
    
    for i in range(nofcomps):
        M[i,i]=gleak[i]-np.sum(M[i])
    
    M=sparse.csr_matrix(M)
        
    return M,memcap
    
if init==1:
    
    M,memcap=calc_Conductance_M()
    
# ------------ end of all initialization -------------------
                
def calc_Vm_tdep(injcomp,curramp):
    
    Vm = np.zeros((nofcomps,maxtime))
    currinj= np.zeros((nofcomps,maxtime))   
    currinj[injcomp,11:(11+injtime)]=curramp
    
    M_tdep=1.0*M
    M_tdep.setdiag(M.diagonal()+memcap/deltat)
    
    print 'timestep ',
    
    for t in range(1, maxtime):
        print t,
        rightsideofeq=Vm[:,t-1]*memcap/deltat+currinj[:,t]
        Vm[:,t] = spsolve(M_tdep,rightsideofeq)
    
    Vm[:,:] = 1000.0*Vm # mV
    
    return Vm
    
def calc_Vm_steady(injcomp,curramp):

    currinj=np.zeros(nofcomps)  
    currinj[injcomp]=curramp

    Vm = spsolve(M,currinj)
    
    Vm = 1000.0*Vm # mV
    
    return Vm
    
# -------- end of claculations -------  
    
def draw_Vm_tree(tree,Vm,DC,dim1=2,dim2=3):

    plt.figure(figsize=(12,7))    
    
    n=tree.shape
    nofcomps=n[0]
    
    MyVm=Vm-DC*np.min(Vm)
    MyVm=MyVm/np.max(MyVm)
    sfac=10
    
    for i in range(1,nofcomps,1):
        intensity=int(256*MyVm[i])
        if tree[i,6]!=-1:
            a=tree[i,0]-1
            b=tree[i,6]-1
            x1=tree[a,dim1]
            x2=tree[b,dim1]
            y1=tree[a,dim2]
            y2=tree[b,dim2]
            diam=tree[i,5]+tree[i-1,5]
            mylw=int(sfac*diam)
            plt.plot([x1,x2],[y1,y2],linewidth=mylw,color=plt.cm.rainbow(intensity))

def syncurrinj(n=12):
      
    VmDistribution=np.zeros((n,3749))
    injcomps=[0,2972,1689,434,1892,3280,1232,282,2275,3492,1089,2699]
    
    for i in range(n):

        print i,        
        injcomp=injcomps[i]
        Vm=calc_Vm_steady(injcomp,curramp)
        VmDistribution[i]=Vm
        
    np.save('VmDistribution', VmDistribution)
    
def scanparameterspace():
    
    global Rm, Ra, M
    
    Vm_Ratio=np.zeros((10,10))

    injcomps=[2972,1689,434,1892,3280,1232,282,2275,3492,1089,2699]
    
    for i in range(10):
        Ra=(i+1)*50.0
        for j in range(10):
            Rm=(j+1)*1000.0
            M,memcap=calc_Conductance_M()
            Vm=calc_Vm_steady(injcomps[10],curramp)
            Vm_Ratio[j,i]=np.mean(Vm[injcomps[0:10]])/Vm[injcomps[10]]*100.0
            print Ra, Rm, Vm_Ratio[j,i]
            
    np.save('Vm_Ratio', Vm_Ratio)
    
def plot_paramscan(interpol=1,myvmax=70):
    
    Vm_Ratio=np.load('Vm_Ratio.npy')
    
    plt.figure(figsize=(10,7))
    
    if interpol==1:       
        plt.imshow(Vm_Ratio, vmin=0, vmax=myvmax, origin='lower')
    if interpol==0:       
        plt.imshow(Vm_Ratio, vmin=0, vmax=myvmax, origin='lower', interpolation='None')
    plt.colorbar()    
    plt.xticks(np.arange(5)*2+1,(np.arange(5)*2+2)*50)
    plt.xlabel('Ra ['+r'$\Omega$'+'cm]')
    plt.yticks(np.arange(10),(np.arange(10)+1))
    plt.ylabel('Rm [k'+r'$\Omega$'+'cm^2]')
    plt.title('Mean Vm [% of max]', fontsize=18)
    
    plt.contour(Vm_Ratio,[20.0])

def Gauss1D(FWHM,RFsize):
    
    myrange=RFsize/2
    sigma=FWHM/(2.0*np.sqrt(2*np.log(2)))
    x=np.arange(-myrange,(myrange+1),1)*1.0
    z=np.exp(-x**2/(2*(sigma**2)))
    z=z/np.sum(z)
    return z
    
def calc_synresponse(synamps=[0,100,0],syncomps=[3492,1089,2699]):
    
    # gsyn in pico Siemens

    gsyn=np.zeros(nofcomps)
    
    for i in range(3):  
        gsyn[syncomps[i]]=synamps[i]*(10**(-12))
        
    M_syn=1.0*M
    M_syn.setdiag(M.diagonal()+gsyn)
    
    Vm = spsolve(M_syn,0.05*gsyn)
    
    Vm = 1000.0*Vm # mV
    
    return Vm
    
def calc_RF(synfac=1000.0):
    
    inp_RF=np.zeros((3,31))
    out_RF=np.zeros(31)
    
    inp_RF[1]=Gauss1D(7,30)
    inp_RF[0]=np.roll(inp_RF[1],-5)
    inp_RF[2]=np.roll(inp_RF[1],+5)
    
    for i in range(31):
        print i
        synamp=[inp_RF[0,i]*synfac,inp_RF[1,i]*synfac,inp_RF[2,i]*synfac]
        Vm=calc_synresponse(synamps=synamp,syncomps=[3492,1089,2699])
        out_RF[i]=Vm[1089]
        
    for i in range(3):
        plt.plot(inp_RF[i]/np.max(inp_RF[i]))
    plt.plot(out_RF/np.max(out_RF))
    
    RF=np.zeros((4,31))
    RF[0:3]=inp_RF
    RF[3]=out_RF
    
    np.save('RecFields',RF)
 
def create_colorswcfile(n=12):
    
    VmDistribution=np.load('VmDistribution.npy')
    mycell=np.loadtxt(swcfile[0])
    dirname='CT1_swc_deposit/'
    
    for i in range(12):
        
        MyVm=VmDistribution[i]
        MyVm=MyVm/np.max(MyVm)
        MyVm=MyVm*240+20
        MyVm=MyVm.astype(int)
        
        mycell[:,1]=MyVm
        
        filename=dirname+'colorBigCT1ci'+str(i)+'.swc'
        
        np.savetxt(filename,mycell.astype(int),fmt='%-5d')
        
def calc_inputsignals():
    
    myimage=np.zeros((320,320))
    spat_freq=5
    sinewave=np.sin(np.arange(0,320)/320.*2.0*np.pi*spat_freq)
    myimage[:,0::]=sinewave
    myimage=myimage*0.5+0.5

    maxtime=500
    maxvelo=0.5   # x 22.5 deg per second
    tf=np.zeros(maxtime)
    tf[100:maxtime-100]=maxvelo

    movie=np.zeros((320,320,maxtime))         
   
    for t in range(maxtime):
        movie[:,:,t]=np.roll(myimage,int(np.sum(tf[0:t])),axis=1)

    inputsignals=np.zeros((12,12,500))

    for y in range(12):
        ypos=22*(y+1)
        for x in range(12):
            xpos=22*(x+1)
            inputsignals[y,x]=movie[ypos,xpos]

    return inputsignals
        
def calc_motionresponses():  

    maxtime=500
    deltat=0.01    
    synamp=10**(-10)
    
    Vm = np.zeros((nofcomps,maxtime))
    
    syncomps=[2972,1689,434,1892,3280,1232,282,2275,3492,1089,2699]
    
    x=[7,8,9,10,7,8,9,10,8,9,10]
    y=[9,9,9,9,8,8,8,8,7,7,7]
    
    M,memcap=calc_Conductance_M()
    
    gsyn=np.zeros((nofcomps,maxtime))   
    
    input_signals=calc_inputsignals()
    
    for i in range(11):
        
        gsyn[syncomps[i]]=input_signals[y[i],x[i]]*synamp
     
    M_tdep=1.0*M
    
    print 'timestep ',
    
    for t in range(1, maxtime):
        print t,
        M_tdep.setdiag(M.diagonal()+memcap/deltat+gsyn[:,t])
        rightsideofeq=Vm[:,t-1]*memcap/deltat+0.05*gsyn[:,t]
        Vm[:,t] = spsolve(M_tdep,rightsideofeq)
    
    Vm[:,:] = 1000.0*Vm # mV
    
    myt=myt=np.arange(500)*0.01
    
    mycolor=np.array(['red','blue'])
    
    plt.subplot(2,1,1)
    
    for i in range(2):
        plt.plot(myt,gsyn[syncomps[i]]/(10**(-12)),linewidth=2,color=mycolor[i],label=str(i+1))
                   
    plt.legend(loc=2,frameon=False,fontsize=10)
    plt.ylim(-10,110)

    plt.ylabel('Synaptic Conductance [pS]')
    
    plt.subplot(2,1,2)
                
    for i in range(2):
        plt.plot(myt,Vm[syncomps[i]],linewidth=2,color=mycolor[i],label=str(i+1))
                       
    plt.legend(loc=2,frameon=False,fontsize=10)
    plt.xlabel('time [sec]')
    plt.ylabel('Membrane Response [mV]')
    
    plt.ylim(0,20)
    
    return Vm

        
    
    
    


    
    
    
        
    
    
    




