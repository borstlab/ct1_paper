# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 14:24:40 2015

@author: aborst
"""

import numpy as np
import matplotlib.pyplot as plt
import blindschleiche as bs

# --- 600 ms per tstep ------------

init=0

Rm=8000.0 # Ohm cm^2
Ra=0400.0 # Ohm cm
Cm=0.6    # microFarad/(cm**2)
tau=Rm*Cm # microsecond
deltat=0.001

maxtime=50
injtime=30      
injcomp=0
curramp=10*(10**(-12))  # picoAmpere

swcfile1='CT1_swc_deposit/211.swc'
swcfile2='CT1_swc_deposit/BigCT1.swc'
swcfile3='othercell_swc_deposit/cable.swc'
swcfile4='othercell_swc_deposit/simpletree.swc'

# ----- artificial test neurons --------------

def build_cable(nofcomps,complength,compdiam):
    
    cable=np.zeros((nofcomps,7))
    
    for i in range(nofcomps):
        
        cable[i,0]=i+1
        cable[i,2]=i*complength
        cable[i,5]=compdiam/2.0
        cable[i,6]=i
        
    cable[0,6]=-1
    
    return cable
    
def build_tree(nofcomps,complength,compradius):
    
    tree=np.zeros((nofcomps,7))
    
    for i in range(nofcomps/2):
        
        tree[i,0]=i+1
        tree[i,2]=i*complength
        tree[i,5]=compradius
        tree[i,6]=i
        
    for i in range(nofcomps/2,nofcomps):
        
        tree[i,0]=i+1
        tree[i,2]=nofcomps/4*complength
        tree[i,3]=(i+2-nofcomps/2)*complength
        tree[i,5]=compradius*0.5
        tree[i,6]=i
        
    tree[nofcomps/2,6]=nofcomps/4+1
    tree[0,6]=-1
    
    return tree
    
# --------- load BigCT1 swc file ----------------------

if init==1:
    
    mycell=np.loadtxt(swcfile2)
    mycell[:,2:6]=0.01*mycell[:,2:6]
    
# --------calculate Adjancy and Conductance Matrices ---

def calc_AdjM():
    
    interim=mycell.shape
    nofcomps=interim[0]
    
    IdenM=np.identity(nofcomps)
    AdjM=np.zeros((nofcomps,nofcomps))
    
    for i in range(nofcomps):
        if mycell[i,6]!=-1:
            AdjM[i,mycell[i,6]-1]=1
            
    AdjM[:,:]=AdjM[:,:]+np.transpose(AdjM[:,:])+IdenM  
    
    return AdjM

def calc_Conductance_M():
    
    # --- CondM_tdep, CondM_steady --------
    
    interim=mycell.shape
    nofcomps=interim[0]
    
    AdjM=calc_AdjM()
    
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
            
    complength[0]=complength[1]

    # --------------------------------------
    
    CondM_tdep=0.0*AdjM
    CondM_steady=0.0*AdjM
        
    # set conductances
      
    absRm=(10**8)*Rm/(compdiam*np.pi*complength)
    gleak=1.0/absRm
    
    absCm=(10**-8)*Cm*(compdiam*np.pi*complength)
    memcap=absCm*(10**-6)
    
    absRa=(10**4)*Ra/(compdiam**2.0/4.0*np.pi)*complength
    glength=1.0/absRa
            
    # set up matrix
    
    for i in range(nofcomps):
        for j in range(nofcomps):
            if i!=j and AdjM[i,j]==1:               
                CondM_tdep[i,j]   = -0.5*(glength[i]+glength[j])
                CondM_steady[i,j] = -0.5*(glength[i]+glength[j])
        CondM_tdep[i,i]=gleak[i]+memcap[i]/deltat-np.sum(CondM_tdep[i])
        CondM_steady[i,i]=gleak[i]-np.sum(CondM_steady[i])
        
    return CondM_tdep,CondM_steady,memcap
    
if init==1:
    
    CondM_tdep,CondM_steady,memcap=calc_Conductance_M()
    
# ---- mycell, CondM_tdep, CondM_steady, and memcap are all global variables ---
    
# ------------ end of all initialization -------------------
                
def calc_Vm_tdep(injcomp,curramp):
    
    interim=CondM_tdep.shape
    nofcomps=interim[0]
    
    Vm = np.zeros((nofcomps,maxtime))
    currinj= np.zeros((nofcomps,maxtime))   
    currinj[injcomp,11:(11+injtime)]=curramp
    
    print 'timestep ',
    
    for t in range(1, maxtime):
        print t,
        rightsideofeq=Vm[:,t-1]*memcap/deltat+currinj[:,t]
        Vm[:,t] = np.linalg.solve(CondM_tdep,rightsideofeq)
    
    Vm[:,:] = 1000.0*Vm # mV
    
    return Vm
    
def calc_Vm_steady(injcomp,curramp):
    
    interim=CondM_steady.shape
    nofcomps=interim[0]

    currinj=np.zeros(nofcomps)  
    currinj[injcomp]=curramp

    Vm = np.linalg.solve(CondM_steady,currinj)
    
    Vm = 1000.0*Vm # mV
    
    return Vm
    
# end of claculations -------
    
def showM(m,x,y):
    
    plt.matshow(m[x:y,x:y],cmap='gray')   
    
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
    
    global Rm, Ra, CondM_steady
    
    Vm_Ratio=np.zeros((10,10))

    injcomps=[2972,1689,434,1892,3280,1232,282,2275,3492,1089,2699]
    
    for i in range(10):
        Ra=(i+1)*50.0
        for j in range(10):
            Rm=(j+1)*1000.0
            print Ra, Rm
            CondM_tdep,CondM_steady,memcap=calc_Conductance_M()
            Vm=calc_Vm_steady(injcomps[10],curramp)
            Vm_Ratio[j,i]=np.mean(Vm[injcomps[0:10]])/Vm[injcomps[10]]*100.0
            
    np.save('Vm_Ratio', Vm_Ratio)
    
def plot_paramscan(interpol=1):
    
    Vm_Ratio=np.load('Vm_Ratio.npy')
    
    plt.figure(figsize=(10,7))
    
    if interpol==1:       
        plt.imshow(Vm_Ratio, vmin=0, vmax=50, origin='lower')
    if interpol==0:       
        plt.imshow(Vm_Ratio, vmin=0, vmax=50, origin='lower', interpolation='None')
    plt.colorbar()    
    plt.xticks(np.arange(10),(np.arange(10)+1)*50)
    plt.xlabel('Ra [Ohm cm]')
    plt.yticks(np.arange(10),(np.arange(10)+1)*1000)
    plt.ylabel('Rm [Ohm cm^2]')
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
    
    interim=CondM_steady.shape
    nofcomps=interim[0]
    gsyn=np.zeros(nofcomps)
    
    for i in range(3):  
        gsyn[syncomps[i]]=synamps[i]*(10**(-12))
    
    synCondM=CondM_steady+np.identity(nofcomps)*gsyn

    Vm = np.linalg.solve(synCondM,0.05*gsyn)
    
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
    mycell=np.loadtxt(swcfile2)
    dirname='CT1_swc_deposit/'
    
    for i in range(12):
        
        MyVm=VmDistribution[i]
        MyVm=MyVm/np.max(MyVm)
        MyVm=MyVm*240+20
        MyVm=MyVm.astype(int)
        
        mycell[:,1]=MyVm
        
        filename=dirname+'colorBigCT1ci'+str(i)+'.swc'
        
        np.savetxt(filename,mycell.astype(int),fmt='%-5d')
    
    inpR=np.zeros(12)
    injcomps=[0,2972,1689,434,1892,3280,1232,282,2275,3492,1089,2699]
    
    print 'inj comp','   inpR [GOhm]'   
    
    for i in range(12):
        inpR[i]=np.max(VmDistribution[i])*10**(-3)/curramp
        print injcomps[i],'      ', inpR[i]*10**(-9)
    
    np.save('inpR',inpR)
         
    
def showlamda(mycell,Vm_space):
    
    diam=mycell[0,5]*2.0
    lamda=np.sqrt(Rm/Ra*diam/4.*(10**(-4)))*10**4
    print 'lamda [micrometer] = ', lamda
    plt.plot(mycell[:,2],Vm_space)
    l=1000.0
    x=np.arange(l)
    Vm_calc=Vm_space[0]*np.cosh((l-x)/lamda)/np.cosh(l/lamda)
    plt.plot(x,Vm_calc)

    
    
    
        
    
    
    




