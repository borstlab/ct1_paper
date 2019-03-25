# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 08:24:10 2019

@author: aborst
"""

import numpy as np
import matplotlib.pyplot as plt

# swc file format

# 0 1 2 3 4 5 6
# n T x y z R P 
# n is an integer label that identifies the current point and increments by one from one line to the next. 
# T is an integer representing the type of neuronal segment, such as soma, axon, apical dendrite, etc.  
# 0 = undefined
# 1 = soma
# 2 = axon
# 3 = dendrite
# 4 = apical dendrite
# 5 = fork point
# 6 = end point
# 7 = custom
# x, y, z gives the cartesian coordinates of each node. 
# R is the radius at that node. 
# P indicates the parent (the integer label) of the current point or -1 to indicate an origin (soma). 

dirname='CT1_swc_deposit/'
fragments=[168,183,204,220,211,83561775,83622464,84618194,86749209,87275759,88675974]
colors=['black','blue','red','green','magenta','brown','yellow','green','black','blue','red']
noffragments=11
nofcomps=np.zeros(noffragments)

for i in range(noffragments):
    
    CT1=np.loadtxt(dirname+np.str(fragments[i])+'.swc')
    interim=CT1.shape
    nofcomps[i]=interim[0]
    print 'nofcomps', nofcomps[i]
    plt.scatter(CT1[:,2],CT1[:,3],color=colors[i], label=str(i)+': '+str(nofcomps[i]))
    
totalnofcomps=int(np.sum(nofcomps))

print 'total # of comps: ', totalnofcomps
    
BigCT1=np.zeros((totalnofcomps,7))
BigCT1[:,:]=BigCT1.astype(int)

connectpoint=np.array([451,858,1256,1642,1911,2297,2724,2991,3265,3518])
anchorpoint=np.array([39,45,27,80,24,10,32,38,25,70])

def build_BigCT1(n=noffragments, plotit=0):
    
    offset=0

    for i in range(n):
        
        interfile=np.loadtxt(dirname+np.str(fragments[i])+'.swc')
        BigCT1[offset:offset+nofcomps[i]]=interfile.astype(int)
        BigCT1[offset:offset+nofcomps[i],0]=BigCT1[offset:offset+nofcomps[i],0]+offset
        BigCT1[offset:offset+nofcomps[i],6]=BigCT1[offset:offset+nofcomps[i],6]+offset
        BigCT1[offset,6]=-1
        print 'offset=', offset
        offset=offset+nofcomps[i]
    
    for i in range(n-1):
        BigCT1[connectpoint[i],6]=anchorpoint[i]

    if plotit==1:

        colors=np.arange(0,totalnofcomps)
        plt.figure()
        plt.scatter(BigCT1[:,2],BigCT1[:,3],c=colors,s=50)
        
def reorder(oldtree,x,y):
    
    print 'x=',x, 'y=',y
    
    n=oldtree.shape
    n=n[0]
    
    newtree=1.0*oldtree

    # invert the segment between x and y, but keep tree[0] and tree[6] indices

    oldseq=1.0*oldtree[x:y+1]
    newtree[x:y+1]=oldseq[::-1] 
    
    # look for connection points beyond y and replace tree[6] with new indices 
    
    for i in range(y+1,n,1):
        for j in range(x,y,1):
            if newtree[i,6]==oldtree[j,0]:
                oldtree[i,6]=newtree[j,0] 
                #print 'in line ',i,'loop ', j, ':  oldindex', oldtree[j,0], 'replaced by ',newtree[j,0]
                
    # copy oldtree[y:n,6] into newtree[y:n,6]
                
    newtree[y:n,6]=oldtree[y:n,6]
                
    newtree[x,0]=x+1
        
    for i in range(x+1,y+1,1):
        newtree[i,0]=i+1
        newtree[i,6]=i
                
    return newtree
    
def reorder_BigCT1():
    
    stitchpoint=np.where(BigCT1[:,6]==-1)
    stitchpoint=stitchpoint[0]
    
    for i in range(noffragments-1):       
        BigCT1[:,:]=reorder(BigCT1,stitchpoint[i+1],connectpoint[i])
    
def draw_tree(tree, test=0):
    
    n=tree.shape
    nofcomps=n[0]
    plt.scatter(tree[:,2],tree[:,3])
    
    for i in range(1,nofcomps,1):
        if tree[i,6]!=-1:
            a=tree[i,0]-1
            b=tree[i,6]-1
            x1=tree[a,2]
            x2=tree[b,2]
            y1=tree[a,3]
            y2=tree[b,3]
            plt.plot([x1,x2],[y1,y2],color='black')
            if abs(x1-x2)>300:
                print 'index', i, tree[i,0], tree[i,6]
                plt.plot([x1,x2],[y1,y2],color='red')
                plt.scatter(tree[i,2],tree[i,3],color='red')
                
            if test==1:
                print i,a,b,x1,x2

build_BigCT1()
BigCT1[:,1]=0
BigCT1=BigCT1.astype(int)
print np.where(BigCT1[:,6]==-1)
reorder_BigCT1()
print np.where(BigCT1[:,6]==-1)
plt.figure()
draw_tree(BigCT1)
np.savetxt(dirname+'BigCT1.swc',BigCT1.astype(int),fmt='%-5d')
BigCT1forNeuron=1.0*BigCT1
BigCT1forNeuron[:,2:6]=0.01*BigCT1forNeuron[:,2:6]
np.savetxt(dirname+'BigCT1forNeuron.swc',BigCT1forNeuron)

