#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%=================================================================================
#                     Modules
#===================================================================================
import Utilities as util
(Sysa,NSysa,Arg)=util.Parseur(['Save','Visu','Temp','Data','Over','Vel','Only','Compa','All'],0,'Arg : Case (Jaravel,Garnier,Sevault)')
(                             [ SAVE , VISU , TEMP , DATA , OVER , VEL , ONLY , COMPA , ALL ])=Arg
if ALL : SAVE,VISU,TEMP,OVER,DATA=True,True,True,True,True
if VEL : VISU,SAVE=True,True

from numpy import *
import os
import sys
import csv
import time
import Fluent as fl

t0=time.time()
(plt,mtp)=util.Plot0()

#%%=================================================================================
#                     Parameters
#===================================================================================

dir0='/mnt/scratch/ZEUS/FLUENT/Sandia-Sevault/RUN-99-Scaling/SCALE/'

#%%=================================================================================
#                     Process
#===================================================================================

Scales=[ s for s in os.listdir(dir0) if s[-4:]=='.out' ]# ; print(Scales)
Ns=len(Scales)

It    =zeros(Ns)
Cpu   =zeros(Ns)
T_iter=zeros(Ns)
T_tot =zeros(Ns)
Nn    =zeros(Ns) 
for n,s in enumerate(Scales) :
    util.Section(s,0,3,'b')
    Lines=[ l for l in os.popen('tail -n1000 '+dir0+s).read().splitlines() if l ]
    if   '1n' in s : Nn[n]=1
    elif '2n' in s : Nn[n]=2
    else : sys.exit('=> Wrong number of node')
    for l in Lines :
        l2=str(l) #; print(l2)
        if 'Performance Timer for'                 in l2 : L0=l2.split(' ') ; it=int(L0[3]) ; cpu=int(L0[6])
        if 'Average wall-clock time per iteration' in l2 : t_iter=float(l2.split(':')[-1].strip()[:-3])
        if 'Total wall-clock time'                 in l2 : t_tot =float(l2.split(':')[-1].strip()[:-3])
    print( '=> It : %.0f  ,  cpu : %.0f  ,  t/iter : %.3f  ,  t_tot : %.3f'%(it,cpu,t_iter,t_tot) )
    It    [n]=it
    Cpu   [n]=cpu
    T_iter[n]=t_iter
    T_tot [n]=t_tot

#%%=================================================================================
#                     Ploting
#===================================================================================
N1=(Nn==1)
N2=(Nn==2)

fig,ax=plt.subplots(figsize=(8,6))
ax.set_xlabel('N tasks'               ,fontsize=30)
ax.set_ylabel('Time per iteration [s]',fontsize=30)
ax.plot(2*[128],[3,9],':k')
ax.plot(2*[256],[3,9],':k')
ax.plot(2*[512],[3,9],':k')
ax.plot(Cpu[N1],T_iter[N1],'ob')
ax.plot(Cpu[N2],T_iter[N2],'or')
util.SaveFig(fig,dir0+'Scaling.pdf')
