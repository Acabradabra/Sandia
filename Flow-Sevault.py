#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from numpy import *
import Utilities as util
(plt,mtp)=util.Plot0()

Re0=15e3
D=5e-3
mu =1.15e-5 #9e-6
rho=3.3e-1  #8e-2
nu=mu/rho
V=nu*Re0/D ; print('=> V = %g m/s'%(V))

H2=array([ 0.55, 0.45, 0.37]) #, 0.55, 0.55,0.55 ])
Re=array([   15,   15,   15]) #,   12,   15,  18 ])
Vm=array([ 98.2, 84.5, 75.8]) #, 78.6, 98.2,117.8])
Vc=array([0.778,0.755,0.739]) #,0.622,0.778,0.933])

Rv=Vm/Vc
P_r=polyfit(H2,Rv,2)
P_v=polyfit(H2,Vm,2)
H2_2=arange(0.37,1.001,0.01)
Rv_2=polyval(P_r,H2_2)
Vm_2=polyval(P_v,H2_2)

print(f'=> 100 %    Vm : {Vm_2[-1]:.3f} [m/s]  ,  Vc : {Vm_2[-1]/Rv_2[-1]:.3f} [m/s]')

R_phi=(4-3*H2)/(0.32*2) ; print(R_phi)

fig_r,ax_r=plt.subplots(figsize=(13,10))
ax_r.plot(H2  ,Rv  ,'.k')
ax_r.plot(H2_2,Rv_2,'-k')
# ax_r.plot(H2  ,1/R_phi,':k')
util.SaveFig(fig_r,'PLOT/Sevault_Rv.png')

fig_v,ax_v=plt.subplots(figsize=(13,10))
ax_v.plot(H2  ,Vm  ,'.k')
ax_v.plot(H2_2,Vm_2,'-k')
ax_v.plot([0.55],[V],'og')
util.SaveFig(fig_v,'PLOT/Sevault_Vm.png')