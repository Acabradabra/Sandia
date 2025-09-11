#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Utilities as util
(Sysa,NSysa,Arg)=util.Parseur(['Save'],0,'Arg : ')
(                             [ SAVE ])=Arg

from numpy import *
from Params import *
import sys
import time
import Fluent as fl

t0=time.time()
(plt,mtp)=util.Plot0()
#===================================================================================
#                     Parameters
#===================================================================================

D=D0
dirr='/mnt/d/Python/SandiaJaravel/REF/'
dird='/mnt/d/FLUENT/Sandia-Jaravel/RUN-01-EBM/DATA/'
dirp=dird+'PLOT/'
Lines=['l0','l1','l2','l3','l75','l15','l30']

#===================================================================================
#                     reading
#===================================================================================
def PlotIm(ax,pic,Ext) :
	im=mtp.image.imread(pic)
	(Ni,Nj,Nk)=im.shape
	ax.imshow( im,extent=Ext )
	ax.set_aspect( (Ext[1]/Ext[3])*(Ni/Nj) )
#-----------------------------------------------------------------------------------
def PlotVar(ax,var,X,M,T,TXT,BD) :
	(Xtxt,txt)=TXT
	if var=='mix' :
		Ic1=T.index('ch4')
		Ic2=T.index('co2')
		Ic3=T.index('co')
		Yc=fl.Yc(M[:,Ic1],M[:,Ic2],M[:,Ic3],Mol_m)
		Var=(Yc-Yc_o)/(Yc_f-Yc_o)
	else :
		I=T.index(var)
		Var=M[:,I]
	ax.plot( X,Var,'r' )
	ax.text( Xtxt[0],Xtxt[1],txt )
	ax.set_xlim((BD[0],BD[1]))
	ax.set_ylim((BD[2],BD[3]))
#===================================================================================
util.MKDIR(dirp)

figA,axA=plt.subplots(ncols=2,        figsize=(20,10))
figT,axT=plt.subplots(ncols=2,nrows=3,figsize=(20,10)) ; figT.suptitle('Radial temperature profiles'     ,fontsize=30)
figU,axU=plt.subplots(ncols=2,nrows=3,figsize=(20,10)) ; figU.suptitle('Radial velocity profiles'        ,fontsize=30)
figZ,axZ=plt.subplots(ncols=2,nrows=3,figsize=(20,10)) ; figZ.suptitle('Radial mixture fraction profiles',fontsize=30)

#=====> Axial Profiles
ym=38
Ym=[1.1,2e3]
Yt=[[0,0.5,1],arange(0,2100,500)]
for n,v in enumerate(['Z','T']) :
	#=====> Read Jaravel
	PlotIm( axA[n],dirr+'0D-{}.png'.format(v),[0,ym,0,Ym[n]] )
	axA[n].set_yticks(Yt[n])
	axA[n].set_xlim((0,ym))
#=====> Read Simulation
(T,M)=fl.ReadSurf(dird+'Data-l0.dat')
([Ix,Iy],[Vx,Vy])=fl.Space(T,M)
#===> Mixture Fraction
Ic1=T.index('ch4')
Ic2=T.index('co2')
Ic3=T.index('co')
Yc_f=fl.Yc(Ych4[0],Yco2[0],Yco[0],Mol_m)
Yc_p=fl.Yc(Ych4[1],Yco2[1],Yco[1],Mol_m)
Yc_o=fl.Yc(Ych4[2],Yco2[2],Yco[2],Mol_m)
Yc  =fl.Yc(M[:,Ic1],M[:,Ic2],M[:,Ic3],Mol_m)
Za=(Yc-Yc_o)/(Yc_f-Yc_o)
axA[0].plot( Vx/D,Za,'r' )
axA[0].set_ylabel('<Z> [-]',fontsize=30)
#===> Temperature
It=T.index('temperature')
axA[1].plot( Vx/D,M[:,It],'r' )
axA[1].set_ylabel('<T> [K]',fontsize=30)

#=====> Radial Profiles
for n,l in enumerate(Lines[1:]) :
	if l[1:]=='75' : txt='7.5'
	else           : txt=l[1:]
	i=n//3
	j=n-(i*3)
	#=====> Read Jaravel
	yu=65
	yz=1.1
	rm=2.5
	if   i==0 : yt=2100 
	elif i==1 : yt=2200
	if j==1 : 
		axT[j,i].set_ylabel('<T> [K]'  ,fontsize=30)
		axU[j,i].set_ylabel('<U> [m/s]',fontsize=30)
		axZ[j,i].set_ylabel('<Z> [-]'  ,fontsize=30)
	PlotIm(axT[j,i],dirr+l[1:]+'D-T.png',[0,rm,0,yt])
	PlotIm(axU[j,i],dirr+l[1:]+'D-U.png',[0,rm,0,yu])
	PlotIm(axZ[j,i],dirr+l[1:]+'D-Z.png',[0,rm,0,yz])
	axT[j,i].set_yticks(arange(0,2100,400))
	#=====> Read simulation
	(T,M)=fl.ReadSurf(dird+'Data-{}.dat'.format(l))
	([Ix,Iy],[Vx,Vy])=fl.Space(T,M)
	PlotVar(axT[j,i],'temperature'       ,Vy/D,M,T,[[2,1800],txt+'D'],[0,rm,0,2250])
	PlotVar(axU[j,i],'velocity-magnitude',Vy/D,M,T,[[2,50  ],txt+'D'],[0,rm,0,62  ])
	PlotVar(axZ[j,i],'mix'               ,Vy/D,M,T,[[2,0.75],txt+'D'],[0,rm,0,yz  ])
axU[0,0].plot(2*[0.5      ],[0,60],':k')
axU[0,0].plot(2*[0.5*D1/D0],[0,60],':k')

if SAVE :
	figA.savefig('PLOT/Profiles-Axis.pdf')
	figT.savefig('PLOT/Profiles-T.pdf')
	figU.savefig('PLOT/Profiles-U.pdf')
	figZ.savefig('PLOT/Profiles-Z.pdf')
else :
	plt.show()