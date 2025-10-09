#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Utilities as util
(Sysa,NSysa,Arg)=util.Parseur(['Save','Visu'],0,'Arg : ')
(                             [ SAVE , VISU ])=Arg

from numpy import *
import sys
import time
import Fluent as fl

t0=time.time()
(plt,mtp)=util.Plot0()
#===================================================================================
#                     Parameters
#===================================================================================
Fuel='H2'
# Fuel='CH4'
if   Fuel=='H2' :
	from ParamsGarnier import *
	dirr='/mnt/d/Python/SandiaJaravel/REF-Garnier/'
	Lines=['l00','l18','l14','l38','l12','l58','l34','l11']
elif Fuel=='CH4' :
	from ParamsJaravel import *
	dirr='/mnt/d/Python/SandiaJaravel/REF-Jaravel/'
	Lines=['l0','l1','l2','l3','l75','l15','l30']
else :
	sys.exit('=> Error : Fuel not recognized')

dird='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-00-Small/DUMP-04-UCSD-EDC/DATA/'
dirp=dird+'PLOT/'

D=D0
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
		Yc=fl.Yc( {'CH4':M[:,Ic1],'CO2':M[:,Ic2],'CO':M[:,Ic3]} , fl.Mol_m )
		Var=(Yc-Yc_o)/(Yc_f-Yc_o)
	elif var=='co' :
		I=T.index('co')
		Var=M[:,I]*1e2
	else :
		I=T.index(var)
		Var=M[:,I]
	ax.plot( X,Var,'r' )
	ax.text( Xtxt[0],Xtxt[1],txt )
	ax.set_xlim((BD[0],BD[1]))
	ax.set_ylim((BD[2],BD[3]))
#===================================================================================
util.MKDIR(dirp)
#===================================================================================
if VISU :
	util.Section( 'Visualisation : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	Vp=[ Lc+n*D0 for n in [1,2,3,7.5,15,30] ]

	fl.Visu(dird+'Data-all.dat','velocity-magnitude','Velocity [m/s]'        ,[0,0.3],[],arange(0,100,10)    ,[],(25,5),'cividis',dirp+'Visu-Velocity.png'   ,[])
	fl.Visu(dird+'Data-all.dat','temperature'       ,'Temperature [K]'       ,[0,0.3],[],arange(250,2500,250),[],(25,5),'inferno',dirp+'Visu-Temperature.png',['LINES',Vp])
	fl.Visu(dird+'Data-all.dat','mixC'              ,'Mixture fraction [-]'  ,[0,0.3],[],arange(0,1.1,0.1)   ,[],(25,5),'viridis',dirp+'Visu-Mix.png'        ,['MIXC',[fl.Mol_m,Y_m,Y_p,Y_o]])
	fl.Visu(dird+'Data-all.dat','co'                ,'Carbon monoxide [x100]',[0,0.3],[],arange(0,11,1)      ,[],(25,5),'viridis',dirp+'Visu-CO.png'         ,['CO',1e2])

	# sys.exit('=> End of visualisation')

#===================================================================================

figA,axA=plt.subplots(ncols=2,        figsize=(20,10)) #; figA.suptitle('Axial profiles'     ,fontsize=30)
figT,axT=plt.subplots(ncols=2,nrows=3,figsize=(20,10)) ; figT.suptitle('Radial temperature profiles'     ,fontsize=30)
figU,axU=plt.subplots(ncols=2,nrows=3,figsize=(20,10)) ; figU.suptitle('Radial velocity profiles'        ,fontsize=30)
figZ,axZ=plt.subplots(ncols=2,nrows=3,figsize=(20,10)) ; figZ.suptitle('Radial mixture fraction profiles',fontsize=30)
figC,axC=plt.subplots(ncols=2,nrows=3,figsize=(20,10)) ; figC.suptitle('Radial CO profiles'              ,fontsize=30)
axA[0].set_title('Mean mixture fraction',fontsize=30)
axA[1].set_title('Mean temperature'     ,fontsize=30)

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
Yc_f=fl.Yc( Y_m , fl.Mol_m )
Yc_p=fl.Yc( Y_p , fl.Mol_m )
Yc_o=fl.Yc( Y_o , fl.Mol_m )
Yc  =fl.Yc( {'CH4':M[:,Ic1],'CO2':M[:,Ic2],'CO':M[:,Ic3]} , fl.Mol_m )
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
	yc=4
	if   i==0 : rtxt=2 ; yt=2100 ; rm=2.5
	elif i==1 : rtxt=4 ; yt=2200 ; rm=5
	if j==1 : 
		axT[j,i].set_ylabel('<T> [K]'  ,fontsize=30)
		axU[j,i].set_ylabel('<U> [m/s]',fontsize=30)
		axZ[j,i].set_ylabel('<Z> [-]'  ,fontsize=30)
	PlotIm(axT[j,i],dirr+l[1:]+'D-T.png' ,[0,rm,0,yt])
	PlotIm(axU[j,i],dirr+l[1:]+'D-U.png' ,[0,rm,0,yu])
	PlotIm(axZ[j,i],dirr+l[1:]+'D-Z.png' ,[0,rm,0,yz])
	PlotIm(axC[j,i],dirr+l[1:]+'D-CO.png',[0,rm,0,yc])
	axT[j,i].set_yticks(arange(0,2100,400))
	#=====> Read simulation
	(T,M)=fl.ReadSurf(dird+'Data-{}.dat'.format(l))
	([Ix,Iy],[Vx,Vy])=fl.Space(T,M)
	PlotVar(axT[j,i],'temperature'       ,Vy/D,M,T,[[rtxt,1800],txt+'D'],[0,rm,0,2250])
	PlotVar(axU[j,i],'velocity-magnitude',Vy/D,M,T,[[rtxt,50  ],txt+'D'],[0,rm,0,62  ])
	PlotVar(axZ[j,i],'mix'               ,Vy/D,M,T,[[rtxt,0.75],txt+'D'],[0,rm,0,yz  ])
	PlotVar(axC[j,i],'co'                ,Vy/D,M,T,[[rtxt,3   ],txt+'D'],[0,rm,0,yc  ])
axU[0,0].plot(2*[0.5      ],[0,60],':k')
axU[0,0].plot(2*[0.5*D1/D0],[0,60],':k')

if SAVE :
	figA.savefig(dirp+'Profiles-Axis.pdf')
	figT.savefig(dirp+'Profiles-T.pdf')
	figU.savefig(dirp+'Profiles-U.pdf')
	figZ.savefig(dirp+'Profiles-Z.pdf')
	figC.savefig(dirp+'Profiles-CO.pdf')
else :
	plt.show()