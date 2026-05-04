#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%=================================================================================
from IPython import get_ipython

ip = get_ipython()
if ip is not None:
    ip.run_line_magic("load_ext", "autoreload")
    ip.run_line_magic("autoreload", "2")
#%%=================================================================================
#                     Modules
#===================================================================================
import Utilities as util
(Sysa,NSysa,Arg)=util.Parseur(['Save','Visu'],0,'Arg : ')
(                             [ SAVE , VISU ])=Arg

from numpy import *
import sys
import time
import Fluent as fl

t0=time.time()
(plt,mtp)=util.Plot0()
#%%=================================================================================
#                     Parameters
#===================================================================================
cmesh=1e3

from ParamsJaravel import *
dirr='/mnt/d/Python/Sandia/REF-Jaravel/'
Lines=['l0','l1','l2','l3','l75','l15','l30']
# dird='/mnt/scratch/ZEUS/FLUENT/Sandia-Jaravel/RUN-D100-01-EDM/DUMP/DATA/'
# dird='/mnt/scratch/ZEUS/FLUENT/Sandia-Jaravel/RUN-D100-02-Laera/DUMP/DATA/'
dird='/mnt/scratch/ZEUS/FLUENT/Sandia-Jaravel/RUN-D100-02-Laera/DUMP-07-EDC-NOx-FD39-O/DATA/'
# dird='/mnt/scratch/ZEUS/FLUENT/Sandia-Jaravel/RUN-D100-01-EDM/DUMP-03-EDC-BFER-NOx-Turb-Phi1-Tp/DATA/'

dirp=dird+'PLOT/'

D=D0
#%%=================================================================================
#                     reading
#===================================================================================
def PlotIm(ax,pic,Ext) :
	im=mtp.image.imread(pic)
	(Ni,Nj,Nk)=im.shape
	ax.imshow( im,extent=Ext )
	ax.set_aspect( (Ext[1]/Ext[3])*(Ni/Nj) )
#===================================================================================
util.MKDIR(dirp)
#===================================================================================
if VISU :
	util.Section( 'Visualisation : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	Vp=[ Lc+n*D0 for n in [1,2,3,7.5,15,30] ]

	fl.Visu(dird+'Data-all.dat','velocity-magnitude','Velocity [m/s]'        ,[0,0.3],[],arange(0,100,10)    ,cmesh,[],(25,5),'cividis',dirp+'Visu-Velocity.png'   ,[])
	fl.Visu(dird+'Data-all.dat','temperature'       ,'Temperature [K]'       ,[0,0.3],[],arange(250,2500,250),cmesh,[],(25,5),'inferno',dirp+'Visu-Temperature.png',['LINES',Vp])
	fl.Visu(dird+'Data-all.dat','mixC'              ,'Mixture fraction [-]'  ,[0,0.3],[],arange(0,1.1,0.1)   ,cmesh,[],(25,5),'viridis',dirp+'Visu-Mix.png'        ,['MIXC',[fl.Mol_m,BC_m,BC_o]])
	fl.Visu(dird+'Data-all.dat','co'                ,'Carbon monoxide [x100]',[0,0.3],[],arange(0,11,1)      ,cmesh,[],(25,5),'viridis',dirp+'Visu-CO.png'         ,['CO',1e2])
	fl.Visu(dird+'Data-all.dat','no'                ,r'NO [$x10^6$]'         ,[0,0.3],[],arange(0,100,10)    ,cmesh,[],(25,5),'viridis',dirp+'Visu-NO.png'         ,['NO',1e6])

	# sys.exit('=> End of visualisation')

#%%=================================================================================

figA,axA=plt.subplots(ncols=2,        figsize=(20,10)) #; figA.suptitle('Axial profiles'     ,fontsize=30)
figB,axB=plt.subplots(ncols=2,        figsize=(20,10)) #; figA.suptitle('Axial profiles'     ,fontsize=30)
figT,axT=plt.subplots(ncols=2,nrows=3,figsize=(20,10)) ; figT.suptitle('Radial temperature profiles'     ,fontsize=30)
figU,axU=plt.subplots(ncols=2,nrows=3,figsize=(20,10)) ; figU.suptitle('Radial velocity profiles'        ,fontsize=30)
figZ,axZ=plt.subplots(ncols=2,nrows=3,figsize=(20,10)) ; figZ.suptitle('Radial mixture fraction profiles',fontsize=30)
figC,axC=plt.subplots(ncols=2,nrows=3,figsize=(20,10)) ; figC.suptitle('Radial CO profiles'              ,fontsize=30)
figN,axN=plt.subplots(ncols=2,nrows=3,figsize=(20,10)) ; figN.suptitle('Radial NO profiles'              ,fontsize=30)
axA[0].set_title('Mean mixture fraction',fontsize=30)
axA[1].set_title('Mean temperature'     ,fontsize=30)
axB[0].set_title('Carbon monoxide',fontsize=30)
axB[1].set_title('Nitroux oxide'  ,fontsize=30)
Col=len(Lines)*['r']

#==========> Read Simulation
(M)=fl.ReadSurfD(dird+'Data-l0.dat')
([Vx,Vy])=fl.SpaceD(M)
#==========> Axial Profiles A
ym=38
Ym=[1.1,2e3]
Yt=[[0,0.5,1],arange(0,2100,500)]
for n,v in enumerate(['Z','T']) :
	#=====> Read Jaravel
	PlotIm( axA[n],dirr+'0D-{}.png'.format(v),[0,ym,0,Ym[n]] )
	axA[n].set_yticks(Yt[n])
	axA[n].set_xlim((0,ym))
#===> Mixture Fraction
Yc_f=fl.Yc( BC_m , fl.Mol_m )
Yc_p=fl.Yc( BC_p , fl.Mol_m )
Yc_o=fl.Yc( BC_o , fl.Mol_m )
Yc  =fl.Yc( {'CH4':M['ch4'],'CO2':M['co2'],'CO':M['co']} , fl.Mol_m )
Za=(Yc-Yc_o)/(Yc_f-Yc_o)
axA[0].plot( Vx/D,Za,'r' )
axA[0].set_ylabel('<Z> [-]',fontsize=30)
#===> Temperature
axA[1].plot( Vx/D,M['temperature'],'r' )
axA[1].set_ylabel('<T> [K]',fontsize=30)
#==========> Axial Profiles B
y0=-0.7
Ym=[6,7]
Yt=[[0,2,4,6],[0,2,4,6]]
for n,v in enumerate(['CO','NO']) :
	#=====> Read Jaravel
	PlotIm( axB[n],dirr+'0D-{}.png'.format(v),[0,ym,y0,Ym[n]] )
	axB[n].set_yticks(Yt[n])
	axB[n].set_xlim(( 0,ym   ))
	axB[n].set_ylim((y0,Ym[n]))
#===> CO
axB[0].plot( Vx/D,M['co']*1e2,'r' )
axB[0].set_ylabel(r'<$Y_{CO}$> [$x100$]',fontsize=30)
#===> NO
axB[1].plot( Vx/D,M['mf-pollut-pollutant-0']*1e6,'r' )
axB[1].set_ylabel(r'<$Y_{NO}$> [$x10^6$]',fontsize=30)

#%%=================================================================================

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
	if   i==0 : rtxt=2 ; yt=2100 ; yn=40           ; rm=2.5
	elif i==1 : rtxt=4 ; yt=2200 ; yn=60+20*(j==0) ; rm=5
	if j==1 : 
		axT[j,i].set_ylabel(r'<T> [K]'       ,fontsize=30)
		axU[j,i].set_ylabel(r'<U> [m/s]'     ,fontsize=30)
		axZ[j,i].set_ylabel(r'<Z> [-]'       ,fontsize=30)
		axC[j,i].set_ylabel(r'<$Y_{CO}$> [%]',fontsize=30)
		axN[j,i].set_ylabel(r'<$Y_{NO}$> [-]',fontsize=30)
	PlotIm(axT[j,i],dirr+l[1:]+'D-T.png' ,[0,rm,0,yt])
	PlotIm(axU[j,i],dirr+l[1:]+'D-U.png' ,[0,rm,0,yu])
	PlotIm(axZ[j,i],dirr+l[1:]+'D-Z.png' ,[0,rm,0,yz])
	PlotIm(axC[j,i],dirr+l[1:]+'D-CO.png',[0,rm,0,yc])
	PlotIm(axN[j,i],dirr+l[1:]+'D-NO.png',[0,rm,0,yn])
	axT[j,i].set_yticks(arange(0,2100,400))
	#=====> Read simulation
	(T,M)=fl.ReadSurf(dird+'Data-{}.dat'.format(l))
	([Ix,Iy],[Vx,Vy])=fl.Space(T,M)
	fl.PlotVar(axT[j,i],Col[n],'temperature'       ,Vy/D,M,T,[[rtxt,1800],txt+'D'],[0,rm,0,2250],[])
	fl.PlotVar(axU[j,i],Col[n],'velocity-magnitude',Vy/D,M,T,[[rtxt,50  ],txt+'D'],[0,rm,0,62  ],[])
	fl.PlotVar(axZ[j,i],Col[n],'mix'               ,Vy/D,M,T,[[rtxt,0.75],txt+'D'],[0,rm,0,yz  ],[BC_m,BC_o])
	fl.PlotVar(axC[j,i],Col[n],'co'                ,Vy/D,M,T,[[rtxt,3   ],txt+'D'],[0,rm,0,yc  ],[])
	fl.PlotVar(axN[j,i],Col[n],'no'                ,Vy/D,M,T,[[rtxt,0   ],     ''],[0,rm,0,yn  ],[])
axU[0,0].plot(2*[0.5      ],[0,60],':k')
axU[0,0].plot(2*[0.5*D1/D0],[0,60],':k')

if SAVE :
	figA.tight_layout() ; figA.savefig(dirp+'Profiles-Axis-A.pdf')
	figB.tight_layout() ; figB.savefig(dirp+'Profiles-Axis-B.pdf')
	figT.tight_layout() ; figT.savefig(dirp+'Profiles-T.pdf')
	figU.tight_layout() ; figU.savefig(dirp+'Profiles-U.pdf')
	figZ.tight_layout() ; figZ.savefig(dirp+'Profiles-Z.pdf')
	figC.tight_layout() ; figC.savefig(dirp+'Profiles-CO.pdf')
	figN.tight_layout() ; figN.savefig(dirp+'Profiles-NO.pdf')
	print('=> Figure in :'+dirp)
# else :
# 	plt.show()

# %%