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

# dird='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-00-Small/DUMP-04-UCSD-EDC/DATA/'
# dird='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-01-Big/DUMP-03-UCSD-EDC/DATA/'
# dird='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-01-Big/DUMP-04-UCSD-EDC-DF4/DATA/'
# dird='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-01-Big/DUMP-05-UCSD-EDC-Tw800/DATA/'
dird='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-01-Big/DUMP-06-UCSD-EDC-Tw500/DATA/'
dirp=dird+'PLOT/'

D=D0
Ld=180
#===================================================================================
#                     reading
#===================================================================================
util.MKDIR(dirp)
#===================================================================================
if VISU :
	util.Section( 'Visualisation : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	Vp=[ Lc+n*D0*Ld for n in Pos ]
	Tiso=list(array([1600,1800])-0)

	# fv=fl.Visu(dird+'Data-all.dat','velocity-magnitude','Velocity [m/s]'      ,[0.01,0.7],[0,0.05],arange(0,350,50)    ,[],(25,5),'cividis',dirp+'Visu-Velocity.png'   ,[])
	fl.Visu(dird+'Data-all.dat','temperature'       ,'Temperature [K]'     ,[0.01,0.7],[0,0.13],arange(250,2500,250),[],(25,5),'inferno',dirp+'Visu-Temperature.png',['LINES',Vp,'ISO',Tiso])
	# fh=fl.Visu(dird+'Data-all.dat','mixH'              ,'Mixture fraction [-]',[0.01,0.7],[0,0.15],arange(0,1.1,0.1)   ,[],(25,5),'viridis',dirp+'Visu-Mix.png'        ,['MIXH',[fl.Mol_m,BC_f,BC_o]])
	# ff=fl.Visu(dird+'Data-all.dat','h2'                ,'Y H2 [-]'            ,[0.01,0.7],[0,0.15],arange(0,1.1,0.1)   ,[],(25,5),'viridis',dirp+'Visu-H2.png'         ,[])
	# fp=fl.Visu(dird+'Data-all.dat','h2o'               ,'Y H2O [-]'           ,[0.01,0.7],[0,0.15],arange(0,1.1,0.1)   ,[],(25,5),'viridis',dirp+'Visu-H2O.png'        ,[])

	sys.exit('=> End of visualisation')

#===================================================================================

figZ,axZ=plt.subplots(figsize=(13,10)) ; figZ.suptitle('Mixture fraction [-]',fontsize=30)
figT,axT=plt.subplots(figsize=(13,10)) ; figT.suptitle('Temperature [K]'     ,fontsize=30)
# figN,axN=plt.subplots(figsize=(13,10)) ; figN.suptitle('X_NO [-]'            ,fontsize=30)

#=====> Parameters
rtxt=14
rs=16
tm=2500
zm=0.12
Col=['c','b','y','gold','orange','m','r']

#=====> Reference
util.PlotIm(axT,dirr+'T.png'  ,[0,rs,0,tm])
util.PlotIm(axZ,dirr+'Mix.png',[0,rs,0,zm])
for n,l in enumerate(Lines[1:]) :
	#=====> Data
	(T,M)=fl.ReadSurf(dird+'Data-{}.dat'.format(l)) ; LM=len(M)
	if LM>0 :
		([Ix,Iy],[Vx,Vy])=fl.Space(T,M)
		fl.PlotVar(axT,Col[n],'temperature',Vy/D,M,T,[[rtxt,2000],''],[0,rs,0,tm],[])
		fl.PlotVar(axZ,Col[n],'mix'        ,Vy/D,M,T,[[rtxt,0.1 ],''],[0,rs,0,zm],[BC_f,BC_o])
	else : print('=> No data for line {}'.format(l))

if SAVE :
	figT.savefig(dirp+'Profiles-T.pdf')
	figZ.savefig(dirp+'Profiles-Z.pdf')
else :
	plt.show()