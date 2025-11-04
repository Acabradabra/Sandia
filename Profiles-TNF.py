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
(Sysa,NSysa,Arg)=util.Parseur(['Save','Visu','Temp','Data','Over','Vel','Only','All'],0,'Arg : ')
(                             [ SAVE , VISU , TEMP , DATA , OVER , VEL , ONLY , ALL ])=Arg
if ALL : SAVE,VISU,TEMP,DATA,OVER=True,True,True,True,True

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
TICKS=True
# TICKS=False

#====================> Burner
# Burner='Pilot'
Burner='Simple'
if   Burner=='Simple' :	from ParamsGarnier import *
elif Burner=='Pilot'  : from ParamsJaravel import *
else : sys.exit('=> Error : Burner not recognized')

#====================> Helium content 
# he='he4'
he=0.4
if   he==0   : Ld=180 ; Umoy=296
elif he==0.2 : Ld=150 ; Umoy=294
elif he==0.4 : Ld=100 ; Umoy=256

#====================> Directories
dir0='/mnt/d/Python/Sandia/'
dirv=dir0+'DATA-H2/ETHZ_H2/%.0fhe/'%(he*100)
dirs=dir0+'DATA-H2/SANDH2_A/he%.0f/he%.0fstatY/'%(he*10,he*10)
# dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-02-Big-40p/DUMP/'
# dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-02-Big-40p/DUMP-00-He-FD3/'
dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-02-Big-40p/DUMP-01-He-Uprof2/'
# dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-02-Big-40p/DUMP-04-He-FD39/'
dird=dirc+'DATA/'
dirp=dirc+'PLOT/'

#====================> Positions Scattering
Files_s=os.popen('ls %s/he%.0f*.fav'%(dirs,he*10)).read().splitlines()
Pos_s=[]
for f in Files_s :
	pos0=f.split('/')[-1].split('.')[0]
	if 'x' in pos0 :
		pos=pos0.split('x')[1]
		xl=float(pos[0])/float(pos[1]) ; Pos_s.append(xl)
		print('=> Position  : x/L_vis={}'.format(xl))
Npos_s=len(Pos_s)
IPos_s=argsort(Pos_s)

#====================> Positions Velocimetry
Files_v=os.popen('ls %s/s*%.0f.dat'%(dirv,he*100)).read().splitlines()
# Pos_v=[0,1/16,1/8,1/4,3/8,1/2,5/8,3/4,1]
Pos_v=[0,1,1/2,1/4,1/16,1/8,3/4,3/8,5/8]
Npos_v=len(Pos_v)
IPos_v=argsort(Pos_v)

#====================> Fields
# Vars=['Vel','k','T','mixH','o2','h2','n2','h2o']
Vars=['Vel','k']
# Vars=['Vel']
# Vars=['T']

#====================> Velocity profile
n=15
Umax=Umoy*(n+2)/n ; print('=> n = %.0f  ,  Umoy = %.3f [m/s]  ,  Umax = %.3f [m/s]'%(n,Umoy,Umax))
Np=int(1e4)
VyV=linspace(0,0.5*D0,Np)
Uth=Umax*(1-(2*VyV/D0)**n)

#%%=================================================================================
#                     reading
#===================================================================================
util.Entete1(104,[dirc,'he {:.0f} %'.format(he*100)],'Profiles TNF')
#===================================================================================
util.MKDIR(dirp)
#===================================================================================
if TEMP :
	util.Section( 'Temporals : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	fl.Probe_plot(dirc+'probe-F.out',dirp+'Probe-F.pdf')
	fl.Probe_plot(dirc+'probe-T.out',dirp+'Probe-T.pdf')
	fl.Probe_plot(dirc+'probe-V.out',dirp+'Probe-V.pdf')
	(D_in,D_ou)=fl.Report_read(dirc+'InletMassFlowRate.out',dirc+'outletMassFlowRate.out')
	It =D_in[:,0]
	Min=D_in[:,1]+D_in[:,2]
	Mou=D_ou[:,1]+D_ou[:,2]
	fig_b,ax_b=fl.Report(It,Min,Mou,0,0,1e9,int(1e2),1e3)
	fig_b.tight_layout()
	fig_b.savefig(dirp+'Balance.pdf')
	if ONLY : sys.exit('=> Stop after temporals')
#%%=================================================================================
if VISU :
	util.Section( 'Visualisation : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	Vp_v=[ Lc+n*D0*Ld for n in Pos_v ]
	Vp_s=[ Lc+n*D0*Ld for n in Pos_s ]
	Tiso=list(array([1600,1800])-0)

	F_int={}
	if 'Vel'  in Vars : F_int['Vel' ]=fl.Visu(dird+'Data-all.dat','velocity-magnitude' ,'Velocity [m/s]'      ,[ 0,0.2],[0,0.04],arange(0,350,50)    ,1e3,[],(25,5),'cividis',dirp+'Visu-Velocity.png'   ,['INTERP','LINES',Vp_v])
	if 'k'    in Vars : F_int['k'   ]=fl.Visu(dird+'Data-all.dat','turb-kinetic-energy','k [$m^2/s^2$]'       ,[ 0,0.2],[0,0.04],arange(250,2500,250),1e3,[],(25,5),'cividis',dirp+'Visu-TKE.png'        ,['INTERP'])
	if 'mixH' in Vars : F_int['mixH']=fl.Visu(dird+'Data-all.dat','mixH'               ,'Mixture fraction [-]',[Lc,0.2],[0,0.04],arange(0,1.1,0.1)   ,1e3,[],(25,5),'viridis',dirp+'Visu-Mix.png'        ,['INTERP','MIXH',[fl.Mol_m,BC_f,BC_o]])
	if 'h2'   in Vars : F_int['h2'  ]=fl.Visu(dird+'Data-all.dat','h2'                 ,'Y H2 [-]'            ,[Lc,0.2],[0,0.04],arange(0,1.1,0.1)   ,1e3,[],(25,5),'viridis',dirp+'Visu-H2.png'         ,['INTERP'])
	if 'n2'   in Vars : F_int['n2'  ]=fl.Visu(dird+'Data-all.dat','n2'                 ,'Y N2 [-]'            ,[Lc,0.2],[0,0.04],arange(0,1.1,0.1)   ,1e3,[],(25,5),'viridis',dirp+'Visu-N2.png'         ,['INTERP'])
	if 'o2'   in Vars : F_int['o2'  ]=fl.Visu(dird+'Data-all.dat','o2'                 ,'Y O2 [-]'            ,[Lc,0.4],[0,0.08],arange(0,1.1,0.1)   ,1e3,[],(25,5),'viridis',dirp+'Visu-O2.png'         ,['INTERP'])
	if 'h2o'  in Vars : F_int['h2o' ]=fl.Visu(dird+'Data-all.dat','h2o'                ,'Y H2O [-]'           ,[Lc,0.4],[0,0.08],arange(0,1.1,0.1)   ,1e3,[],(25,5),'viridis',dirp+'Visu-H2O.png'        ,['INTERP'])
	if 'T'    in Vars : F_int['T'   ]=fl.Visu(dird+'Data-all.dat','temperature'        ,'Temperature [K]'     ,[Lc,0.4],[0,0.08],arange(250,2500,250),1e3,[],(25,5),'inferno',dirp+'Visu-Temperature.png',['INTERP','LINES',Vp_s,'ISO',Tiso])

	D_int={}
	Vy=linspace(0,0.5*D2,Np)
	for v in Vars :
		if v in ['Vel','k'] : D_int[v]=[ F_int[v](Np*[p],Vy) for p in Vp_v ]
		else                : D_int[v]=[ F_int[v](Np*[p],Vy) for p in Vp_s ]

	if DATA :
		util.Section( 'Data extraction : {:.3f} s'.format(time.time()-t0),1,5,'r' )
		for v in Vars :
			St_data=dird+'Profiles-%s.csv'%(v)
			if os.path.exists(St_data) and not OVER : util.Section('Warning : Data file %s already exists'%(St_data),0,1,'y')
			else :
				if os.path.exists(St_data) : util.Section('Overwriting : %s '%(St_data),0,1,'y')
				with open(St_data,'w') as file :
					writer=csv.writer(file)
					writer.writerow( ['r(mm)'] + [ 'x/L_vis=%.2f'%(p) for p in Pos_s ] )
					writer.writerows( [ [ Vy[i]*1e3 ] + [ D_int[v][j][i] for j in IPos_s ] for i in range(Np) ] )

	if VEL and 'Vel' in Vars :
		util.Section( 'Velocity profiles : {:.3f} s'.format(time.time()-t0),1,5,'r' )
		PosV=array([0,0.1,0.5,0.9,1])*Lc
		figVp,axVp=plt.subplots(figsize=(10,7)) ; figVp.suptitle('Velocity profiles'       ,fontsize=30)
		figKp,axKp=plt.subplots(figsize=(10,7)) ; figKp.suptitle('Turbulent kinetic energy',fontsize=30)
		axVp.set_xlabel('r [mm]'      ,fontsize=30) ; axKp.set_xlabel('r [mm]'       ,fontsize=30)
		axVp.set_ylabel('V [m/s]'     ,fontsize=30) ; axKp.set_ylabel('k [$m^2/s^2$]',fontsize=30)
		for p in PosV :
			Vprof=F_int['Vel'](Np*[p],VyV)
			Kprof=F_int['k'  ](Np*[p],VyV)
			print('=> x/Lx=%.2f  ,  Max V=%.3f [m/s]  ,  Mean V=%.3f [m/s]'%(p/Lc,max(Vprof),fl.Umoy(VyV,Vprof)) ) 
			axVp.plot( VyV*1e3,Vprof , label='x/Lc=%.2f'%(p/Lc) )
			axKp.plot( VyV*1e3,Kprof , label='x/Lc=%.2f'%(p/Lc) )
		print( '=> Poiseuil  ,  Max V=%.3f [m/s]  ,  Mean V=%.3f [m/s]'%(max(Uth),fl.Umoy(VyV,Uth)) ) 
		axVp.plot( VyV*1e3,Uth , 'k--' , label='Theoretical' )
		axVp.legend(fontsize=20) ; util.SaveFig(figVp,dirp+'Inlet-Velocity.pdf')
		axKp.legend(fontsize=20) ; util.SaveFig(figKp,dirp+'Inlet-TKE.pdf'     )
#%%=================================================================================
#                     reading
#===================================================================================
def ReadTNF(name) :
	sep=' ' ; skip=3
	with open(name) as file :
		Lines=file.readlines()
		T=[ s.strip() for s in Lines[skip][:-1].split(sep) if s ]
		M=array([ [float(v) for v in L.split(sep) if v ] for L in Lines[skip+1:] ])
	file.closed
	return({ s:M[:,i] for i,s in enumerate(T) })
#===================================================================================
def ReadLDV(name) :
	sep='\t' ; skip=11
	with open(name) as file :
		Lines=file.readlines()
		T=[ s.strip() for s in Lines[skip][2:-1].split(sep) if s ]
		M=array([ [float(v) for v in L.split(sep) if v ] for L in Lines[skip+1:] ])
	file.closed
	return({ s:M[:,i] for i,s in enumerate(T) })
#%%===================================================================================
nc=2 ; nr=Npos_s//nc+int(Npos_s%nc>0)

if 'Vel'  in Vars : figV,axV=plt.subplots(ncols=3 ,nrows=3 ,figsize=(13,10)) ; figV.suptitle('Velocity [m/s]'      ,fontsize=30)
if 'k'    in Vars : figK,axK=plt.subplots(ncols=3 ,nrows=3 ,figsize=(13,10)) ; figK.suptitle('TKE [$m^2/s^2$]'     ,fontsize=30)
if 'T'    in Vars : figT,axT=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figT.suptitle('Temperature [K]'     ,fontsize=30)
if 'mixH' in Vars : figZ,axZ=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figZ.suptitle('Mixture fraction [-]',fontsize=30)
if 'o2'   in Vars : figO,axO=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figO.suptitle('Y O2 [-]'            ,fontsize=30)
if 'h2'   in Vars : figH,axH=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figH.suptitle('Y H2 [-]'            ,fontsize=30)
if 'n2'   in Vars : figN,axN=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figN.suptitle('Y N2 [-]'            ,fontsize=30)
if 'h2o'  in Vars : figP,axP=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figP.suptitle('Y H2O [-]'           ,fontsize=30)
if Npos_s%nc>0 :
	# if 'Vel'  in Vars : axV[-1,-1].axis('off')
	# if 'k'    in Vars : axK[-1,-1].axis('off')
	if 'T'    in Vars : axT[-1,-1].axis('off')
	if 'mixH' in Vars : axZ[-1,-1].axis('off')
	if 'o2'   in Vars : axO[-1,-1].axis('off')
	if 'h2'   in Vars : axH[-1,-1].axis('off')
	if 'n2'   in Vars : axN[-1,-1].axis('off')
	if 'h2o'  in Vars : axP[-1,-1].axis('off')

#=====> Velocity profiles
Rlims=[3,8,10,20,20,20,40,40,40]
Pv=[ '0','1/16','1/8','1/4','3/8','1/2','5/8','3/4','1' ]
for n,p in enumerate(IPos_v) :
	i= n//3
	j= n-3*i
	Data=ReadLDV(Files_v[p])
	Vel=hypot( Data['u'] , Data['v'] )
	Tke=hypot( Data['varu'] , Data['varv'] )
	if 'Vel' in Vars : axV[i,j].plot( Vy*1e3,D_int['Vel' ][p],'r' ) ; axV[i,j].plot( Data['y'],Vel,'ok' ) ; axV[i,j].set_title('x/L_vis='+Pv[n],fontsize=20) ; axV[i,j].ticklabel_format(style='sci',axis='y',scilimits=(-2,4)) ; axV[i,j].set_xlim((0,Rlims[n]))
	if 'k'   in Vars : axK[i,j].plot( Vy*1e3,D_int['k'   ][p],'r' ) ; axK[i,j].plot( Data['y'],Tke,'ok' ) ; axK[i,j].set_title('x/L_vis='+Pv[n],fontsize=20) ; axK[i,j].ticklabel_format(style='sci',axis='y',scilimits=(-2,4)) ; axK[i,j].set_xlim((0,Rlims[n]))
	if i==2 : 
		if 'Vel' in Vars : axV[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'k'   in Vars : axK[i,j].set_xlabel('r [mm]',fontsize=20)
if 'Vel' in Vars : axV[0,0].plot( [0,Rlims[0]],2*[Umoy],':k' )

#=====> Scattering
for n,p in enumerate(IPos_s) :
	i= n//nc
	j= n-nc*i
	Data=ReadTNF(Files_s[p])
	if 'T'    in Vars : axT[i,j].plot(Vy*1e3,D_int['T'   ][p],'r') ; axT[i,j].errorbar(Data['r(mm)'],Data['T(K)' ],yerr=0.03*Data['T(K)' ],ecolor='k',color='k') ; axT[i,j].set_title('x/L_vis=%.2f'%(Pos_s[p]),fontsize=20) ; axT[i,j].ticklabel_format(style='sci',axis='y',scilimits=(-2,4))
	if 'mixH' in Vars : axZ[i,j].plot(Vy*1e3,D_int['mixH'][p],'r') ; axZ[i,j].errorbar(Data['r(mm)'],Data['Fblgr'],yerr=0.04*Data['Fblgr'],ecolor='k',color='k') ; axZ[i,j].set_title('x/L_vis=%.2f'%(Pos_s[p]),fontsize=20) ; axZ[i,j].ticklabel_format(style='sci',axis='y',scilimits=(-2,4))
	if 'o2'   in Vars : axO[i,j].plot(Vy*1e3,D_int['o2'  ][p],'r') ; axO[i,j].errorbar(Data['r(mm)'],Data['YO2'  ],yerr=0.04*Data['YO2'  ],ecolor='k',color='k') ; axO[i,j].set_title('x/L_vis=%.2f'%(Pos_s[p]),fontsize=20) ; axO[i,j].ticklabel_format(style='sci',axis='y',scilimits=(-2,4))
	if 'h2'   in Vars : axH[i,j].plot(Vy*1e3,D_int['h2'  ][p],'r') ; axH[i,j].errorbar(Data['r(mm)'],Data['YH2'  ],yerr=0.04*Data['YH2'  ],ecolor='k',color='k') ; axH[i,j].set_title('x/L_vis=%.2f'%(Pos_s[p]),fontsize=20) ; axH[i,j].ticklabel_format(style='sci',axis='y',scilimits=(-2,4))
	if 'n2'   in Vars : axN[i,j].plot(Vy*1e3,D_int['n2'  ][p],'r') ; axN[i,j].errorbar(Data['r(mm)'],Data['YN2'  ],yerr=0.04*Data['YN2'  ],ecolor='k',color='k') ; axN[i,j].set_title('x/L_vis=%.2f'%(Pos_s[p]),fontsize=20) ; axN[i,j].ticklabel_format(style='sci',axis='y',scilimits=(-2,4))
	if 'h2o'  in Vars : axP[i,j].plot(Vy*1e3,D_int['h2o' ][p],'r') ; axP[i,j].errorbar(Data['r(mm)'],Data['YH2O' ],yerr=0.04*Data['YH2O' ],ecolor='k',color='k') ; axP[i,j].set_title('x/L_vis=%.2f'%(Pos_s[p]),fontsize=20) ; axP[i,j].ticklabel_format(style='sci',axis='y',scilimits=(-2,4))
	if i==nr-1 : 
		if 'T'    in Vars : axT[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'mixH' in Vars : axZ[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'o2'   in Vars : axO[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'h2'   in Vars : axH[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'n2'   in Vars : axN[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'h2o'  in Vars : axP[i,j].set_xlabel('r [mm]',fontsize=20)

if SAVE :
	if 'Vel'  in Vars : util.SaveFig(figV,dirp+'Profiles-V.pdf')
	if 'k'    in Vars : util.SaveFig(figK,dirp+'Profiles-K.pdf')
	if 'T'    in Vars : util.SaveFig(figT,dirp+'Profiles-T.pdf')
	if 'mixH' in Vars : util.SaveFig(figZ,dirp+'Profiles-Z.pdf')
	if 'o2'   in Vars : util.SaveFig(figO,dirp+'Profiles-O.pdf')
	if 'h2'   in Vars : util.SaveFig(figH,dirp+'Profiles-H.pdf')
	if 'n2'   in Vars : util.SaveFig(figN,dirp+'Profiles-N.pdf')
	if 'h2o'  in Vars : util.SaveFig(figP,dirp+'Profiles-P.pdf')
else :
	plt.show()
# %%