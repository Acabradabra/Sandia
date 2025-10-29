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
(Sysa,NSysa,Arg)=util.Parseur(['Save','Visu','Temp','Data','Over','Vel','All'],0,'Arg : ')
(                             [ SAVE , VISU , TEMP , DATA , OVER , VEL , ALL ])=Arg
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

#====================> Burner
# Burner='Pilot'
Burner='Simple'
if   Burner=='Simple' :	from ParamsGarnier import *
elif Burner=='Pilot'  : from ParamsJaravel import *
else : sys.exit('=> Error : Burner not recognized')

#====================> Helium content 
he='he4'
if   he=='he0' : Ld=180
elif he=='he2' : Ld=150
elif he=='he4' : Ld=100

#====================> Directories
dir0='/mnt/d/Python/Sandia/'
dirr=dir0+'DATA-H2/SANDH2_A/%s/%sstatY/'%(he,he)
# dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-02-Big-40p/DUMP-00-He-FD3/'
dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-02-Big-40p/DUMP-01-He-Uprof/'
# dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-02-Big-40p/DUMP-04-He-FD39/'
dird=dirc+'DATA/'
dirp=dirc+'PLOT/'

#====================> Positions
Files=os.popen('ls %s/%s*.fav'%(dirr,he)).read().splitlines()
Pos=[]
for f in Files :
	pos0=f.split('/')[-1].split('.')[0]
	if 'x' in pos0 :
		pos=pos0.split('x')[1]
		xl=float(pos[0])/float(pos[1]) ; Pos.append(xl)
		print('=> Position  : x/L_vis={}'.format(xl))
Npos=len(Pos)
IPos=argsort(Pos)

#====================> Fields
Vars=['Vel','T','mixH','o2','h2','n2','h2o']
# Vars=['Vel']
# Vars=['T']

#====================> Velocity profile
n=15
Umoy=256
Umax=Umoy*(n+2)/n ; print('=> n = %.0f  ,  Umoy = %.3f [m/s]  ,  Umax = %.3f [m/s]'%(n,Umoy,Umax))
Np=int(1e4)
VyV=linspace(0,0.5*D0,Np)
Uth=Umax*(1-(2*VyV/D0)**n)

#%%=================================================================================
#                     reading
#===================================================================================
util.Entete1(104,[dirc,he],'Profiles TNF')
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
#%%=================================================================================
if VISU :
	util.Section( 'Visualisation : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	Vp=[ Lc+n*D0*Ld for n in Pos ]
	Tiso=list(array([1600,1800])-0)

	F_int={}
	if 'Vel'  in Vars : F_int['Vel' ]=fl.Visu(dird+'Data-all.dat','velocity-magnitude','Velocity [m/s]'      ,[Lc,0.7],[0,0.15],arange(0,350,50)    ,[],(25,5),'cividis',dirp+'Visu-Velocity.png'   ,['INTERP'])
	if 'T'    in Vars : F_int['T'   ]=fl.Visu(dird+'Data-all.dat','temperature'       ,'Temperature [K]'     ,[Lc,0.7],[0,0.15],arange(250,2500,250),[],(25,5),'inferno',dirp+'Visu-Temperature.png',['INTERP','LINES',Vp,'ISO',Tiso])
	if 'mixH' in Vars : F_int['mixH']=fl.Visu(dird+'Data-all.dat','mixH'              ,'Mixture fraction [-]',[Lc,0.7],[0,0.15],arange(0,1.1,0.1)   ,[],(25,5),'viridis',dirp+'Visu-Mix.png'        ,['INTERP','MIXH',[fl.Mol_m,BC_f,BC_o]])
	if 'o2'   in Vars : F_int['o2'  ]=fl.Visu(dird+'Data-all.dat','o2'                ,'Y O2 [-]'            ,[Lc,0.7],[0,0.15],arange(0,1.1,0.1)   ,[],(25,5),'viridis',dirp+'Visu-O2.png'         ,['INTERP'])
	if 'h2'   in Vars : F_int['h2'  ]=fl.Visu(dird+'Data-all.dat','h2'                ,'Y H2 [-]'            ,[Lc,0.7],[0,0.15],arange(0,1.1,0.1)   ,[],(25,5),'viridis',dirp+'Visu-H2.png'         ,['INTERP'])
	if 'n2'   in Vars : F_int['n2'  ]=fl.Visu(dird+'Data-all.dat','n2'                ,'Y N2 [-]'            ,[Lc,0.7],[0,0.15],arange(0,1.1,0.1)   ,[],(25,5),'viridis',dirp+'Visu-N2.png'         ,['INTERP'])
	if 'h2o'  in Vars : F_int['h2o' ]=fl.Visu(dird+'Data-all.dat','h2o'               ,'Y H2O [-]'           ,[Lc,0.7],[0,0.15],arange(0,1.1,0.1)   ,[],(25,5),'viridis',dirp+'Visu-H2O.png'        ,['INTERP'])

	D_int={}
	Vy=linspace(0,0.5*D2,Np)
	for v in Vars : D_int[v]=[ F_int[v](Np*[p],Vy) for p in Vp ]

	if DATA :
		util.Section( 'Data extraction : {:.3f} s'.format(time.time()-t0),1,5,'r' )
		for v in Vars :
			St_data=dird+'Profiles-%s.csv'%(v)
			if os.path.exists(St_data) and not OVER : util.Section('Warning : Data file %s already exists'%(St_data),0,1,'y')
			else :
				if os.path.exists(St_data) : util.Section('Overwriting : %s '%(St_data),0,1,'y')
				with open(St_data,'w') as file :
					writer=csv.writer(file)
					writer.writerow( ['r(mm)'] + [ 'x/L_vis=%.2f'%(p) for p in Pos ] )
					writer.writerows( [ [ Vy[i]*1e3 ] + [ D_int[v][j][i] for j in IPos ] for i in range(Np) ] )

	if VEL and 'Vel' in Vars :
		util.Section( 'Velocity profiles : {:.3f} s'.format(time.time()-t0),1,5,'r' )
		PosV=array([0,0.1,0.5,0.9,1])*Lc
		figVp,axVp=plt.subplots(figsize=(10,7)) ; figVp.suptitle('Velocity profiles',fontsize=30)
		axVp.set_xlabel('r [mm]'      ,fontsize=30)
		axVp.set_ylabel('V [m/s]'     ,fontsize=30)
		for p in PosV :
			Vprof=F_int['Vel'](Np*[p],VyV)
			print('=> x/Lx=%.2f  ,  Max V=%.3f [m/s]  ,  Mean V=%.3f [m/s]'%(p/Lc,max(Vprof),fl.Umoy(VyV,Vprof)) ) 
			axVp.plot( VyV*1e3,Vprof , label='x/Lc=%.2f'%(p/Lc) )
		print( '=> Poiseuil  ,  Max V=%.3f [m/s]  ,  Mean V=%.3f [m/s]'%(max(Uth),fl.Umoy(VyV,Uth)) ) 
		axVp.plot( VyV*1e3,Uth , 'k--' , label='Theoretical' )
		axVp.legend(fontsize=20)
		figVp.tight_layout()
		figVp.savefig(dirp+'Profiles-Velocity.pdf')

#%%=================================================================================
#                     reading
#===================================================================================
def ReadTNF(name) :
	sep=' ' ; skip=3
	with open(name) as file :
		Lines=file.readlines()
		T=[ s.strip() for s in Lines[skip][:-1].split(' ') if s ]
		M=array([ [float(v) for v in L.split(sep) if v ] for L in Lines[skip+1:] ])
	file.closed
	return({ s:M[:,i] for i,s in enumerate(T) })
#===================================================================================
nc=2 ; nr=Npos//nc+int(Npos%nc>0)

if 'Vel'  in Vars : figV,axV=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figV.suptitle('Velocity [m/s]'      ,fontsize=30)
if 'T'    in Vars : figT,axT=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figT.suptitle('Temperature [K]'     ,fontsize=30)
if 'mixH' in Vars : figZ,axZ=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figZ.suptitle('Mixture fraction [-]',fontsize=30)
if 'o2'   in Vars : figO,axO=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figO.suptitle('Y O2 [-]'            ,fontsize=30)
if 'h2'   in Vars : figH,axH=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figH.suptitle('Y H2 [-]'            ,fontsize=30)
if 'n2'   in Vars : figN,axN=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figN.suptitle('Y N2 [-]'            ,fontsize=30)
if 'h2o'  in Vars : figP,axP=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figP.suptitle('Y H2O [-]'           ,fontsize=30)
if Npos%nc>0 :
	if 'Vel'  in Vars : axV[-1,-1].axis('off')
	if 'T'    in Vars : axT[-1,-1].axis('off')
	if 'mixH' in Vars : axZ[-1,-1].axis('off')
	if 'o2'   in Vars : axO[-1,-1].axis('off')
	if 'h2'   in Vars : axH[-1,-1].axis('off')
	if 'n2'   in Vars : axN[-1,-1].axis('off')
	if 'h2o'  in Vars : axP[-1,-1].axis('off')

for n,p in enumerate(IPos) :
	i= n//nc
	j= n-nc*i
	Data=ReadTNF(Files[p])
	# if 'Vel'  in Vars : axV[i,j].plot( Vy*1e3,D_int['Vel' ][n],'k' , Data['r(mm)',Data['']] )
	if 'T'    in Vars : axT[i,j].plot(Vy*1e3,D_int['T'   ][p],'k') ; axT[i,j].errorbar(Data['r(mm)'],Data['T(K)' ],yerr=0.03*Data['T(K)' ],ecolor='k',color='k') ; axT[i,j].text(0.5*max(Vy)*1e3,0.95*max(D_int['T'   ][p]),'x/L_vis=%.2f'%(Pos[p]),fontsize=20)
	if 'mixH' in Vars : axZ[i,j].plot(Vy*1e3,D_int['mixH'][p],'k') ; axZ[i,j].errorbar(Data['r(mm)'],Data['Fblgr'],yerr=0.04*Data['Fblgr'],ecolor='k',color='k') ; axZ[i,j].text(0.5*max(Vy)*1e3,0.95*max(D_int['mixH'][p]),'x/L_vis=%.2f'%(Pos[p]),fontsize=20)
	if 'o2'   in Vars : axO[i,j].plot(Vy*1e3,D_int['o2'  ][p],'k') ; axO[i,j].errorbar(Data['r(mm)'],Data['YO2'  ],yerr=0.04*Data['YO2'  ],ecolor='k',color='k') ; axO[i,j].text(0.5*max(Vy)*1e3,0.80*max(D_int['o2'  ][p]),'x/L_vis=%.2f'%(Pos[p]),fontsize=20)
	if 'h2'   in Vars : axH[i,j].plot(Vy*1e3,D_int['h2'  ][p],'k') ; axH[i,j].errorbar(Data['r(mm)'],Data['YH2'  ],yerr=0.04*Data['YH2'  ],ecolor='k',color='k') ; axH[i,j].text(0.5*max(Vy)*1e3,0.95*max(D_int['h2'  ][p]),'x/L_vis=%.2f'%(Pos[p]),fontsize=20)
	if 'n2'   in Vars : axN[i,j].plot(Vy*1e3,D_int['n2'  ][p],'k') ; axN[i,j].errorbar(Data['r(mm)'],Data['YN2'  ],yerr=0.04*Data['YN2'  ],ecolor='k',color='k') ; axN[i,j].text(0.5*max(Vy)*1e3,0.95*max(D_int['n2'  ][p]),'x/L_vis=%.2f'%(Pos[p]),fontsize=20)
	if 'h2o'  in Vars : axP[i,j].plot(Vy*1e3,D_int['h2o' ][p],'k') ; axP[i,j].errorbar(Data['r(mm)'],Data['YH2O' ],yerr=0.04*Data['YH2O' ],ecolor='k',color='k') ; axP[i,j].text(0.5*max(Vy)*1e3,0.95*max(D_int['h2o' ][p]),'x/L_vis=%.2f'%(Pos[p]),fontsize=20)
	if i==nr-1 : 
		if 'Vel'  in Vars : axV[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'T'    in Vars : axT[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'mixH' in Vars : axZ[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'o2'   in Vars : axO[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'h2'   in Vars : axH[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'n2'   in Vars : axN[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'h2o'  in Vars : axP[i,j].set_xlabel('r [mm]',fontsize=20)

if SAVE :
	if 'Vel'  in Vars : figV.savefig(dirp+'Profiles-V.pdf') ; print('=> Profiles-V.pdf saved')
	if 'T'    in Vars : figT.savefig(dirp+'Profiles-T.pdf') ; print('=> Profiles-T.pdf saved')
	if 'mixH' in Vars : figZ.savefig(dirp+'Profiles-Z.pdf') ; print('=> Profiles-Z.pdf saved')
	if 'o2'   in Vars : figO.savefig(dirp+'Profiles-O.pdf') ; print('=> Profiles-O.pdf saved')
	if 'h2'   in Vars : figH.savefig(dirp+'Profiles-H.pdf') ; print('=> Profiles-H.pdf saved')
	if 'n2'   in Vars : figN.savefig(dirp+'Profiles-N.pdf') ; print('=> Profiles-N.pdf saved')
	if 'h2o'  in Vars : figP.savefig(dirp+'Profiles-P.pdf') ; print('=> Profiles-P.pdf saved')
else :
	plt.show()