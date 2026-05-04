#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%=================================================================================
# from unittest import case
# from IPython import get_ipython

# ip = get_ipython()
# if ip is not None:
#     ip.run_line_magic("load_ext", "autoreload")
#     ip.run_line_magic("autoreload", "2")
#%%=================================================================================
#                     Modules
#===================================================================================
import Utilities as util
(Sysa,NSysa,Arg)=util.Parseur(['Save','Visu','Temp','Data','Over','Vel','Only','NoProf','Compa','TProf','Zoom','Struct','Report','Outlet','All'],1,'Arg : Case (Jaravel,Garnier,Sevault)')
(                             [ SAVE , VISU , TEMP , DATA , OVER , VEL , ONLY , NOPROF , COMPA , TPROF , ZOOM , STRUCT , REPORT , OUTLET , ALL ])=Arg ; Case=Sysa[0]
if   ALL and Case=='PRECIZE' : REPORT,VISU,NOPROF=True,True,True
elif ALL                     : SAVE,VISU,TEMP,OVER,DATA=True,True,True,True,True # Over : Overwrite
if VEL   : VISU,SAVE=True,True
if COMPA : SAVE,DATA=True,True

from numpy import *
from scipy.interpolate import interp1d
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
# cmesh=0 
cmesh=1e3
type='svg'

#====================> Burner
if   Case=='Jaravel' : from ParamsJaravel import * ; title=''
elif Case=='Garnier' : from ParamsGarnier import * ; title='he {:.0f} %'.format(he*100)
elif Case=='Sevault' : from ParamsSevault import * ; title='re {:.0f} k  ,  h2 {:.0f} %'.format(re,h2)
elif Case=='PRECIZE' : from ParamsPRECIZE import * ; title='PRECIZE'
else : sys.exit('=> Error : Case not recognized')

#====================> Fields
Vars=['T','o2','ch4','co2','h2o','co','h2']
# Vars=['Vel','k','T','o2','ch4','co2','h2o','co','h2']
# Vars=['Vel','k','T','o2','ch4','co2','h2o','co']
# Vars=['Vel','k','T','mixH','o2','h2','n2','h2o']
# Vars=['T']

#====================> Params
Titre={}
Titre['T'  ]='Temperature [K]'
Titre['Vel']='Velocity [m/s]'
Titre['k'  ]=r'k [$m^2/s^2$]'
Titre['co' ]=r'$Y_{CO}$ [%]'
Titre['o2' ]=r'$Y_{O_2}$ [-]'
Titre['n2' ]=r'$Y_{N_2}$ [-]'
Titre['h2' ]=r'$Y_{H_2}$ [-]'
Titre['ch4']=r'$Y_{CH_4}$ [-]'
Titre['co2']=r'$Y_{CO_2}$ [-]'
Titre['h2o']=r'$Y_{H_2O}$ [-]'

#====================> Visu domain
if   Case=='Jaravel' and dirc.split('/')[-2]=='DUMP-01-EDM' :
	fsize=(15,5)
	RY_d=RY_t=[ 0,0.014]
	RX_d=RX_t=[ 0,0.050]
elif Case=='Jaravel' and 'EDC' in dirc.split('/')[-2] :
	fsize=(25,4)
	RY_d=RY_t=[  0,0.05]
	RX_d=RX_t=[ Lc,0.40]
elif Case=='Sevault' and ZOOM :
	fsize=(25,5)
	RY_d=RY_t=[ 0,0.08]
	RX_d=RX_t=[ 0,0.40]
elif Case=='Sevault' :
	fsize=(25,5)
	RY_d=RY_t=[  0,0.08]
	RX_d=RX_t=[  0,0.80]
else :
	fsize=(25,5)
	# RX_t=[Lc,0.50]
	RY_t=[ 0,0.05]
	RY_d=[ 0,0.04]
	RX_t=[ 0,0.40]
	RX_d=[ 0,0.2 ]

#====================> Velocity profile
n=15
Np=int(1e4)

#%%=================================================================================
#                     Process
#===================================================================================
util.Entete1(104,[dirc,title],'Profiles TNF')
#===================================================================================
util.MKDIR(dirp)
#===================================================================================
def List_Jaravel(dirs,flame) :
	# lim=8
	lim=1e5
	Files_s=os.popen('ls %s/pmCDEF/pm%s.stat/*.Yfav'%(dirs,flame)).read().splitlines()
	Files_v=os.popen('ls %s/TUD_LDV_DEF/TUD_LDV_%s.*'%(dirs,flame)).read().splitlines()
	Pos_s=[] ; F_s=[] ; A_s=[]
	Pos_v=[] ; F_v=[] ; A_v=[]
	#====================> Positions Scattering
	for f in Files_s :
		pos0=f.split('/')[-1].split('.')[0][1:] #; print(pos0)
		if   pos0[:2]=='CL'  : A_s.append(f)
		elif pos0    =='075' : F_s.append(f) ; Pos_s.append( 7.5 )
		elif int(pos0) <lim  : F_s.append(f) ; Pos_s.append( int(pos0) )
	#====================> Positions Velocimetry
	for f in Files_v :
		pos0=f.split('/')[-1].split('.')[-1] #; print(pos0)
		if   pos0=='exit'      : pass
		elif pos0=='axial'     : A_v.append(f)
		elif pos0[1:]=='075'   : F_v.append(f) ; Pos_v.append( 7.5 )
		elif int(pos0[1:])<lim : F_v.append(f) ; Pos_v.append( int(pos0[1:]) )
	return(F_s+A_s,Pos_s,F_v+A_v,Pos_v)
#-----------------------------------------------------------------------------------
def List_Garnier(dirs,he) :
	#====================> Positions Scattering
	Files_s=os.popen('ls %s/he%.0f*.fav'%(dirs,he*10)).read().splitlines()
	Pos_s=[]
	for f in Files_s :
		pos0=f.split('/')[-1].split('.')[0]
		if 'x' in pos0 :
			pos=pos0.split('x')[1]
			xl=float(pos[0])/float(pos[1]) ; Pos_s.append(xl)
			print('=> Position  : x/L_vis={}'.format(xl))
	#====================> Positions Velocimetry
	Files_v=os.popen('ls %s/s*%.0f.dat'%(dirv,he*100)).read().splitlines()
	# Pos_v=[0,1/16,1/8,1/4,3/8,1/2,5/8,3/4,1]
	Pos_v=[0,1,1/2,1/4,1/16,1/8,3/4,3/8,5/8]
	return(Files_s,Pos_s,Files_v,Pos_v)
#-----------------------------------------------------------------------------------
def List_Sevault(dirs,f) :
	#====================> Positions Scattering
	Files_s=os.popen('ls %s/%s_*_mean_favre.txt'%(dirs,f)).read().splitlines()
	Pos_s=[]
	for f in Files_s :
		pos0=f.split('/')[-1].split('_')[1] ; print(pos0)
		Pos_s.append( int(pos0) )
	return(Files_s,Pos_s,[],[])
#===================================================================================
util.Section( 'TNF reading : {:.3f} s'.format(time.time()-t0),0,3,'b' )
#===================================================================================
if   Case=='Jaravel' : (Files_s,Pos_s,Files_v,Pos_v)=List_Jaravel(dirs,flame)
elif Case=='Garnier' : (Files_s,Pos_s,Files_v,Pos_v)=List_Garnier(dirs,he)
elif Case=='Sevault' : (Files_s,Pos_s,Files_v,Pos_v)=List_Sevault(dirs,flame)
#===================================================================================
util.Section( 'Inlet profile : {:.3f} s'.format(time.time()-t0),0,3,'b' )
#===================================================================================
if Case in ['Jaravel','Garnier','Sevault'] :
	Umax=Umoy*(n+2)/n ; print('=> n = %.0f  ,  Umoy = %.3f [m/s]  ,  Umax = %.3f [m/s]'%(n,Umoy,Umax))
	VyV=linspace(0,0.5*D0,Np)
	Uth=Umax*(1-(2*VyV/D0)**n)
	Npos_s=len(Pos_s) ; IPos_s=argsort(Pos_s)
	Npos_v=len(Pos_v) ; IPos_v=argsort(Pos_v)
#===================================================================================
if REPORT :
	util.Section( 'Report reading : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	# for rf in os.popen('ls %s/report-*-rfile.out'%(dirc)).read().splitlines() :
	for rf in os.popen('ls %s/report-*-rfile.out'%(dirc)).read().splitlines() :
		r_name=rf.split('/')[-1][7:-10]
		Dr=fl.Report_read(rf) ; Keys=list(Dr.keys()) #; print(Keys)
		if '(' in Keys[1] : Labels=['Iteration']+[ k.split('(')[1][:-1] for k in Keys[1:] if '(' in k ]
		else              : Labels=Keys
		print('=> \033[31m%s\033[0m : '%(r_name) , Labels )
		figr,axr=plt.subplots(figsize=(10,7)) #; bxr=axr.twinx()
		if r_name == 'mass-balance' :
			# (Mf_f,Mf_o,Mf_b,Mf_s,Mb)=fl.Mf_sep(Dr,Keys) ; print('=> fuel : %.3f g/s  ,  oxid : %.3f g/s  ,  out : %.3f g/s  ,  slope : %.3f g/s  ,  balance : %.3f g/s'%(mean(Mf_f)*1e3,mean(Mf_o)*1e3,mean(Mf_b)*1e3,mean(Mf_s)*1e3,mean(Mb)*1e3))
			(Mf_f,Mf_o,Mf_b,Mf_l,Mf_s,Mb)=fl.Mf_sep2(Dr,Keys)
			print('=> fuel : %.3f g/s  ,  oxid : %.3f g/s  ,  out : %.3f g/s  ,   leak : %.3f g/s  ,  slope : %.3f g/s  ,  balance : %.3f g/s'%(mean(Mf_f[-1])*1e3,mean(Mf_o[-1])*1e3,mean(Mf_b[-1])*1e3,mean(Mf_f[-1])*1e3,mean(Mf_s[-1])*1e3,mean(Mb[-1])*1e3))
			axr.plot( Dr['Iteration'],1e3*Mf_f  ,label='inlet-fuel'  )
			axr.plot( Dr['Iteration'],1e3*Mf_o  ,label='inlet-oxid'  )
			axr.plot( Dr['Iteration'],1e3*Mf_b  ,label='outlet'      )
			axr.plot( Dr['Iteration'],1e3*Mf_l  ,label='leak'        )			
			axr.plot( Dr['Iteration'],1e3*Mf_s  ,label='slope-zone'  )			
			# bxr.plot( Dr['Iteration'],1e3*Mb,'k',label='mass-balance')
			axr.plot( Dr['Iteration'],1e3*Mb,'k',label='mass-balance')
			axr.set_ylabel('Mass flow rate [g/s]',fontsize=25)
			figr.legend(fontsize=15,loc='center',bbox_to_anchor=(0.7,0.3))
		else :
			for n,k in enumerate(Keys[1:]) : axr.plot( Dr['Iteration'],Dr[k],label=Labels[n+1] )
			if len(Keys)>2 : axr.legend(fontsize=15) #,loc='center',bbox_to_anchor=(0.7,0.3))
		# if len(Keys)>2 : figr.legend(fontsize=15,loc='center',bbox_to_anchor=(0.7,0.3))
		if r_name in ['heat-release','mass-balance','shell-loss','co2'] :
			axr.ticklabel_format( axis='y' , scilimits=(-3,3) )
		util.SaveFig(figr,dirp+'Report-%s.pdf'%(r_name))
	if ONLY : sys.exit('=> Stop after report reading')
#===================================================================================
if TEMP :
	util.Section( 'Temporals : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	for v in ['F','T','U'] : 
		f_temp=dirc+f'probe-{v:s}.out'
		if os.path.exists(f_temp) : fl.Probe_plot(f_temp,dirp+f'Probe-{v:s}.pdf')
	(D_in,D_ou)=fl.Mass_read(dirc+'InletMassFlowRate.out',dirc+'outletMassFlowRate.out')
	It =D_in[:,0]
	Min=D_in[:,1]+D_in[:,2]
	Mou=D_ou[:,1]+D_ou[:,2]
	fig_b,ax_b=fl.Report(It,Min,Mou,0,0,1e9,int(1e2),1e3)
	fig_b.tight_layout()
	fig_b.savefig(dirp+'Balance.pdf')
	if ONLY : sys.exit('=> Stop after temporals')
#%%=================================================================================
if TPROF :
	h=16.5
	l=16.2
	T0=300
	T1=400
	L=1e-2
	om=sqrt(h/l)
	Xw=linspace(0,L,int(1e3))
	Tw=T0+(T1-T0)*exp( (Xw-L)*om ) ; print('=> Tw0 : %.2f'%(Tw[0]))
	figw,axw=plt.subplots()
	axw.plot(Xw,Tw,'k')
	util.SaveFig(figw,dirp+'TProf_Wall.pdf')
	if ONLY : sys.exit('=> Stop after wall profile')
#%%=================================================================================
if STRUCT :
	Tiso=list(array([1600,1800])-0)
	alp=BC_f['CH4']/BC_f['H2']
	XH2_s =1/(1.5+3*alp)
	XCH4_s=       alp *XH2_s
	XO2_s =(0.5+2*alp)*XH2_s
	Mav_s=fl.Mol_m['H2']*XH2_s+fl.Mol_m['O2']*XO2_s+fl.Mol_m['CH4']*XCH4_s
	YH_s =fl.Mol_m['H']*(4*XCH4_s+2*XH2_s)/Mav_s
	if BC_f['unit']=='X' : BC_f=fl.ConvBC_XY(BC_f)
	if BC_o['unit']=='X' : BC_o=fl.ConvBC_XY(BC_o)
	YH_o=fl.Yh(BC_o,fl.Mol_m)
	YH_f=fl.Yh(BC_f,fl.Mol_m)
	Zst=(YH_s-YH_o)/(YH_f-YH_o) ; print('=> Stoichiometric mixture fraction Zst = %.3f'%(Zst))
	# fl.Visu(dird+slice,'velocity-magnitude' ,r'Velocity [m/s]'        ,[0,0.3],[0,0.05],arange(0,350,50)    ,cmesh,[],(25,5),'cividis',dirp+'Struct-V.png',[])
	# fl.Visu(dird+slice,'turb-kinetic-energy',r'k [$m^2/s^2$]'         ,[0,0.3],[0,0.05],arange(250,2500,250),cmesh,[],(25,5),'cividis',dirp+'Struct-K.png',[])
	# fl.Visu(dird+slice,'turb-diss-rate'     ,r'$\epsilon$ [$m^2/s^3$]',[0,0.3],[0,0.05],arange(0,2e7,1e6)   ,cmesh,[],(25,5),'cividis',dirp+'Struct-E.png',[])
	fl.Visu(dird+slice,'tc'                 ,r'$1/t_c$ [$1/s$]'       ,[0,0.3],[0,0.05],logspace(1,4,5)     ,cmesh,[],(25,5),'cividis',dirp+'Struct-tc.png',['Zst',[Zst,'b'],'MIXH',[fl.Mol_m,BC_f,BC_o],'Tiso',[1600,1800,2000]] )
	# fl.Visu(dird+slice,'tt'                 ,r'$1/t_t$ [$1/s$]'       ,[0,0.3],[0,0.05],logspace(1,4,5)     ,cmesh,[],(25,5),'cividis',dirp+'Struct-tt.png',['Zst',[Zst,'b'],'MIXH',[fl.Mol_m,BC_f,BC_o],'Tiso',[1600,1800,2000],'tt',['dyn',4]] )
	# fl.Visu(dird+slice,'mixH'               ,r'Mixture fraction [-]'  ,[0,0.3],[0,0.05],arange(0,1.1,0.1)   ,cmesh,[],(25,5),'viridis',dirp+'Struct-Z.png',['ISO',[Zst],'MIXH',[fl.Mol_m,BC_f,BC_o]])
	# fl.Visu(dird+slice,'temperature'        ,r'Temperature [K]'       ,[0,0.3],[0,0.05],arange(250,2500,250),cmesh,[],(25,5),'inferno',dirp+'Struct-T.png',['ISO',[2000],'RECIRC','b','Zst',[Zst,'b'],'MIXH',[fl.Mol_m,BC_f,BC_o]])
	if ONLY : sys.exit('=> Stop after mixture fraction')
#%%=================================================================================
if VISU :
	slice_m=slice_f=slice
	BD_Vars={}
	Ticks={}
	Param={'co':[] ,'o2':[] ,'h2':[] ,'ch4':[],'co2':[],'h2o':[],'T':[],'Vel':[]}
	util.Section( 'Visualisation : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	name0='Visu-'
	#=========================> Positions
	if Case=='Garnier' :
		Vp_v=[ Lc+n*D0*Ld for n in Pos_v ]
		Vp_s=[ Lc+n*D0*Ld for n in Pos_s ]
		Vp_h=[]
	elif Case=='Sevault' :
		Vp_v=[ Lc+z*1e-3       for z in Pos_s ]
		Vp_s=[ Lc+z*1e-3       for z in Pos_s ]
		Vp_h=[]
	elif Case=='Jaravel' :
		Vp_v=[ Lc+n*D0 for n in Pos_v ]
		Vp_s=[ Lc+n*D0 for n in Pos_s ]
		Vp_h=linspace(Lc,80*D0,1000)
	elif Case=='PRECIZE' :
		Vp_v=[]
		Vp_s=[]
	#=========================> Struct
	if STRUCT : Param['T']=['MIX',[fl.Mol_m,BC_f,BC_o]]
	else      : Param['T']=[]
	#=========================> Params
	if   Case in ['Garnier','Sevault','Jaravel'] : 
		Tiso=list(array([1600,1800])-0)
		Param['T']+=['RECIRC','b']
		Param['co']=['CO',1e2]
		BD_Vars['T' ]=[]
		BD_Vars['co']=[]
		BD_Vars['o2']=[]
		BD_Vars['co2']=[]
		BD_Vars['h2o']=[]
		BD_Vars['k'  ]=[]
		Ticks['T'  ]=arange(250,3001,250)
		Ticks['co' ]=arange(0,1.1,0.1)
		Ticks['o2' ]=arange(0,1.1,0.1)
		Ticks['h2' ]=arange(0,1.1,0.1)
		Ticks['ch4']=arange(0,1.1,0.1)
		Ticks['co2']=arange(0,1.1,0.1)
		Ticks['h2o']=arange(0,1.1,0.1)
		Ticks['k'  ]=arange(250,2500,250)
		slice_m=slice_f=slice
		if   Case=='Jaravel' :
			Ticks['o2' ]=arange(0,0.25,0.02)
			Ticks['ch4']=arange(0,0.2 ,0.02)
			Ticks['k'  ]=arange(0,100,10)
			Ticks['co' ]=arange(0,10,1)
			BD_Vars['k'  ]=[Ticks['k'][0],Ticks['k'][-1]]
		elif Case=='Sevault' :
			Ticks['co' ]=arange(0,50,2)
	if ZOOM : name0+='Zoom-'
	F_int={}
	if 'T'    in Vars : F_int['T'   ]=fl.Visu(dird+slice_m,'xy','temperature'        ,Titre['T' ]           ,RX_t,RY_t,Ticks['T']          ,cmesh,BD_Vars['T']  ,fsize,'inferno',dirp+name0+'Temperature.png',['INTERP','LINES',Vp_s,'ISO',Tiso]+Param['T']  )
	if 'Vel'  in Vars : F_int['Vel' ]=fl.Visu(dird+slice_f,'xy','velocity-magnitude' ,Titre['Vel']          ,RX_d,RY_d,arange(0,350,50)    ,cmesh,[]            ,fsize,'cividis',dirp+name0+'Velocity.png'   ,['INTERP','LINES',Vp_v]           +Param['Vel'])
	if 'k'    in Vars : F_int['k'   ]=fl.Visu(dird+slice_f,'xy','turb-kinetic-energy',Titre['k'  ]          ,RX_d,RY_d,Ticks['k']          ,cmesh,BD_Vars['k']  ,fsize,'cividis',dirp+name0+'TKE.png'        ,['INTERP'])
	if 'tt'   in Vars : F_int['tt'  ]=fl.Visu(dird+slice_f,'xy','tt'                 ,'tt [s]'              ,RX_d,RY_d,arange(0,1e-3,1e-4) ,cmesh,[]            ,fsize,'cividis',dirp+name0+'tt.png'         ,['INTERP'])
	if 'mixH' in Vars : F_int['mixH']=fl.Visu(dird+slice_f,'xy','mixH'               ,'Mixture fraction [-]',RX_t,RY_t,arange(0,1.1,0.1  ) ,cmesh,[]            ,fsize,'viridis',dirp+name0+'Mix.png'        ,['INTERP','MIXH',[fl.Mol_m,BC_f,BC_o]])
	if 'mix'  in Vars : F_int['mix' ]=fl.Visu(dird+slice_f,'xy','fmean'              ,'Mixture fraction [-]',RX_t,RY_t,arange(0,1.1,0.1  ) ,cmesh,[]            ,fsize,'viridis',dirp+name0+'Fmean.png'      ,['INTERP'])
	if 'h2'   in Vars : F_int['h2'  ]=fl.Visu(dird+slice_f,'xy','h2'                 ,Titre['h2']           ,RX_t,RY_t,Ticks['h2']         ,cmesh,[]            ,fsize,'viridis',dirp+name0+'H2.png'         ,['INTERP']+Param['h2'])
	if 'n2'   in Vars : F_int['n2'  ]=fl.Visu(dird+slice_f,'xy','n2'                 ,Titre['n2']           ,RX_t,RY_t,arange(0,1.1,0.1  ) ,cmesh,[]            ,fsize,'viridis',dirp+name0+'N2.png'         ,['INTERP'])
	if 'o2'   in Vars : F_int['o2'  ]=fl.Visu(dird+slice_f,'xy','o2'                 ,Titre['o2']           ,RX_t,RY_t,Ticks['o2']         ,cmesh,BD_Vars['o2'] ,fsize,'viridis',dirp+name0+'O2.png'         ,['INTERP']+Param['o2'])
	if 'co'   in Vars : F_int['co'  ]=fl.Visu(dird+slice_f,'xy','co'                 ,Titre['co']           ,RX_t,RY_t,Ticks['co']         ,cmesh,BD_Vars['co'] ,fsize,'viridis',dirp+name0+'CO.png'         ,['INTERP']+Param['co'])
	if 'ch4'  in Vars : F_int['ch4' ]=fl.Visu(dird+slice_f,'xy','ch4'                ,Titre['ch4']          ,RX_t,RY_t,Ticks['ch4']        ,cmesh,[]            ,fsize,'viridis',dirp+name0+'CH4.png'        ,['INTERP']+Param['ch4'])
	if 'h2o'  in Vars : F_int['h2o' ]=fl.Visu(dird+slice_f,'xy','h2o'                ,Titre['h2o']          ,RX_t,RY_t,Ticks['h2o']        ,cmesh,BD_Vars['h2o'],fsize,'viridis',dirp+name0+'H2O.png'        ,['INTERP']+Param['h2o'])
	if 'co2'  in Vars : F_int['co2' ]=fl.Visu(dird+slice_f,'xy','co2'                ,Titre['co2']          ,RX_t,RY_t,Ticks['co2']        ,cmesh,BD_Vars['co2'],fsize,'viridis',dirp+name0+'CO2.png'        ,['INTERP']+Param['co2'])
	if 'Xo2'  in Vars : F_int['Xo2' ]=fl.Visu(dird+slice_f,'xy','molef-o2'           ,Titre['Xo2']          ,RX_t,RY_t,Ticks['o2']         ,cmesh,BD_Vars['o2'] ,fsize,'viridis',dirp+name0+'XO2.png'        ,['INTERP']+Param['o2'])
	if 'Xco'  in Vars : F_int['Xco' ]=fl.Visu(dird+slice_f,'xy','molef-co'           ,Titre['Xco']          ,RX_t,RY_t,Ticks['co']         ,cmesh,BD_Vars['co'] ,fsize,'viridis',dirp+name0+'XCO.png'        ,['INTERP']+Param['co'])
	if 'Xco2' in Vars : F_int['Xco2']=fl.Visu(dird+slice_f,'xy','molef-co2'          ,Titre['Xco2']         ,RX_t,RY_t,Ticks['co2']        ,cmesh,BD_Vars['co2'],fsize,'viridis',dirp+name0+'XCO2.png'       ,['INTERP']+Param['co2'])
	if 'Xh2o' in Vars : F_int['Xh2o']=fl.Visu(dird+slice_f,'xy','molef-h2o'          ,Titre['Xh2o']         ,RX_t,RY_t,Ticks['h2o']        ,cmesh,BD_Vars['h2o'],fsize,'viridis',dirp+name0+'XH2O.png'       ,['INTERP']+Param['h2o'])

	if not NOPROF :
		D_int={}
		Vy=linspace(0,0.5*D2,Np)
		Vx=linspace(0,Lc+Lt,Np)
		for v in Vars :
			if v in ['Vel','k'] : D_int[v]=[ F_int[v](Np*[p],Vy) for p in Vp_v ]+[ F_int[v](Vx,Np*[0]) ]
			else                : D_int[v]=[ F_int[v](Np*[p],Vy) for p in Vp_s ]+[ F_int[v](Vx,Np*[0]) ]
			# print( len(D_int[v][-1]) )
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
				file.closed

	if VEL and 'Vel' in Vars and 'k' in Vars :
		util.Section( 'Velocity profiles : {:.3f} s'.format(time.time()-t0),1,5,'r' )
		PosV=array([0.01,0.1,0.5,0.9,1])*Lc
		# PosV=array([0])*Lc
		figVp,axVp=plt.subplots(figsize=(10,7)) ; figVp.suptitle('Velocity profiles'       ,fontsize=30)
		figKp,axKp=plt.subplots(figsize=(10,7)) ; figKp.suptitle('Turbulent kinetic energy',fontsize=30)
		axVp.set_xlabel('r [mm]'      ,fontsize=30) ; axKp.set_xlabel('r [mm]'       ,fontsize=30)
		axVp.set_ylabel('V [m/s]'     ,fontsize=30) ; axKp.set_ylabel('k [$m^2/s^2$]',fontsize=30)
		print( '=> Poiseuil  ,  Max V=%.3f [m/s]  ,  Mean V=%.3f [m/s]'%(max(Uth),fl.Umoy(VyV,Uth)) ) 
		axVp.plot( VyV*1e3,Uth , 'k--' , label='Theoretical' )
		MVprof=zeros((Np,len(PosV)))
		MKprof=zeros((Np,len(PosV)))
		for n,p in enumerate(PosV) :
			Vprof=F_int['Vel'](Np*[p],VyV) ; MVprof[:,n]=Vprof
			Kprof=F_int['k'  ](Np*[p],VyV) ; MKprof[:,n]=Kprof
			print('=> x/Lx=%.2f  ,  Max V=%.3f [m/s]  ,  Mean V=%.3f [m/s]'%(p/Lc,max(Vprof),fl.Umoy(VyV,Vprof)) ) 
			axVp.plot( VyV*1e3,Vprof , label='x/Lc=%.2f'%(p/Lc) )
			axKp.plot( VyV*1e3,Kprof , label='x/Lc=%.2f'%(p/Lc) )
		axVp.legend(fontsize=20) ; util.SaveFig(figVp,dirp+'Inlet-Velocity.pdf')
		axKp.legend(fontsize=20) ; util.SaveFig(figKp,dirp+'Inlet-TKE.pdf'     )
		if DATA :
			St_data=dirp+'Inlet-Profiles-Velocity.csv'
			with open(St_data,'w') as file :
				writer=csv.writer(file)
				writer.writerow( ['r(mm)'] + [ 'x/Lc=%.2f'%(p/Lc) for p in PosV ] )
				writer.writerows( column_stack( ( VyV*1e3 , MVprof ) ) )
			file.closed
			St_data=dirp+'Inlet-Profiles-TKE.csv'
			with open(St_data,'w') as file :
				writer=csv.writer(file)
				writer.writerow( ['r(mm)'] + [ 'x/Lc=%.2f'%(p/Lc) for p in PosV ] )
				writer.writerows( column_stack( ( VyV*1e3 , MKprof ) ) )
			file.closed
else : 
	Vy=[]
	D_int={}
if NOPROF and not COMPA : sys.exit('=> Stop after visu')
#%%=================================================================================
#                     Profiles
#===================================================================================
def ReadTNF(name,Params) :
	# print(name)
	sep,skip,i0=Params[:3]
	with open(name) as file :
		Lines=file.readlines()
		if i0>=0 : T=[ s.strip() for s in Lines[skip][i0:-1].split(sep) if s ]
		else     : T=Params[3]
		M=array([ [float(v) for v in L.split(sep) if v ] for L in Lines[skip+1:] ]) ; (np,nv)=M.shape
	file.closed
	return({ s:M[:,i] for i,s in enumerate(T[:nv]) })
#=====================================================================================
def tpos_garnier(pos) : return('z/L_vis=%.2f'%(pos))
def tpos_sevault(pos) : return('z(mm)=%.0f'%(pos))
def tpos_jaravel(pos) : 
	if pos==7.5 : return('z/d=%.1f'%(pos))
	else        : return('z/d=%.0f'%(pos))
#=====================================================================================
def Var_dyn(Data,Case) :
	if   Case=='Garnier' :
		Vel=hypot( Data['u']    , Data['v']    )
		Tke=hypot( Data['varu'] , Data['varv'] )
	elif Case=='Jaravel' :
		Vel=hypot( Data['u']    , Data['v']    )
		Tke=       Data['varu'] + Data['varv']
	return(Vel,Tke)
#=====================================================================================
def plot_exp( ax , X,Y,E ) : ax.errorbar(X,Y,yerr=E*abs(Y),ecolor='k',color='k',marker='o') #,linestyle='none')
def Plot_Exp( ax , Data,var , Cor,Err,titre,xlim) :
	X=Data[Cor['r']]*(D0*1e3)**(Case=='Jaravel')
	Y=Data[Cor[var]]*(1e2**(var=='co'))
	if len(Data.keys()) > 0 : plot_exp( ax , X,Y,Err[var] )
	ax.set_title(titre,fontsize=20)
	ax.ticklabel_format(style='sci',axis='y',scilimits=(-2,4))
	if xlim[1]>0 : ax.set_xlim(xlim)
#=====================================================================================
def PlotFile(ax,dir0,var,nc) :
	D=util.ReadCSV(dir0+'Profiles-%s.csv'%(var)) ; Keys=list(D.keys()) #; print('D keys : ',Keys)
	Rad=D[Keys[0]]
	for n,k in enumerate(Keys[1:]) :
		i= n//nc
		j= n-nc*i
		ax[i,j].plot(Rad,D[k])
	ax[-1,-1].plot(1,1,'.',label=dir0.split('/')[-2][5:])
	ax[-1,-1].legend(loc='center',bbox_to_anchor=(0.5,0.5))
#=====================================================================================
def Plot_Axial( ax , X_sim,Y_sim , X_exp,Y_exp,E ) :
	ax.plot( X_sim , Y_sim , 'r'  )
	plot_exp( ax , X_exp , Y_exp , E )
#=====================================================================================
def Profile(ax,Vy,D_int,Data,Cor,Err,p,var,tpos,xlim) :
	if len(Vy)>0 : grid=Vy*1e3
	# if var=='co' : coef=1e2
	# else         : coef=1
	coef=1
	if var in D_int.keys()  : ax.plot(grid,coef*D_int[var][p],'r')
	Plot_Exp( ax , Data,var , Cor,Err,tpos(Pos_s[p]),xlim )
	# if len(Data.keys()) > 0 : ax.errorbar(Data[Cor['r']]*D0*1e3,coef*Data[Cor[var]],yerr=Err[var]*abs(Data[Cor[var]]),ecolor='k',color='k',marker='.') #,linestyle='none')
	# ax.set_title(tpos(Pos_s[p]),fontsize=20)
	# ax.ticklabel_format(style='sci',axis='y',scilimits=(-2,4))
	# if xlim[1]>0 : ax.set_xlim(xlim)
#=====================================================================================
def Saving(dirp,Vars,type) :
	if 'Vel'  in Vars : util.SaveFig(figV,dirp+'Profiles-V.'+type)
	if 'k'    in Vars : util.SaveFig(figK,dirp+'Profiles-K.'+type)
	if 'T'    in Vars : util.SaveFig(figT,dirp+'Profiles-T.'+type)
	if 'mixH' in Vars : util.SaveFig(figZ,dirp+'Profiles-Z.'+type)
	if 'o2'   in Vars : util.SaveFig(figO,dirp+'Profiles-O.'+type)
	if 'h2'   in Vars : util.SaveFig(figH,dirp+'Profiles-H.'+type)
	if 'n2'   in Vars : util.SaveFig(figN,dirp+'Profiles-N.'+type)
	if 'co'   in Vars : util.SaveFig(figA,dirp+'Profiles-A.'+type)
	if 'h2o'  in Vars : util.SaveFig(figP,dirp+'Profiles-P.'+type)
	if 'co2'  in Vars : util.SaveFig(figC,dirp+'Profiles-C.'+type)
	if 'ch4'  in Vars : util.SaveFig(figM,dirp+'Profiles-M.'+type)
#%%===================================================================================
#                     Figures
#=====================================================================================
if   Case=='Garnier' : Params_s=( ' ',3, 0 )         ; tpos=tpos_garnier ; Rlims=Npos_s*[0]                   ; Abs=['y'  ,1     ] ; Params_v=('\t' ,11,2)
elif Case=='Jaravel' : Params_s=('  ',3, 0 )         ; tpos=tpos_jaravel ; Rlims=[14,14,16,16,20,40,50,70,80] ; Abs=['X/D',D0*1e3] ; Params_v=(4*' ',12,-1,['X/D','u','varu','v','varv'])
elif Case=='Sevault' : Params_s=('\t',0,-1,Vars_TNF) ; tpos=tpos_sevault ; Rlims=[5,8,11,15,20]
#=====================================================================================
nc=2 ; nr=Npos_s//nc+int(Npos_s%nc>0)
if Case=='Garnier' : Nv=[3,3]
else               : Nv=[2,0] ; Nv[1]=Npos_v//Nv[0]+int(Npos_v%Nv[0]>0)
Nv[1]+=(Nv[1]==0)

if 'Vel'  in Vars : figV,axV=plt.subplots(ncols=Nv[0],nrows=Nv[1],figsize=(13,10)) ; figV.suptitle('Velocity [m/s]'      ,fontsize=30)
if 'k'    in Vars : figK,axK=plt.subplots(ncols=Nv[0],nrows=Nv[1],figsize=(13,10)) ; figK.suptitle('TKE [$m^2/s^2$]'     ,fontsize=30)
if 'T'    in Vars : figT,axT=plt.subplots(ncols=nc   ,nrows=nr   ,figsize=(13,10)) ; figT.suptitle('Temperature [K]'     ,fontsize=30)
if 'mixH' in Vars : figZ,axZ=plt.subplots(ncols=nc   ,nrows=nr   ,figsize=(13,10)) ; figZ.suptitle('Mixture fraction [-]',fontsize=30)
if 'o2'   in Vars : figO,axO=plt.subplots(ncols=nc   ,nrows=nr   ,figsize=(13,10)) ; figO.suptitle(r'$Y_{O_2}$ [-]'      ,fontsize=30)
if 'h2'   in Vars : figH,axH=plt.subplots(ncols=nc   ,nrows=nr   ,figsize=(13,10)) ; figH.suptitle(r'$Y_{H_2}$ [-]'      ,fontsize=30)
if 'n2'   in Vars : figN,axN=plt.subplots(ncols=nc   ,nrows=nr   ,figsize=(13,10)) ; figN.suptitle(r'$Y_{N_2}$ [-]'      ,fontsize=30)
if 'co'   in Vars : figA,axA=plt.subplots(ncols=nc   ,nrows=nr   ,figsize=(13,10)) ; figA.suptitle(r'$Y_{CO}$ [-]'       ,fontsize=30)
if 'h2o'  in Vars : figP,axP=plt.subplots(ncols=nc   ,nrows=nr   ,figsize=(13,10)) ; figP.suptitle(r'$Y_{H_2O}$ [-]'     ,fontsize=30)
if 'co2'  in Vars : figC,axC=plt.subplots(ncols=nc   ,nrows=nr   ,figsize=(13,10)) ; figC.suptitle(r'$Y_{CO_2}$ [-]'     ,fontsize=30)
if 'ch4'  in Vars : figM,axM=plt.subplots(ncols=nc   ,nrows=nr   ,figsize=(13,10)) ; figM.suptitle(r'$Y_{CH_4}$ [-]'     ,fontsize=30)
if Npos_s%nc>0 : #and not COMPA :
	if 'T'    in Vars : axT[-1,-1].axis('off')
	if 'mixH' in Vars : axZ[-1,-1].axis('off')
	if 'o2'   in Vars : axO[-1,-1].axis('off')
	if 'h2'   in Vars : axH[-1,-1].axis('off')
	if 'n2'   in Vars : axN[-1,-1].axis('off')
	if 'co'   in Vars : axA[-1,-1].axis('off')
	if 'h2o'  in Vars : axP[-1,-1].axis('off')
	if 'co2'  in Vars : axC[-1,-1].axis('off')
	if 'ch4'  in Vars : axM[-1,-1].axis('off')
#=====================================================================================
if COMPA :
	util.MKDIR(dirp)
	D_compa=[ d for d in os.listdir(dirc) if d[:5]=='DATA-' ]
	D_compa.sort()
	#=====> Axial
	Data_v=ReadTNF(Files_v[-1],Params_v) ; (Vel,Tke)=Var_dyn(Data_v,Case)
	Data_s=ReadTNF(Files_s[-1],Params_s)
	FIG_a,AX_a=[],[]
	for v in Vars :
		fig_a,ax_a=plt.subplots(figsize=(10,7)) ; FIG_a.append(fig_a) ; AX_a.append(ax_a)
		ax_a.set_ylabel(Titre[v]             ,fontsize=30)
		ax_a.set_xlabel('Axial position [mm]',fontsize=30)
		if   v=='Vel' : plot_exp( ax_a , Data_v[Abs[0]]*Abs[1]  ,Vel                , Err[v] )
		elif v=='k'   : plot_exp( ax_a , Data_v[Abs[0]]*Abs[1]  ,Tke                , Err[v] )
		elif v=='co'  : plot_exp( ax_a , Data_s[Cor['r']]*D0*1e3,Data_s[Cor[v]]*1e2 , Err[v] )
		else          : plot_exp( ax_a , Data_s[Cor['r']]*D0*1e3,Data_s[Cor[v]]     , Err[v] )
	for d in D_compa :
		dc=dirc+d+'/' #; print('=> Compa with : %s '%(d))
		D=util.ReadCSV(dc+'Axial-Profiles.csv') ; Keys=list(D.keys()) #; print('D keys : ',Keys)
		Rad=D[Keys[0]]
		for n,v in enumerate(Vars) :
			AX_a[n].plot(Rad,D[v],label=d[5:])
	for n,v in enumerate(Vars) :
		AX_a[n].legend(fontsize=20)
		util.SaveFig(FIG_a[n],dirp+'Axial-%s.'%(v)+type)
	#=====> Expe
	for n,p in enumerate(IPos_s) :
		i= n//nc
		j= n-nc*i
		Data=ReadTNF(Files_s[p],Params_s)
		if 'T'    in Vars : Plot_Exp(axT[i,j],Data,'T'   ,Cor,Err,tpos(Pos_s[p]),[0,Rlims[p]])
		if 'mixH' in Vars : Plot_Exp(axZ[i,j],Data,'mixH',Cor,Err,tpos(Pos_s[p]),[0,Rlims[p]])
		if 'o2'   in Vars : Plot_Exp(axO[i,j],Data,'o2'  ,Cor,Err,tpos(Pos_s[p]),[0,Rlims[p]])
		if 'h2'   in Vars : Plot_Exp(axH[i,j],Data,'h2'  ,Cor,Err,tpos(Pos_s[p]),[0,Rlims[p]])
		if 'n2'   in Vars : Plot_Exp(axN[i,j],Data,'n2'  ,Cor,Err,tpos(Pos_s[p]),[0,Rlims[p]])
		if 'co'   in Vars : Plot_Exp(axA[i,j],Data,'co'  ,Cor,Err,tpos(Pos_s[p]),[0,Rlims[p]])
		if 'h2o'  in Vars : Plot_Exp(axP[i,j],Data,'h2o' ,Cor,Err,tpos(Pos_s[p]),[0,Rlims[p]])
		if 'co2'  in Vars : Plot_Exp(axC[i,j],Data,'co2' ,Cor,Err,tpos(Pos_s[p]),[0,Rlims[p]])
		if 'ch4'  in Vars : Plot_Exp(axM[i,j],Data,'ch4' ,Cor,Err,tpos(Pos_s[p]),[0,Rlims[p]])
	#=====> Simu
	for d in D_compa :
		dc=dirc+d+'/' #; print('=> Compa with : %s '%(d))
		if 'T'    in Vars : PlotFile(axT,dc,'T'   ,nc)
		if 'mixH' in Vars : PlotFile(axZ,dc,'mixH',nc)
		if 'o2'   in Vars : PlotFile(axO,dc,'o2'  ,nc)
		if 'h2'   in Vars : PlotFile(axH,dc,'h2'  ,nc)
		if 'n2'   in Vars : PlotFile(axN,dc,'n2'  ,nc)
		if 'co'   in Vars : PlotFile(axA,dc,'co'  ,nc)
		if 'h2o'  in Vars : PlotFile(axP,dc,'h2o' ,nc)
		if 'co2'  in Vars : PlotFile(axC,dc,'co2' ,nc)
		if 'ch4'  in Vars : PlotFile(axM,dc,'ch4' ,nc)
	#=====> Saving
	if NOPROF : Saving(dirp,Vars,type) ; sys.exit('=> Stop after compa')
#=====================================================================================

#=====> Velocity profiles
if Case in ['Garnier','Jaravel'] :
	for n,p in enumerate(IPos_v) :
		i= n//Nv[0]
		j= n-Nv[0]*i
		Data=ReadTNF(Files_v[p],Params_v)
		(Vel,Tke)=Var_dyn(Data,Case)
		if 'Vel' in Vars : axV[i,j].plot( Vy*1e3,D_int['Vel' ][p],'r' ) ; axV[i,j].plot( Data[Abs[0]]*Abs[1],Vel,'-ok' ) ; axV[i,j].set_title(tpos(Pos_v[n]),fontsize=20) ; axV[i,j].ticklabel_format(style='sci',axis='y',scilimits=(-2,4)) ; axV[i,j].set_xlim((0,Rlims[n]))
		if 'k'   in Vars : axK[i,j].plot( Vy*1e3,D_int['k'   ][p],'r' ) ; axK[i,j].plot( Data[Abs[0]]*Abs[1],Tke,'-ok' ) ; axK[i,j].set_title(tpos(Pos_v[n]),fontsize=20) ; axK[i,j].ticklabel_format(style='sci',axis='y',scilimits=(-2,4)) ; axK[i,j].set_xlim((0,Rlims[n]))
		if i==Nv[1]-1 : 
			if 'Vel' in Vars : axV[i,j].set_xlabel('r [mm]',fontsize=20)
			if 'k'   in Vars : axK[i,j].set_xlabel('r [mm]',fontsize=20)
	if 'Vel' in Vars : axV[0,0].plot( [0,Rlims[0]],2*[Umoy],':k' )

#=====> Scattering
for n,p in enumerate(IPos_s) :
	i= n//nc
	j= n-nc*i
	Data=ReadTNF(Files_s[p],Params_s)
	if 'T'    in Vars : Profile(axT[i,j],Vy,D_int,Data,Cor,Err,p,'T'   ,tpos,[0,Rlims[p]])
	if 'mixH' in Vars : Profile(axZ[i,j],Vy,D_int,Data,Cor,Err,p,'mixH',tpos,[0,Rlims[p]])
	if 'o2'   in Vars : Profile(axO[i,j],Vy,D_int,Data,Cor,Err,p,'o2'  ,tpos,[0,Rlims[p]])
	if 'h2'   in Vars : Profile(axH[i,j],Vy,D_int,Data,Cor,Err,p,'h2'  ,tpos,[0,Rlims[p]])
	if 'n2'   in Vars : Profile(axN[i,j],Vy,D_int,Data,Cor,Err,p,'n2'  ,tpos,[0,Rlims[p]])
	if 'co'   in Vars : Profile(axA[i,j],Vy,D_int,Data,Cor,Err,p,'co'  ,tpos,[0,Rlims[p]])
	if 'h2o'  in Vars : Profile(axP[i,j],Vy,D_int,Data,Cor,Err,p,'h2o' ,tpos,[0,Rlims[p]])
	if 'co2'  in Vars : Profile(axC[i,j],Vy,D_int,Data,Cor,Err,p,'co2' ,tpos,[0,Rlims[p]])
	if 'ch4'  in Vars : Profile(axM[i,j],Vy,D_int,Data,Cor,Err,p,'ch4' ,tpos,[0,Rlims[p]])
	if i==nr-1 : 
		if 'T'    in Vars : axT[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'mixH' in Vars : axZ[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'o2'   in Vars : axO[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'h2'   in Vars : axH[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'n2'   in Vars : axN[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'co'   in Vars : axA[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'h2o'  in Vars : axP[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'co2'  in Vars : axC[i,j].set_xlabel('r [mm]',fontsize=20)
		if 'ch4'  in Vars : axM[i,j].set_xlabel('r [mm]',fontsize=20)
	if j==0 :
		if 'T'    in Vars : axT[i,j].set_ylabel(r'$T~[K]$'       ,fontsize=20)
		if 'mixH' in Vars : axZ[i,j].set_ylabel(r'$Z~[-]$'       ,fontsize=20)
		if 'o2'   in Vars : axO[i,j].set_ylabel(r'$Y_{O_2}~[-]$' ,fontsize=20)
		if 'h2'   in Vars : axH[i,j].set_ylabel(r'$Y_{H_2}~[-]$' ,fontsize=20)
		if 'n2'   in Vars : axN[i,j].set_ylabel(r'$Y_{N_2}~[-]$' ,fontsize=20)
		if 'co'   in Vars : axA[i,j].set_ylabel(r'$Y_{CO}~[\%]$' ,fontsize=20)
		if 'h2o'  in Vars : axP[i,j].set_ylabel(r'$Y_{H_2O}~[-]$',fontsize=20)
		if 'co2'  in Vars : axC[i,j].set_ylabel(r'$Y_{CO_2}~[-]$',fontsize=20)
		if 'ch4'  in Vars : axM[i,j].set_ylabel(r'$Y_{CH_4}~[-]$',fontsize=20)

#=====> Axial profiles
if Case=='Jaravel' :
	Data_v=ReadTNF(Files_v[-1],Params_v) ; (Vel,Tke)=Var_dyn(Data_v,Case)
	Data_s=ReadTNF(Files_s[-1],Params_s)
	M_d=zeros((Np,len(Vars)))
	for n,v in enumerate(Vars) :
		M_d[:,n]=D_int[v][-1]
		figA,axA=plt.subplots(figsize=(10,7))
		# figA.suptitle(Titre[v],fontsize=30)
		axA.set_ylabel(Titre[v]             ,fontsize=30)
		axA.set_xlabel('Axial position [mm]',fontsize=30)
		if   v=='Vel' : Plot_Axial( axA , Vx*1e3,D_int[v][-1] , Data_v[Abs[0]]*Abs[1]  ,Vel                , Err[v] )
		elif v=='k'   : Plot_Axial( axA , Vx*1e3,D_int[v][-1] , Data_v[Abs[0]]*Abs[1]  ,Tke                , Err[v] )
		elif v=='co'  : Plot_Axial( axA , Vx*1e3,D_int[v][-1] , Data_s[Cor['r']]*D0*1e3,Data_s[Cor[v]]*1e2 , Err[v] )
		else          : Plot_Axial( axA , Vx*1e3,D_int[v][-1] , Data_s[Cor['r']]*D0*1e3,Data_s[Cor[v]]     , Err[v] )
		if SAVE : util.SaveFig(figA,dirp+'Axial-%s.'%(v)+type)
	if DATA :
		St_data=dird+'Axial-Profiles.csv'
		with open(St_data,'w') as file :
			writer=csv.writer(file)
			writer.writerow( ['x(mm)'] + Vars )
			writer.writerows( column_stack( ( Vx*1e3 , M_d ) ) )
		file.closed

#=====================================================================================
if SAVE : Saving(dirp,Vars,type)
else    : util.Section('No Saving',0,1,'y')
#=====================================================================================