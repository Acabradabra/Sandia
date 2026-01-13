#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#%%=================================================================================
from unittest import case
from IPython import get_ipython

ip = get_ipython()
if ip is not None:
    ip.run_line_magic("load_ext", "autoreload")
    ip.run_line_magic("autoreload", "2")
#%%=================================================================================
#                     Modules
#===================================================================================
import Utilities as util
(Sysa,NSysa,Arg)=util.Parseur(['Save','Visu','Temp','Data','Over','Vel','Only','NoProf','Compa','TProf','Zoom','Struct','Report','All'],1,'Arg : Case (Jaravel,Garnier,Sevault)')
(                             [ SAVE , VISU , TEMP , DATA , OVER , VEL , ONLY , NOPROF , COMPA , TPROF , ZOOM , STRUCT , REPORT , ALL ])=Arg
if ALL : SAVE,VISU,TEMP,OVER,DATA=True,True,True,True,True # Over : Overwrite
if VEL : VISU,SAVE=True,True
if COMPA : SAVE,DATA=True,True

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
cmesh=1e3

#====================> Burner
Case=Sysa[0]
if   Case=='Jaravel' : from ParamsJaravel import * ; title=''
elif Case=='Garnier' : from ParamsGarnier import * ; title='he {:.0f} %'.format(he*100)
elif Case=='Sevault' : from ParamsSevault import * ; title='re {:.0f} k  ,  h2 {:.0f} %'.format(re,h2)
elif Case=='PRECIZE' : from ParamsPRECIZE import * ; title='PRECIZE'
else : sys.exit('=> Error : Case not recognized')

#====================> Fields
# Vars=['Vel','k','T','mix']
Vars=['Vel','k','T','o2','h2','ch4','co2','co']
# Vars=['Vel','k','T','mixH','o2','h2','n2','h2o']
# Vars=['T','o2','h2','ch4','co2','co']
# Vars=['Vel','T','o2','h2']
# Vars=['Vel','k','tt']
# Vars=['Vel','k']
# Vars=['Vel']
# Vars=['T']

#====================> Visu domain
if Case=='PRECIZE' :
	RY_t=[]
	RX_t=[]
	RY_d=[]
	RX_d=[]	
elif ZOOM :
	RY_t=[ 0,0.010]
	RX_t=[ 0,0.025]
	RY_d=[ 0,0.010]
	RX_d=[ 0,0.025]
elif VISU :
	RY_t=[ 0,0.10]
	RX_t=[Lc,0.50]
	RY_d=[ 0,0.04]
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
if   Case=='Jaravel' : sys.exit('=> Error : Jaravel case not implemented yet')
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
	for rf in os.popen('ls %s/report-*-rfile.out'%(dirc)).read().splitlines() :
		r_name=rf.split('/')[-1][7:-10]
		Dr=fl.Report_read(rf) ; Keys=list(Dr.keys())
		print('=> Reading report file : %s  ->  Keys : '%(r_name) , Keys )
		figr,axr=plt.subplots(figsize=(10,7))
		for k in Keys[1:] : axr.plot( Dr['Iteration'],Dr[k],label=k[7:] )
		if len(Keys)>2 : axr.legend(fontsize=15) #,loc='center',bbox_to_anchor=(0.5,1.5))
		if r_name in ['heat-release','mass-balance','shell-loss'] :
			axr.ticklabel_format( axis='y' , scilimits=(-2,2) )
		util.SaveFig(figr,dirp+'Report-%s.pdf'%(r_name))
	if ONLY : sys.exit('=> Stop after report reading')
#===================================================================================
if TEMP :
	# dirp=dirc+'PLOT/'
	util.Section( 'Temporals : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	fl.Probe_plot(dirc+'probe-F.out',dirp+'Probe-F.pdf')
	fl.Probe_plot(dirc+'probe-T.out',dirp+'Probe-T.pdf')
	fl.Probe_plot(dirc+'probe-V.out',dirp+'Probe-V.pdf')
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
	util.Section( 'Visualisation : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	if Case=='Garnier' :
		Vp_v=[ Lc+n*D0*Ld for n in Pos_v ]
		Vp_s=[ Lc+n*D0*Ld for n in Pos_s ]
	elif Case=='Sevault' :
		Vp_v=[ Lc+z*1e-3       for z in Pos_s ]
		Vp_s=[ Lc+z*1e-3       for z in Pos_s ]
	elif Case=='PRECIZE' :
		Vp_v=[]
		Vp_s=[]
	if STRUCT : ParT=['MIX',[fl.Mol_m,BC_f,BC_o]]
	else      : ParT=[]
	if   Case in ['Garnier','Sevault'] : 
		Tiso=list(array([1600,1800])-0)
		ParT+=['RECIRC','b']
		fsize=(25,5)
	elif Case=='PRECIZE' : 
		Tiso=[2000,2500]
		fsize=(15,7)

	F_int={}
	if 'T'    in Vars : F_int['T'   ]=fl.Visu(dird+slice,'temperature'        ,'Temperature [K]'     ,RX_t,RY_t,arange(250,3000,250),cmesh,[],fsize,'inferno',dirp+'Visu-Temperature.png',['INTERP','LINES',Vp_s,'ISO',Tiso]+ParT)
	if 'Vel'  in Vars : F_int['Vel' ]=fl.Visu(dird+slice,'velocity-magnitude' ,'Velocity [m/s]'      ,RX_d,RY_d,arange(0,350,50)    ,cmesh,[],fsize,'cividis',dirp+'Visu-Velocity.png'   ,['INTERP','LINES',Vp_v])
	if 'k'    in Vars : F_int['k'   ]=fl.Visu(dird+slice,'turb-kinetic-energy','k [$m^2/s^2$]'       ,RX_d,RY_d,arange(250,2500,250),cmesh,[],fsize,'cividis',dirp+'Visu-TKE.png'        ,['INTERP'])
	if 'tt'   in Vars : F_int['tt'  ]=fl.Visu(dird+slice,'tt'                 ,'tt [s]'              ,RX_d,RY_d,arange(0,1e-3,1e-4) ,cmesh,[],fsize,'cividis',dirp+'Visu-tt.png'         ,['INTERP'])
	if 'mixH' in Vars : F_int['mixH']=fl.Visu(dird+slice,'mixH'               ,'Mixture fraction [-]',RX_t,RY_t,arange(0,1.1,0.1)   ,cmesh,[],fsize,'viridis',dirp+'Visu-Mix.png'        ,['INTERP','MIXH',[fl.Mol_m,BC_f,BC_o]])
	if 'mix'  in Vars : F_int['mix' ]=fl.Visu(dird+slice,'fmean'              ,'Mixture fraction [-]',RX_t,RY_t,arange(0,1.1,0.1)   ,cmesh,[],fsize,'viridis',dirp+'Visu-Fmean.png'      ,['INTERP'])
	if 'h2'   in Vars : F_int['h2'  ]=fl.Visu(dird+slice,'h2'                 ,'Y H2 [-]'            ,RX_t,RY_t,arange(0,1.1,0.1)   ,cmesh,[],fsize,'viridis',dirp+'Visu-H2.png'         ,['INTERP'])
	if 'n2'   in Vars : F_int['n2'  ]=fl.Visu(dird+slice,'n2'                 ,'Y N2 [-]'            ,RX_t,RY_t,arange(0,1.1,0.1)   ,cmesh,[],fsize,'viridis',dirp+'Visu-N2.png'         ,['INTERP'])
	if 'o2'   in Vars : F_int['o2'  ]=fl.Visu(dird+slice,'o2'                 ,'Y O2 [-]'            ,RX_t,RY_t,arange(0,1.1,0.1)   ,cmesh,[],fsize,'viridis',dirp+'Visu-O2.png'         ,['INTERP'])
	if 'co'   in Vars : F_int['co'  ]=fl.Visu(dird+slice,'co'                 ,'Y CO [-]'            ,RX_t,RY_t,arange(0,1.1,0.1)   ,cmesh,[],fsize,'viridis',dirp+'Visu-CO.png'         ,['INTERP'])
	if 'ch4'  in Vars : F_int['ch4' ]=fl.Visu(dird+slice,'ch4'                ,'Y CH4 [-]'           ,RX_t,RY_t,arange(0,1.1,0.1)   ,cmesh,[],fsize,'viridis',dirp+'Visu-CH4.png'        ,['INTERP'])
	if 'h2o'  in Vars : F_int['h2o' ]=fl.Visu(dird+slice,'h2o'                ,'Y H2O [-]'           ,RX_t,RY_t,arange(0,1.1,0.1)   ,cmesh,[],fsize,'viridis',dirp+'Visu-H2O.png'        ,['INTERP'])
	if 'co2'  in Vars : F_int['co2' ]=fl.Visu(dird+slice,'co2'                ,'Y CO2 [-]'           ,RX_t,RY_t,arange(0,1.1,0.1)   ,cmesh,[],fsize,'viridis',dirp+'Visu-CO2.png'        ,['INTERP'])

	if not NOPROF :
		D_int={}
		Vy=linspace(0,0.5*D2,Np)
		for v in Vars :
			if v in ['Vel','k'] : D_int[v]=[ F_int[v](Np*[p],Vy) for p in Vp_v ]
			else                : D_int[v]=[ F_int[v](Np*[p],Vy) for p in Vp_s ]

	if DATA :
		util.Section( 'Data extraction : {:.3f} s'.format(time.time()-t0),1,5,'r' )
		for v in Vars :
			St_data=dirp+'Profiles-%s.csv'%(v)
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
if NOPROF : sys.exit('=> Stop after visu')
#%%=================================================================================
#                     Profiles
#===================================================================================
def ReadTNF(name,Params) :
	sep,skip,i0=Params[:3]
	with open(name) as file :
		Lines=file.readlines()
		if i0>=0 : T=[ s.strip() for s in Lines[skip][i0:-1].split(sep) if s ]
		else     : T=Params[3]
		M=array([ [float(v) for v in L.split(sep) if v ] for L in Lines[skip+1:] ]) ; (np,nv)=M.shape
	file.closed
	return({ s:M[:,i] for i,s in enumerate(T[:nv]) })
#=====================================================================================
def tpos_garnier(pos) : return('x/L_vis=%.2f'%(pos))
def tpos_sevault(pos) : return('x(mm)=%.0f'%(pos))
#=====================================================================================
def PlotFile(ax,dir0,var,nc) :
	D=util.ReadCSV(dir0+'Profiles-%s.csv'%(var)) ; Keys=list(D.keys()) #; print('D keys : ',Keys)
	Rad=D[Keys[0]]
	for n,k in enumerate(Keys[1:]) :
		i= n//nc
		j= n-nc*i
		ax[i,j].plot(Rad,D[k]) #,label=dir0.split('/')[-2])
	ax[-1,-1].plot(1,1,'.',label=dir0.split('/')[-2][5:])
	ax[-1,-1].legend(loc='center',bbox_to_anchor=(0.5,0.5))
#=====================================================================================
def Profile(ax,Vy,D_int,Data,Cor,Err,p,var,tpos,xlim) :
	if len(Vy)>0 : grid=Vy*1e3
	if var in D_int.keys() : ax.plot(grid,D_int[var][p],'r')
	ax.errorbar(Data[Cor['r']],Data[Cor[var]],yerr=Err[var]*abs(Data[Cor[var]]),ecolor='k',color='k')
	ax.set_title(tpos(Pos_s[p]),fontsize=20)
	ax.ticklabel_format(style='sci',axis='y',scilimits=(-2,4))
	if xlim[1]>0 : ax.set_xlim(xlim)
#%%===================================================================================
nc=2 ; nr=Npos_s//nc+int(Npos_s%nc>0)

if 'Vel'  in Vars : figV,axV=plt.subplots(ncols=3 ,nrows=3 ,figsize=(13,10)) ; figV.suptitle('Velocity [m/s]'      ,fontsize=30)
if 'k'    in Vars : figK,axK=plt.subplots(ncols=3 ,nrows=3 ,figsize=(13,10)) ; figK.suptitle('TKE [$m^2/s^2$]'     ,fontsize=30)
if 'T'    in Vars : figT,axT=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figT.suptitle('Temperature [K]'     ,fontsize=30)
if 'mixH' in Vars : figZ,axZ=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figZ.suptitle('Mixture fraction [-]',fontsize=30)
if 'o2'   in Vars : figO,axO=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figO.suptitle('Y O2 [-]'            ,fontsize=30)
if 'h2'   in Vars : figH,axH=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figH.suptitle('Y H2 [-]'            ,fontsize=30)
if 'n2'   in Vars : figN,axN=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figN.suptitle('Y N2 [-]'            ,fontsize=30)
if 'co'   in Vars : figA,axA=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figA.suptitle('Y CO [-]'            ,fontsize=30)
if 'h2o'  in Vars : figP,axP=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figP.suptitle('Y H2O [-]'           ,fontsize=30)
if 'co2'  in Vars : figC,axC=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figC.suptitle('Y CO2 [-]'           ,fontsize=30)
if 'ch4'  in Vars : figM,axM=plt.subplots(ncols=nc,nrows=nr,figsize=(13,10)) ; figM.suptitle('Y CH4 [-]'           ,fontsize=30)
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

if Case=='Garnier' :
	#=====> Velocity profiles
	for n,p in enumerate(IPos_v) :
		i= n//3
		j= n-3*i
		# Data=ReadLDV(Files_v[p])
		Data=ReadTNF(Files_v[p],('\t',11,2))
		Vel=hypot( Data['u'] , Data['v'] )
		Tke=hypot( Data['varu'] , Data['varv'] )
		if 'Vel' in Vars : axV[i,j].plot( Vy*1e3,D_int['Vel' ][p],'r' ) ; axV[i,j].plot( Data['y'],Vel,'ok' ) ; axV[i,j].set_title('x/L_vis='+Pv[n],fontsize=20) ; axV[i,j].ticklabel_format(style='sci',axis='y',scilimits=(-2,4)) ; axV[i,j].set_xlim((0,Rlims[n]))
		if 'k'   in Vars : axK[i,j].plot( Vy*1e3,D_int['k'   ][p],'r' ) ; axK[i,j].plot( Data['y'],Tke,'ok' ) ; axK[i,j].set_title('x/L_vis='+Pv[n],fontsize=20) ; axK[i,j].ticklabel_format(style='sci',axis='y',scilimits=(-2,4)) ; axK[i,j].set_xlim((0,Rlims[n]))
		if i==2 : 
			if 'Vel' in Vars : axV[i,j].set_xlabel('r [mm]',fontsize=20)
			if 'k'   in Vars : axK[i,j].set_xlabel('r [mm]',fontsize=20)
	if 'Vel' in Vars : axV[0,0].plot( [0,Rlims[0]],2*[Umoy],':k' )

if   Case=='Garnier' : Params=( ' ',3, 0 )         ; tpos=tpos_garnier ; Rlims=Npos_s*[0]
elif Case=='Sevault' : Params=('\t',0,-1,Vars_TNF) ; tpos=tpos_sevault ; Rlims=[5,8,11,15,20]
else : sys.exit('=> Error : Case not built yet')
#=====> Scattering
for n,p in enumerate(IPos_s) :
	i= n//nc
	j= n-nc*i
	Data=ReadTNF(Files_s[p],Params)
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
	if j==0 and i==1 : 
		if 'T'    in Vars : axT[i,j].set_ylabel(r'$T~[K]$'       ,fontsize=20)
		if 'mixH' in Vars : axZ[i,j].set_ylabel(r'$Z~[-]$'       ,fontsize=20)
		if 'o2'   in Vars : axO[i,j].set_ylabel(r'$Y_{O_2}~[-]$' ,fontsize=20)
		if 'h2'   in Vars : axH[i,j].set_ylabel(r'$Y_{H_2}~[-]$' ,fontsize=20)
		if 'n2'   in Vars : axN[i,j].set_ylabel(r'$Y_{N_2}~[-]$' ,fontsize=20)
		if 'co'   in Vars : axA[i,j].set_ylabel(r'$Y_{CO}~[-]$'  ,fontsize=20)
		if 'h2o'  in Vars : axP[i,j].set_ylabel(r'$Y_{H_2O}~[-]$',fontsize=20)
		if 'co2'  in Vars : axC[i,j].set_ylabel(r'$Y_{CO_2}~[-]$',fontsize=20)
		if 'ch4'  in Vars : axM[i,j].set_ylabel(r'$Y_{CH_4}~[-]$',fontsize=20)

if COMPA :
	dirp=dirc+'PLOT/'
	D_compa=[ d for d in os.listdir(dirc) if d[:5]=='PLOT-' ]
	D_compa.sort()
	for d in D_compa :
		dc=dirc+d+'/'
		if 'T'    in Vars : PlotFile(axT,dc,'T'   ,nc)
		if 'mixH' in Vars : PlotFile(axZ,dc,'mixH',nc)
		if 'o2'   in Vars : PlotFile(axO,dc,'o2'  ,nc)
		if 'h2'   in Vars : PlotFile(axH,dc,'h2'  ,nc)
		if 'n2'   in Vars : PlotFile(axN,dc,'n2'  ,nc)
		if 'co'   in Vars : PlotFile(axA,dc,'co'  ,nc)
		if 'h2o'  in Vars : PlotFile(axP,dc,'h2o' ,nc)
		if 'co2'  in Vars : PlotFile(axC,dc,'co2' ,nc)
		if 'ch4'  in Vars : PlotFile(axM,dc,'ch4' ,nc)

if SAVE :
	if 'Vel'  in Vars : util.SaveFig(figV,dirp+'Profiles-V.pdf')
	if 'k'    in Vars : util.SaveFig(figK,dirp+'Profiles-K.pdf')
	if 'T'    in Vars : util.SaveFig(figT,dirp+'Profiles-T.pdf')
	if 'mixH' in Vars : util.SaveFig(figZ,dirp+'Profiles-Z.pdf')
	if 'o2'   in Vars : util.SaveFig(figO,dirp+'Profiles-O.pdf')
	if 'h2'   in Vars : util.SaveFig(figH,dirp+'Profiles-H.pdf')
	if 'n2'   in Vars : util.SaveFig(figN,dirp+'Profiles-N.pdf')
	if 'co'   in Vars : util.SaveFig(figA,dirp+'Profiles-A.pdf')
	if 'h2o'  in Vars : util.SaveFig(figP,dirp+'Profiles-P.pdf')
	if 'co2'  in Vars : util.SaveFig(figC,dirp+'Profiles-C.pdf')
	if 'ch4'  in Vars : util.SaveFig(figM,dirp+'Profiles-M.pdf')
else :
	util.Section('No Saving',0,1,'y')
	# plt.show()
# %%