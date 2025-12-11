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
(Sysa,NSysa,Arg)=util.Parseur(['Plot','Out'],1,'Arg : Case [ Jaravel | Garnier | Sevault ] ')
(                             [ PLOT , OUT ])=Arg

from numpy import *
import sys
import time
import Fluent as fl
import FluentWally as fw

t0=time.time()
(plt,mtp)=util.Plot0()
#%%=================================================================================
#                     Parameters
#===================================================================================
Case='PRECIZE' # Sysa[0]
if   Case=='Jaravel' : import ParamsJaravel as pm ; Species=fl.Spe_Laera
elif Case=='Garnier' : import ParamsGarnier as pm ; Species=fl.Spe_UCSD
elif Case=='Sevault' : import ParamsSevault as pm ; Species=fl.Spe_Laera_l1
elif Case=='PRECIZE' : import ParamsPRECIZE as pm ; Species=fl.Spe_Walter
else : sys.exit('=> Error : Case not recognized')

# file_in=pm.dird+'Interp-Big-Laera.ip'
file_in=pm.dird+'Interp-PRECIZE-TXT.ip'
file_ou=file_in[:-3]+'-New.ip'

# Adds=['mix']
# Adds=['mix','mix2']
# Adds=['H2Rad','Ignit']
# Adds=['Ignit']
# Adds=['Laera','Ignit']
Adds=['Walter']

Fuel='H2'
dil=10
# Tad=2400
Tad=3000
Tu=300
T0=300

#%%=================================================================================
util.Section('Reading : {:.3f}'.format(time.time()-t0),1,5,'r')
#===================================================================================
data_in=fw.FluentInterpolationParser(file_in)

#%%=================================================================================
util.Section('Processing : {:.3f}'.format(time.time()-t0),1,5,'r')
#===================================================================================
BC_f=pm.BC_f
BC_o=pm.BC_o
if BC_f['unit']=='X' : BC_f=fl.ConvBC_XY(BC_f)
if BC_o['unit']=='X' : BC_o=fl.ConvBC_XY(BC_o)

#===============> Integer
Nd=data_in.n_dimensions
Np=data_in.n_points
Nv=data_in.n_variables

#===============> Space
X0=data_in.get_data( 'x0' )
X1=data_in.get_data( 'x1' )

#===============> Species
if 'mix' in Adds or 'mix2' in Adds :
    if 'N2'  in Species : Y_n2 =data_in.get_data( 'species-'+str(Species.index('N2' )) ) ; print('=> n2  : {:.3f}  ,  {:.3f}'.format(min(Y_n2 ),max(Y_n2 )))
    if 'H2'  in Species : Y_h2 =data_in.get_data( 'species-'+str(Species.index('H2' )) ) ; print('=> h2  : {:.3f}  ,  {:.3f}'.format(min(Y_h2 ),max(Y_h2 )))
    if 'O2'  in Species : Y_o2 =data_in.get_data( 'species-'+str(Species.index('O2' )) ) ; print('=> o2  : {:.3f}  ,  {:.3f}'.format(min(Y_o2 ),max(Y_o2 )))
    if 'H2O' in Species : Y_h2o=data_in.get_data( 'species-'+str(Species.index('H2O')) ) ; print('=> h2o : {:.3f}  ,  {:.3f}'.format(min(Y_h2o),max(Y_h2o)))
    if 'CO'  in Species : Y_co =data_in.get_data( 'species-'+str(Species.index('CO' )) ) ; print('=> co  : {:.3f}  ,  {:.3f}'.format(min(Y_co ),max(Y_co )))
    if 'CO2' in Species : Y_co2=data_in.get_data( 'species-'+str(Species.index('CO2')) ) ; print('=> co2 : {:.3f}  ,  {:.3f}'.format(min(Y_co2),max(Y_co2)))
    if 'CH4' in Species : Y_ch4=data_in.get_data( 'species-'+str(Species.index('CH4')) ) ; print('=> ch4 : {:.3f}  ,  {:.3f}'.format(min(Y_ch4),max(Y_ch4)))
    Temp =data_in.get_data( 'temperature' )
    Vx =data_in.get_data( 'x-velocity' )
    Vy =data_in.get_data( 'y-velocity' )
    Vel=hypot(Vx,Vy)
    Yc_f=fl.Yc(pm.BC_f,fl.Mol_m) ; Yo_f=fl.Yo(pm.BC_f,fl.Mol_m) ; Yh_f=fl.Yh(pm.BC_f,fl.Mol_m) #; Yn_m=pm.BC_f['N2']
    Yc_o=fl.Yc(pm.BC_o,fl.Mol_m) ; Yo_o=fl.Yo(pm.BC_o,fl.Mol_m) ; Yh_o=fl.Yh(pm.BC_o,fl.Mol_m) #; Yn_o=pm.BC_o['N2']
    if 'CH4' in Species : Yc=fl.Yc( {'CH4':Y_ch4,'CO2':Y_co2,'CO':Y_co            } ,fl.Mol_m )
    if 'CO2' in Species : Yo=fl.Yo( {'O2' :Y_o2 ,'CO2':Y_co2,'CO':Y_co,'H2O':Y_h2o} ,fl.Mol_m )
    if 'H2'  in Species : Yh=fl.Yh( {'CH4':Y_ch4,'H2O':Y_h2o,'H2':Y_h2            } ,fl.Mol_m )
    if 'N2'  in Species : Yn=Y_n2
    Mix_f=(Yh-Yh_o)/(Yh_f-Yh_o) ; print('=> Mix_fh : {:.3f}  ,  {:.3f}'.format(min(Mix_f),max(Mix_f)))
    Mix_o=(Yo-Yo_f)/(Yo_o-Yo_f) ; print('=> Mix_o : {:.3f}  ,  {:.3f}'.format(min(Mix_o),max(Mix_o)))
    if 'mix2' in Adds : 
        # Ys=1-(Yc+Yo)
        # Yt=Y_o2+Y_h2o+Y_ch4+Y_co+Y_co2+Y_n2 ; print('=> Yt  : {:.3f}  ,  {:.3f}'.format(min(Yt),max(Yt)))
        # Mix_f=(Yc-Yc_o)/(Yc_f-Yc_o) ; print('=> Mix_fc : {:.3f}  ,  {:.3f}'.format(min(Mix_f),max(Mix_f)))
        # Mix_n=(Yn-Yn_f)/(Yn_o-Yn_f) ; print('=> Mix_n : {:.3f}  ,  {:.3f}'.format(min(Mix_n),max(Mix_n)))
        # Mix_s=1-(Mix_f+Mix_n)       ; print('=> Mix_s : {:.3f}  ,  {:.3f}'.format(min(Mix_s),max(Mix_s)))
        # Zf=Y_ch4+Y_o2+Y_n2
        # Zp=Y_o2+Y_h2o+Y_co2+Y_n2
        # Zo=Y_o2+Y_n2
        # Zs=1-(Zf+Zo)
        Zf=Y_ch4
        Zf_m=pm.Y_m['CH4']
        Zf_p=pm.Y_p['CH4']
        Zf_o=pm.Y_o['CH4']
        Zp= Y_h2o+Y_co2+Y_co #+ Y_o2
        Zp_m=pm.Y_m['H2O']+pm.Y_m['CO2']+pm.Y_m['CO'] #+pm.Y_m['O2']
        Zp_p=pm.Y_p['H2O']+pm.Y_p['CO2']+pm.Y_p['CO'] #+pm.Y_p['O2']
        Zp_o=pm.Y_o['H2O']+pm.Y_o['CO2']+pm.Y_o['CO'] #+pm.Y_o['O2']
        Zo=1-(Zf+Zp)
        F_f=(Zf-Zf_o)/(Zf_m-Zf_o) ; print('=> F_f : {:.3f}  ,  {:.3f}'.format(min(F_f),max(F_f)))
        F_p=(Zp-Zp_o)/(Zp_p-Zp_o) ; print('=> F_p : {:.3f}  ,  {:.3f}'.format(min(F_p),max(F_p)))
        F_o=1-(F_f+F_p)           ; print('=> F_o : {:.3f}  ,  {:.3f}'.format(min(F_o),max(F_o)))
if 'H2Rad' in Adds :
    Spe_in=fl.Spe_H2Air
    Spe_ou=fl.Spe_UCSD
    Nspe=[]
    Ns_in,Ns_ou=len(Spe_in),len(Spe_ou)
    DSPE_ou=zeros((Np,Ns_ou))
    for n,s in enumerate(Spe_ou) :
        Nspe.append( 'species-'+str(n) )
        if s in Spe_in :
            i=Spe_in.index(s)
            DSPE_ou[:,n]=data_in.get_data( 'species-'+str(i) )
        else :
            print('=> Species {} not in input file'.format(s))
elif 'Laera' in Adds or 'Walter' in Adds :
    if 'Laera' in Adds :
        Spe_in=fl.Spe_Laera_l0
        Spe_ou=fl.Spe_Laera_l1
    elif 'Walter' in Adds :
        # Spe_in=fl.Spe_Laera_l1
        # Spe_ou=fl.Spe_Walter
        Spe_in=fl.Spe_Walter
        Spe_ou=fl.Spe_Laera_l1
    Nspe=[]
    Ns_in,Ns_ou=len(Spe_in),len(Spe_ou)
    DSPE_ou=zeros((Np,Ns_ou))
    for n,s in enumerate(Spe_ou) :
        Nspe.append( 'species-'+str(n) )
        if s in Spe_in :
            i=Spe_in.index(s)
            DSPE_ou[:,n]=data_in.get_data( 'species-'+str(i) )
        else :
            print('=> Species {} not in input file'.format(s))
else :
    Nspe=[ data_in.variable_names[n] for n in range(Nd,Nv+Nd) if 'species' in data_in.variable_names[n] ] ; print('=> Species : ',Nspe)
    Ns_in=len(Nspe) ; Ns_ou=Ns_in
    DSPE_ou=zeros((Np,Ns_in))
    for n,s in enumerate(Nspe) : 
        DSPE_ou[:,n]=data_in.get_data( s )
if 'Ignit' in Adds :
    # If=Spe_in.index(Fuel)
    # Yf=data_in.get_data( 'species-'+str(If) ) ; Yf0=min(Yf) ; Yf1=max(Yf)
    # C=1-(Yf-Yf0)/(Yf1-Yf0) ; C=clip(C,0,1)
    T_in=data_in.get_data( 'temperature' )
    C=dil*(T_in-T0)/(max(T_in)-T0) ; C=clip(C,0,1)
    Temp=Tu+C*(Tad-Tu)

#===============> Visualisation
if PLOT :
    c_inf=mtp.colormaps['inferno']
    c_civ=mtp.colormaps['cividis']
    c_vir=mtp.colormaps['viridis']
    tri=mtp.tri.Triangulation(X0,X1)
    r0=0.5*pm.D0
    r1=0.5*pm.D1
    r2=0.5*pm.D2
    re=0.5*pm.ep
    if Case=='Jaravel' :
        CMask=[
            [0 ,pm.Lc,   pm.Lc,       0, 0,pm.Lc,   pm.Lc,   pm.Ls,pm.Ls,0 ],
            [r0,   r0,r0+pm.ep,r0+pm.ep,r1,   r1,r1+pm.ep,r1+pm.ep,   r2,r2]]
    elif Case=='Garnier' : pass
    elif Case=='Sevault' :
        CMask=[
            [0 ,pm.Lc,pm.Lc, 0],
            [r0,   r0,   r1,r1]]
    if 'mix' in Adds or 'mix2' in Adds :
        fl.Field2(                      tri,Vel  ,'Vel [m/s]',False,[0,0.5],[0,0.75*r2],0,arange(0  ,500 ,10 ),0,c_civ,CMask,True,pm.dirp+'Interp-Vel.png'  ,(20,5))
        fl.Field2(                      tri,Temp ,'Temp [K]' ,False,[0,0.5],[0,0.75*r2],0,arange(0  ,3000,100),0,c_inf,CMask,True,pm.dirp+'Interp-Temp.png' ,(20,5))
        if 'CH4' in Species : fl.Field2(tri,Y_ch4,'Y_ch4 [-]',False,[0,0.5],[0,0.75*r2],0,arange(0  ,1.1,0.1 ),0,c_vir,CMask,True,pm.dirp+'Interp-CH4.png'  ,(20,5))
        if 'CO2' in Species : fl.Field2(tri,Y_co2,'Y_co2 [-]',False,[0,0.5],[0,0.75*r2],0,arange(0  ,1.1,0.1 ),0,c_vir,CMask,True,pm.dirp+'Interp-CO2.png'  ,(20,5))
        if 'H2O' in Species : fl.Field2(tri,Y_co2,'Y_h2o [-]',False,[0,0.5],[0,0.75*r2],0,arange(0  ,1.1,0.1 ),0,c_vir,CMask,True,pm.dirp+'Interp-H2O.png'  ,(20,5))
        if 'N2'  in Species : fl.Field2(tri,Y_n2 ,'Y_n2 [-]' ,False,[0,0.5],[0,0.75*r2],0,arange(0  ,1.1,0.1 ),0,c_vir,CMask,True,pm.dirp+'Interp-N2.png'   ,(20,5))
        if 'O2'  in Species : fl.Field2(tri,Y_o2 ,'Y_o2 [-]' ,False,[0,0.5],[0,0.75*r2],0,arange(0  ,1.1,0.1 ),0,c_vir,CMask,True,pm.dirp+'Interp-O2.png'   ,(20,5))
        if 'H2'  in Species : fl.Field2(tri,Y_h2 ,'Y_h2 [-]' ,False,[0,0.5],[0,0.75*r2],0,arange(0  ,1.1,0.1 ),0,c_vir,CMask,True,pm.dirp+'Interp-H2.png'   ,(20,5))
        if 'H2'  in Species : fl.Field2(tri,Yh   ,'Yh [-]'   ,False,[0,0.5],[0,0.75*r2],0,arange(0  ,1.1,0.1 ),0,c_vir,CMask,True,pm.dirp+'Interp-Yh.png'   ,(20,5))
        if 'CH4' in Species : fl.Field2(tri,Yc   ,'Yc [-]'   ,False,[0,0.5],[0,0.75*r2],0,arange(0  ,1.1,0.1 ),0,c_vir,CMask,True,pm.dirp+'Interp-Yc.png'   ,(20,5))
        if 'O2'  in Species : fl.Field2(tri,Yo   ,'Yo [-]'   ,False,[0,0.5],[0,0.75*r2],0,arange(0  ,1.1,0.1 ),0,c_vir,CMask,True,pm.dirp+'Interp-Yo.png'   ,(20,5))
        fl.Field2(                      tri,Mix_f,'Mix_f [-]',False,[0,0.5],[0,0.75*r2],0,arange(0  ,1.1,0.1 ),0,c_vir,CMask,True,pm.dirp+'Interp-Mix_f.png',(20,5))
        fl.Field2(                      tri,Mix_o,'Mix_o [-]',False,[0,0.5],[0,0.75*r2],0,arange(0  ,1.1,0.1 ),0,c_vir,CMask,True,pm.dirp+'Interp-Mix_o.png',(20,5))
        # fl.Field2(tri,Ys   ,'Ys [-]'   ,False,[0,0.5],[0,0.75*r2],0,arange(0  ,1.1,0.1 ),0,c_vir,CMask,True,pm.dirp+'Interp-Ys.png'   ,(20,5))
        # fl.Field2(tri,Yt   ,'Yt [-]'   ,False,[0,0.5],[0,0.75*r2],0,arange(0.9,1.1,0.01),0,c_vir,CMask,True,pm.dirp+'Interp-Yt.png'   ,(20,5))
        # fl.Field2(tri,Mix_n,'Mix_n [-]',False,[0,0.5],[0,0.75*r2],0,arange(0  ,1.1,0.1 ),0,c_vir,CMask,True,pm.dirp+'Interp-Mix_n.png',(20,5))
        # fl.Field2(tri,Mix_s,'Mix_s [-]',False,[0,0.5],[0,0.75*r2],0,arange(0  ,1.1,0.1 ),0,c_vir,CMask,True,pm.dirp+'Interp-Mix_s.png',(20,5))
    if 'mix2' in Adds :
        fl.Field2(tri,Zf ,'Zf [-]' ,False,[0,0.5],[],0,arange(0,1.1,0.1),0,c_vir,CMask,True,'Plot/Interp-Zf.png',(20,5))
        fl.Field2(tri,Zp ,'Zp [-]' ,False,[0,0.5],[],0,arange(0,1.1,0.1),0,c_vir,CMask,True,'Plot/Interp-Zp.png',(20,5))
        fl.Field2(tri,Zo ,'Zo [-]' ,False,[0,0.5],[],0,arange(0,1.1,0.1),0,c_vir,CMask,True,'Plot/Interp-Zo.png',(20,5))
        # fl.Field2(tri,Zs ,'Zs [-]' ,False,[0,0.5],[],0,arange(0,1.1,0.1),0,c_vir,CMask,True,'Plot/Interp-Zs.png',(20,5))
        fl.Field2(tri,F_f ,'F_f [-]' ,False,[0,0.5],[],0,arange(0,1.1,0.1),0,c_vir,CMask,True,'Plot/Interp-F_f.png',(20,5))
        fl.Field2(tri,F_p ,'F_p [-]' ,False,[0,0.5],[],0,arange(0,1.1,0.1),0,c_vir,CMask,True,'Plot/Interp-F_p.png',(20,5))
        fl.Field2(tri,F_o ,'F_o [-]' ,False,[0,0.5],[],0,arange(0,1.1,0.1),0,c_vir,CMask,True,'Plot/Interp-F_o.png',(20,5))
    if 'Ignit' in Adds :
        fl.Field2(tri,Temp ,'Temp [K]' ,False,[0,0.5],[0,r2],0,[],0,c_inf,CMask,True,'Plot/Visu-Temp.png' ,(20,5))

#%%=================================================================================
util.Section('Writing : {:.3f}'.format(time.time()-t0),1,5,'r')
#===================================================================================

#===============> Output
V_mesh=     [ v for v in data_in.variable_names if v in ['x0','x1','x2'] ]
V_spec=sort([ v for v in data_in.variable_names if 'species' in v ])
V_flow=     [ v for v in data_in.variable_names if not v in V_mesh and not v in V_spec ]

# Fields_ou=data_in.variable_names[Nd:-Ns_in]#+Nspe
Fields_ou=V_flow[:]
if    'mix'  in Adds : Fields_ou+=['fmean']
if    'mix2' in Adds : Fields_ou+=['fmean2']
if not 'mix' in Adds : Fields_ou+=Nspe
Nv_ou=len(Fields_ou) ; Nv0=Nd+Nv-Ns_in
print('=> Fields out : ',Fields_ou)

DATA_ou=zeros((Np,Nd+Nv_ou))
# for n,v in enumerate(data_in.variable_names[:-Ns_in]) : 
for n,v in enumerate(V_mesh+V_flow) : 
    print(n,v)
    DATA_ou[:,n]=data_in.get_data(v)
if   'mix'  in Adds : DATA_ou[:,Nv0 ]= Mix_f 
elif 'mix2' in Adds : DATA_ou[:,Nv0:]=[Yf,Yp]
else                : DATA_ou[:,Nv0:Nv0+Ns_ou+1]=DSPE_ou
if 'Ignit'  in Adds : DATA_ou[:,Nd+Fields_ou.index('temperature')]=Temp 

if not OUT : sys.exit('=> End of processing')
fou=open(file_ou,'w')
fou.write( '3\n' )
fou.write( '{}\n'.format(Nd) )
fou.write( '{}\n'.format(Np) )
fou.write( '{}\n'.format(Nv_ou) )
for f in Fields_ou :
    fou.write( '{}\n'.format(f) )
for v in range(Nv_ou+Nd) :
    fou.write('(')
    for p in range(Np) : fou.write( ' {:.12e}\n'.format(DATA_ou[p,v]) )
    fou.write(')\n')
fou.close()

#%%=================================================================================
util.Section('Programme completed : {:.3f}'.format(time.time()-t0),1,5,'r')
#===================================================================================