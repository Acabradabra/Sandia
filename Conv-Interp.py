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
(Sysa,NSysa,Arg)=util.Parseur(['Plot'],0,'Arg : ')
(                             [ PLOT ])=Arg

from numpy import *
import sys
import time
import Fluent as fl
import FluentWally as fw
import Params as pm

t0=time.time()
(plt,mtp)=util.Plot0()
#%%=================================================================================
#                     Parameters
#===================================================================================

# file_in='/mnt/scratch/PRECIZE/Sandia-Jaravel/RUN-D100-02-Laera/INIT-Fluent/Interp-00-Fluent.ip'
file_in='/mnt/scratch/PRECIZE/Sandia-Jaravel/RUN-D100-02-Laera/DUMP-02-EDC-PB-OD2-FD3/DATA/Interp-EDC-FD3.ip'
file_ou=file_in[:-3]+'-New.ip'

Adds=['mix','mix2']

#%%=================================================================================
util.Section('Reading : {:.3f}'.format(time.time()-t0),1,5,'r')
#===================================================================================
data_in=fw.FluentInterpolationParser(file_in)

#%%=================================================================================
util.Section('Processing : {:.3f}'.format(time.time()-t0),1,5,'r')
#===================================================================================

#===============> Integer
Nd=data_in.n_dimensions
Np=data_in.n_points
Nv=data_in.n_variables

#===============> Space
X0=data_in.get_data( 'x0' )
X1=data_in.get_data( 'x1' )

#===============> Species
if 'mix' or 'mix2' in Adds :
    Species=['O2','H2O','CH4','CO','CO2','H2','H','O','OH','HO2','H2O2','CH3','CH2O','CH3O','CH3OH','C2H2','C2H4','C2H6','CH2CO','CH','CH2','CH2(S)','HCO','CH2OH','C2H3','C2H5','HCCO','CH2CHO','N2']
    I_h2 =Species.index('H2' )
    I_n2 =Species.index('N2' )
    I_o2 =Species.index('O2' )
    I_h2o=Species.index('H2O')
    I_co =Species.index('CO' )
    I_co2=Species.index('CO2')
    I_ch4=Species.index('CH4')
    Y_h2 =data_in.get_data( 'species-'+str(I_h2 ) ) ; print('=> h2  : {:.3f}  ,  {:.3f}'.format(min(Y_h2 ),max(Y_h2 )))
    Y_n2 =data_in.get_data( 'species-'+str(I_n2 ) ) ; print('=> n2  : {:.3f}  ,  {:.3f}'.format(min(Y_n2 ),max(Y_n2 )))
    Y_o2 =data_in.get_data( 'species-'+str(I_o2 ) ) ; print('=> o2  : {:.3f}  ,  {:.3f}'.format(min(Y_o2 ),max(Y_o2 )))
    Y_h2o=data_in.get_data( 'species-'+str(I_h2o) ) ; print('=> h2o : {:.3f}  ,  {:.3f}'.format(min(Y_h2o),max(Y_h2o)))
    Y_co =data_in.get_data( 'species-'+str(I_co ) ) ; print('=> co  : {:.3f}  ,  {:.3f}'.format(min(Y_co ),max(Y_co )))
    Y_co2=data_in.get_data( 'species-'+str(I_co2) ) ; print('=> co2 : {:.3f}  ,  {:.3f}'.format(min(Y_co2),max(Y_co2)))
    Y_ch4=data_in.get_data( 'species-'+str(I_ch4) ) ; print('=> ch4 : {:.3f}  ,  {:.3f}'.format(min(Y_ch4),max(Y_ch4)))
    Temp =data_in.get_data( 'temperature' )
    Vx =data_in.get_data( 'x-velocity' )
    Vy =data_in.get_data( 'y-velocity' )
    Vel=hypot(Vx,Vy)
    Yc_m=fl.Yc(pm.Y_m,fl.Mol_m) ; Yo_m=fl.Yo(pm.Y_m,fl.Mol_m) ; Yn_m=pm.Y_m['N2']
    Yc_p=fl.Yc(pm.Y_p,fl.Mol_m) ; Yo_p=fl.Yo(pm.Y_p,fl.Mol_m) ; Yn_p=pm.Y_p['N2']
    Yc_o=fl.Yc(pm.Y_o,fl.Mol_m) ; Yo_o=fl.Yo(pm.Y_o,fl.Mol_m) ; Yn_o=pm.Y_o['N2']
    Dic_c={'CH4':Y_ch4,'CO2':Y_co2,'CO':Y_co}
    Dic_o={'O2' :Y_o2 ,'CO2':Y_co2,'CO':Y_co,'H2O':Y_h2o}
    Dic_h={'CH4':Y_ch4,'H2O':Y_h2o,'H2':Y_h2}
    Yc=fl.Yc( Dic_c , fl.Mol_m)
    Yo=fl.Yo( Dic_o , fl.Mol_m)
    Yh=fl.Yh( Dic_h , fl.Mol_m)
    Yn=Y_n2
    Ys=1-(Yc+Yo)
    Yt=Y_o2+Y_h2o+Y_ch4+Y_co+Y_co2+Y_n2 ; print('=> Yt  : {:.3f}  ,  {:.3f}'.format(min(Yt),max(Yt)))
    Mix_f=(Yc-Yc_o)/(Yc_m-Yc_o) ; print('=> Mix_f : {:.3f}  ,  {:.3f}'.format(min(Mix_f),max(Mix_f)))
    Mix_o=(Yo-Yo_m)/(Yo_o-Yo_m) ; print('=> Mix_o : {:.3f}  ,  {:.3f}'.format(min(Mix_o),max(Mix_o)))
    Mix_n=(Yn-Yn_m)/(Yn_o-Yn_m) ; print('=> Mix_n : {:.3f}  ,  {:.3f}'.format(min(Mix_n),max(Mix_n)))
    Mix_s=1-(Mix_f+Mix_n)       ; print('=> Mix_s : {:.3f}  ,  {:.3f}'.format(min(Mix_s),max(Mix_s)))
    if 'mix2' in Adds : 
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

#===============> Visualisation
PLOT=True
if PLOT :
    c_inf=mtp.colormaps['inferno']
    c_civ=mtp.colormaps['cividis']
    c_vir=mtp.colormaps['viridis']
    tri=mtp.tri.Triangulation(X0,X1)
    r0=0.5*pm.D0
    r1=0.5*pm.D1
    r2=0.5*pm.D2
    re=0.5*pm.ep
    CMask=[
        [0, pm.Lc,pm.Lc , 0,0 , pm.Lc,pm.Lc , pm.Ls,pm.Ls,0],
        [r0,r0 , r0+pm.ep,r0+pm.ep , r1,r1, r1+pm.ep,r1+pm.ep, r2,r2]
    ]
    if 'mix' or 'mix2' in Adds :
        fl.Field2(tri,Vel  ,'Vel [m/s]',False,[0,0.5],[0,r2],0,arange(0  ,100,10  ),c_civ,CMask,True,'Plot/Visu-Vel.png'  ,(20,5))
        fl.Field2(tri,Temp ,'Temp [K]' ,False,[0,0.5],[0,r2],0,arange(0  ,1.1,0.1 ),c_inf,CMask,True,'Plot/Visu-Temp.png' ,(20,5))
        fl.Field2(tri,Y_ch4,'Y_ch4 [-]',False,[0,0.5],[0,r2],0,arange(0  ,1.1,0.1 ),c_vir,CMask,True,'Plot/Visu-CH4.png'  ,(20,5))
        fl.Field2(tri,Y_o2 ,'Y_o2 [-]' ,False,[0,0.5],[0,r2],0,arange(0  ,1.1,0.1 ),c_vir,CMask,True,'Plot/Visu-O2.png'   ,(20,5))
        fl.Field2(tri,Y_n2 ,'Y_n2 [-]' ,False,[0,0.5],[0,r2],0,arange(0  ,1.1,0.1 ),c_vir,CMask,True,'Plot/Visu-N2.png'   ,(20,5))
        fl.Field2(tri,Y_h2 ,'Y_h2 [-]' ,False,[0,0.5],[0,r2],0,arange(0  ,1.1,0.1 ),c_vir,CMask,True,'Plot/Visu-H2.png'   ,(20,5))
        fl.Field2(tri,Yh   ,'Yh [-]'   ,False,[0,0.5],[0,r2],0,arange(0  ,1.1,0.1 ),c_vir,CMask,True,'Plot/Visu-Yh.png'   ,(20,5))
        fl.Field2(tri,Yc   ,'Yc [-]'   ,False,[0,0.5],[0,r2],0,arange(0  ,1.1,0.1 ),c_vir,CMask,True,'Plot/Visu-Yc.png'   ,(20,5))
        fl.Field2(tri,Yo   ,'Yo [-]'   ,False,[0,0.5],[0,r2],0,arange(0  ,1.1,0.1 ),c_vir,CMask,True,'Plot/Visu-Yo.png'   ,(20,5))
        fl.Field2(tri,Ys   ,'Ys [-]'   ,False,[0,0.5],[0,r2],0,arange(0  ,1.1,0.1 ),c_vir,CMask,True,'Plot/Visu-Ys.png'   ,(20,5))
        fl.Field2(tri,Yt   ,'Yt [-]'   ,False,[0,0.5],[0,r2],0,arange(0.9,1.1,0.01),c_vir,CMask,True,'Plot/Visu-Yt.png'   ,(20,5))
        fl.Field2(tri,Mix_f,'Mix_f [-]',False,[0,0.5],[0,r2],0,arange(0  ,1.1,0.1 ),c_vir,CMask,True,'Plot/Visu-Mix_f.png',(20,5))
        fl.Field2(tri,Mix_o,'Mix_o [-]',False,[0,0.5],[0,r2],0,arange(0  ,1.1,0.1 ),c_vir,CMask,True,'Plot/Visu-Mix_o.png',(20,5))
        fl.Field2(tri,Mix_n,'Mix_n [-]',False,[0,0.5],[0,r2],0,arange(0  ,1.1,0.1 ),c_vir,CMask,True,'Plot/Visu-Mix_n.png',(20,5))
        fl.Field2(tri,Mix_s,'Mix_s [-]',False,[0,0.5],[0,r2],0,arange(0  ,1.1,0.1 ),c_vir,CMask,True,'Plot/Visu-Mix_s.png',(20,5))
    if 'mix2' in Adds :
        fl.Field2(tri,Zf ,'Zf [-]' ,False,[0,0.5],[],0,arange(0,1.1,0.1),c_vir,CMask,True,'Plot/Visu-Zf.png',(20,5))
        fl.Field2(tri,Zp ,'Zp [-]' ,False,[0,0.5],[],0,arange(0,1.1,0.1),c_vir,CMask,True,'Plot/Visu-Zp.png',(20,5))
        fl.Field2(tri,Zo ,'Zo [-]' ,False,[0,0.5],[],0,arange(0,1.1,0.1),c_vir,CMask,True,'Plot/Visu-Zo.png',(20,5))
        # fl.Field2(tri,Zs ,'Zs [-]' ,False,[0,0.5],[],0,arange(0,1.1,0.1),c_vir,CMask,True,'Plot/Visu-Zs.png',(20,5))
        fl.Field2(tri,F_f ,'F_f [-]' ,False,[0,0.5],[],0,arange(0,1.1,0.1),c_vir,CMask,True,'Plot/Visu-F_f.png',(20,5))
        fl.Field2(tri,F_p ,'F_p [-]' ,False,[0,0.5],[],0,arange(0,1.1,0.1),c_vir,CMask,True,'Plot/Visu-F_p.png',(20,5))
        fl.Field2(tri,F_o ,'F_o [-]' ,False,[0,0.5],[],0,arange(0,1.1,0.1),c_vir,CMask,True,'Plot/Visu-F_o.png',(20,5))

#%%=================================================================================
util.Section('Writing : {:.3f}'.format(time.time()-t0),1,5,'r')
#===================================================================================

#===============> Output
DATA_ou=zeros((Np,Nv+Nd+len(Adds)))
for n,v in enumerate(data_in.variable_names) : DATA_ou[:,n]=data_in.get_data(v)
if   'mix'  in Adds : DATA_ou[:,Nv+Nd ]=Mix ; Fields_ou=data_in.variable_names[Nd:]+['fmean']
elif 'mix2' in Adds : DATA_ou[:,Nv+Nd:]=[Yf,Yp] ; Fields_ou=data_in.variable_names[Nd:]+['fmean','fmean2']
else                : Fields_ou=data_in.variable_names[Nd:]

Nv_ou=len(Fields_ou)

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