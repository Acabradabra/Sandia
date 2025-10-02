#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Utilities as util
(Sysa,NSysa,Arg)=util.Parseur(['BD','Conv','Mix','VelProf','Pos','Flow'],0,'Arg : ')
(                             [ BD , CONV , MIX , VELPROF , POS , FLOW ])=Arg

from numpy import *
import sys
import time
import pygmsh as pg
import gmsh   as gm
import Meshing as me
import meshio as mo

import matplotlib.pyplot as plt
import matplotlib

from Params import *

t0=time.time()
#===================================================================================
#                     Parameters
#===================================================================================

gdim = 2  # Geometric dimension of the mesh

#=====> files
# d_mesh='MESH/'
d_mesh='MESH-D100/'
d_plot='PLOT/'
name='Sandia-Jaravel'

#=====> Mesh Size
h0=1e-4
h1=1e-4
h2=2e-4
hf=2e-3

hfc=2e-4
hfs=5e-4

#=====> Flame refinement params
Lr=0.2
Lf=4*ep
Dr=1e-2

#=====> Boundary layer params
rb=1.1 # Ratio slices boundary layer
rt=1   # Ratio transition
rs=1   # Ratio surface caré/triangle
ra=0.7 # Aspect ratio last slice
rj=0.45
N0=10
N1=10
N2=20

#===================================================================================
#                     Processes
#===================================================================================

if FLOW :

	Umax=61
	n=round( 2/(Umax/Um-1),0 )
	Umax=Um*(n+2)/n

	print('\n=> n : {:.0f}  ,  Umax : {:.12f} \n'.format(n,Umax))

	sys.exit('=> Flow rates computed')

if POS :
	#===================================================================================
	util.Section( 'Profile Positions : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	#===================================================================================

	Vd=[1,2,3,7.5,15,30]
	for d in Vd : print('=> Pos : X/D {}   ,   X {:.5f}'.format(d,Lc+d*D0))

	sys.exit('=> Profile created')

if VELPROF :
	#===================================================================================
	util.Section( 'Velocity Profile : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	#===================================================================================
	def Pois(Um,r,r0,R,np) : return( Um*(1-(abs(r-r0)/R)**np) )
	#===================================================================================

	np=25
	nd=int(1e4)

	r0=0.5*D0
	r1=0.5*D1 ; Dr1=r1-(r0+ep)
	r2=0.5*D2 ; Dr2=r2-(r1+ep)

	#=====> Computation
	Vr_m=linspace( 0   ,r0,nd)
	Vr_p=linspace(r0+ep,r1,nd)
	Vr_c=linspace(r1+ep,r2,nd)
	Vu_m=Pois( Um,Vr_m, 0        , r0    ,np )
	Vu_p=Pois( Up,Vr_p,r1-0.5*Dr1,0.5*Dr1,np )
	Vu_c=Pois( Uc,Vr_c,r2        ,    Dr2,np )

	#=====> Ploting
	fig,ax=plt.subplots(ncols=3,figsize=(18,6))
	ax[0].plot(Vr_m*1e3,Vu_m,'k')
	ax[1].plot(Vr_p*1e3,Vu_p,'k')
	ax[2].plot(Vr_c*1e3,Vu_c,'k')
	fig.savefig('Plot/VelocityProfiles.pdf')

	#=====> Writing
	savetxt( 'DATA/Profile-Velocity-Main.dat'  ,transpose(array([Vr_m,Vu_m])),fmt='%.12e',delimiter=',',newline='\n',header='Radius,Velocity' )
	savetxt( 'DATA/Profile-Velocity-Pilot.dat' ,transpose(array([Vr_p,Vu_p])),fmt='%.12e',delimiter=',',newline='\n',header='Radius,Velocity' )
	savetxt( 'DATA/Profile-Velocity-Coflow.dat',transpose(array([Vr_c,Vu_c])),fmt='%.12e',delimiter=',',newline='\n',header='Radius,Velocity' )

	sys.exit('=> Profile created')

if MIX :
	#===================================================================================
	util.Section( 'Mixture fraction : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	#===================================================================================
	Yc_f=fl.Yc(Ych4[0],Yco2[0],Yco[0],Mol_m)
	Yc_p=fl.Yc(Ych4[1],Yco2[1],Yco[1],Mol_m)
	Yc_o=fl.Yc(Ych4[2],Yco2[2],Yco[2],Mol_m)
	Zp  =(Yc_p-Yc_o)/(Yc_f-Yc_o)
	
	print('=> Yc_f : {:.5f}'.format(Yc_f))
	print('=> Yc_p : {:.5f}'.format(Yc_p))
	print('=> Yc_o : {:.5f}'.format(Yc_o))
	print('=> Zp   : {:.5f}'.format(Zp  ))

	a=0.79/(0.21+0.79)
	Nt=3+2*a
	Xco2=1/Nt
	Xh2o=2/Nt
	Xn2 =2*a/Nt
	Mt=Xco2*Mol_m[2]+Xh2o*Mol_m[4]+Xn2*Mol_m[5]
	Yco2_0=Xco2*Mol_m[2]/Mt
	Cp=Yco2[1]/Yco2_0

	print('=> Cp : {:.5f}'.format(Cp))

	sys.exit()

if False :
	#===================================================================================
	util.Section( 'Conversion : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	#===================================================================================
	d_fluent='/mnt/d/FLUENT/Sandia-Jaravel/MESH/'
	# f_in=d_fluent+'burner-axisymmetric-2d.msh'
	# f_ou=d_fluent+'burner-axisymmetric-2d-ASCII.msh'

	mesh=mo.gmsh.read(d_mesh+name+'-MESH.msh')
	# mesh=mo.ansys.read(f_in)
	# mo.ansys.write(f_ou,binary=False)
	mo.ansys.write(d_mesh+name+'-ANSYS.msh',mesh,binary=False)
	# mo.cgns.write(d_mesh+name+'-MESH.cgns',mesh)
	sys.exit('=> File converted')

#===================================================================================
util.Section( 'Geometrie : {:.3f} s'.format(time.time()-t0),1,5,'r' )
#===================================================================================

#=====> X dim
rp=0.5*ep
r0=0.5*D0+rp
r1=0.5*D1+rp
r2=0.5*D2

#=====> Boundary layer
rat=0.25*sqrt(3)*rs ; print( '=> Aspect Triangle=Square : {:.3f}'.format(rat) )
(h00,Lt0)=me.h_smooth(rb,N0,rt*h0) ; print('=> Main   : h0 {:.2e} , Lt {:.2e} , ra*h/h0 {:.2e}'.format(h00,Lt0,ra*h0/h00))
(h10,Lt1)=me.h_smooth(rb,N1,rt*h1) ; print('=> Pilot  : h0 {:.2e} , Lt {:.2e} , ra*h/h0 {:.2e}'.format(h10,Lt1,ra*h1/h10))
(h20,Lt2)=me.h_smooth(rb,N2,rt*h2) ; print('=> Coflow : h0 {:.2e} , Lt {:.2e} , ra*h/h0 {:.2e}'.format(h20,Lt2,ra*h2/h20))
N0z=int(round(Lc/(h0*ra),0))+1
N1z=int(round(Lc/(h1*ra),0))+1
N2z=int(round(Ls/(h2*ra),0))+1
N0r=int(round( (0.5*D0        -  Lt0)/h0 ,0))+1
N1r=int(round( (0.5*D1-(r0+rp)-2*Lt1)/h1 ,0))+1
N2r=int(round( (r2    -(r1+rp)-  Lt2)/h2 ,0))+1
(Nt,rj2)=me.Nbl2(rj,ra,rb) ; print('=> Transition : N {}  ,  r2 {:.3f}'.format(Nt,rj2))
Ltj0=me.ht(h0*rj,rj2,Nt) ; print('=> Main   : L transition jet {:.3f} [mm]'.format(Ltj0*1e3))
Ltj1=me.ht(h1*rj,rj2,Nt) ; print('=> Pilot  : L transition jet {:.3f} [mm]'.format(Ltj1*1e3))
Ltj2=me.ht(h2*rj,rj2,Nt) ; print('=> Coflow : L transition jet {:.3f} [mm]'.format(Ltj2*1e3))

if BD : sys.exit('=> Boundary layer')

#=====> Y dim
y1=Lc-Ls
y2=Lc
y3=Lt+Lc

gm.initialize()
#=====> Points
Points=[
gm.model.geo.add_point(y2     ,0        ,0,meshSize=h0 ), #===> Main Jet
gm.model.geo.add_point(y2-Ltj0,0        ,0,meshSize=h0 ),
gm.model.geo.add_point( 0     ,0        ,0,meshSize=h0 ),
gm.model.geo.add_point( 0     ,r0-rp-Lt0,0,meshSize=h0 ),
gm.model.geo.add_point( 0     ,r0-rp    ,0,meshSize=h00),
gm.model.geo.add_point(y2-Ltj0,r0-rp    ,0,meshSize=h00),
gm.model.geo.add_point(y2     ,r0-rp    ,0,meshSize=h00),
gm.model.geo.add_point(y2     ,r0-rp-Lt0,0,meshSize=h0 ),
gm.model.geo.add_point(y2     ,r0+rp    ,0,meshSize=h10), #===> Pilot Jet
gm.model.geo.add_point(y2     ,r0+rp+Lt1,0,meshSize=h1 ), 
gm.model.geo.add_point(y2-Ltj1,r0+rp    ,0,meshSize=h10),
gm.model.geo.add_point( 0     ,r0+rp    ,0,meshSize=h10),
gm.model.geo.add_point( 0     ,r0+rp+Lt1,0,meshSize=h1 ),
gm.model.geo.add_point( 0     ,r1-rp-Lt1,0,meshSize=h1 ),
gm.model.geo.add_point( 0     ,r1-rp    ,0,meshSize=h10),
gm.model.geo.add_point(y2-Ltj1,r1-rp    ,0,meshSize=h10),
gm.model.geo.add_point(y2     ,r1-rp    ,0,meshSize=h10),
gm.model.geo.add_point(y2     ,r1-rp-Lt1,0,meshSize=h1 ),
gm.model.geo.add_point(y2     ,r1+rp    ,0,meshSize=h20), #===> Coflow
gm.model.geo.add_point(y2     ,r1+rp+Lt2,0,meshSize=h2 ),
gm.model.geo.add_point(y2-Ltj2,r1+rp    ,0,meshSize=h20),
gm.model.geo.add_point(y1     ,r1+rp    ,0,meshSize=h20),
gm.model.geo.add_point(y1     ,r1+rp+Lt2,0,meshSize=h2 ),
gm.model.geo.add_point(y1     ,r2       ,0,meshSize=h2 ),
gm.model.geo.add_point(y2-Ltj2,r2       ,0,meshSize=h2 ),
gm.model.geo.add_point(y2     ,r2       ,0,meshSize=h2 ),
gm.model.geo.add_point(y3     ,r2       ,0,meshSize=hf ), #===> Outlet
gm.model.geo.add_point(y3     ,0        ,0,meshSize=hf ) ]
Np=len(Points)
gm.model.geo.synchronize()

#=====> Lines
lma=gm.model.geo.add_line(Points[ 0],Points[ 1]) #===> Main left transition
lml=gm.model.geo.add_line(Points[ 1],Points[ 2]) #===> Main left
lmc=gm.model.geo.add_line(Points[ 2],Points[ 3]) #===> Main Inlet flow
lm1=gm.model.geo.add_line(Points[ 3],Points[ 4]) #===> Main Inlet BL
lmr=gm.model.geo.add_line(Points[ 4],Points[ 5]) #===> Main right
lmb=gm.model.geo.add_line(Points[ 5],Points[ 6]) #===> Main right transition
lm2=gm.model.geo.add_line(Points[ 6],Points[ 7]) #===> Main Outlet BL
lmo=gm.model.geo.add_line(Points[ 7],Points[ 0]) #===> Main Outlet flow
lme=gm.model.geo.add_line(Points[ 6],Points[ 8]) #===> Main leap

lp2=gm.model.geo.add_line(Points[ 8],Points[ 9]) #===> Pilot Outlet BL L
lpa=gm.model.geo.add_line(Points[ 8],Points[10]) #===> Pilot left transition
lpl=gm.model.geo.add_line(Points[10],Points[11]) #===> Pilot left
lp1=gm.model.geo.add_line(Points[11],Points[12]) #===> Pilot inlet BL L
lpc=gm.model.geo.add_line(Points[12],Points[13]) #===> Pilot inlet flow
lp3=gm.model.geo.add_line(Points[13],Points[14]) #===> Pilot inlet BL R
lpr=gm.model.geo.add_line(Points[14],Points[15]) #===> Pilot right
lpb=gm.model.geo.add_line(Points[15],Points[16]) #===> Pilot right transition
lp4=gm.model.geo.add_line(Points[16],Points[17]) #===> Pilot Outlet BL R
lpo=gm.model.geo.add_line(Points[17],Points[ 9]) #===> Pilot Outlet
lpe=gm.model.geo.add_line(Points[16],Points[18]) #===> Pilot leap

lc2=gm.model.geo.add_line(Points[18],Points[19]) #===> Coflow Outlet BL 
lca=gm.model.geo.add_line(Points[18],Points[20]) #===> Coflow Left transition
lcl=gm.model.geo.add_line(Points[20],Points[21]) #===> Coflow Left
lc1=gm.model.geo.add_line(Points[21],Points[22]) #===> Coflow Inlet BL
lcc=gm.model.geo.add_line(Points[22],Points[23]) #===> Coflow Inlet
lcr=gm.model.geo.add_line(Points[23],Points[24]) #===> Coflow right
lcb=gm.model.geo.add_line(Points[24],Points[25]) #===> Coflow right transition
lco=gm.model.geo.add_line(Points[25],Points[19]) #===> Coflow outlet

lor=gm.model.geo.add_line(Points[25],Points[26]) #===> Outlet right
loo=gm.model.geo.add_line(Points[26],Points[27]) #===> Outlet
lol=gm.model.geo.add_line(Points[27],Points[ 0]) #===> Outlet left
gm.model.geo.synchronize()

#=====> transfinite
gm.model.geo.mesh.setTransfiniteCurve(lm1,N0+1,meshType='Progression',coef=1/rb) #===> Main Jet
gm.model.geo.mesh.setTransfiniteCurve(lm2,N0+1,meshType='Progression',coef=  rb)
gm.model.geo.mesh.setTransfiniteCurve(lma,Nt+1,meshType='Progression',coef=  rj2)
gm.model.geo.mesh.setTransfiniteCurve(lmb,Nt+1,meshType='Progression',coef=1/rj2)
gm.model.geo.mesh.setTransfiniteCurve(lmc,N0r)
gm.model.geo.mesh.setTransfiniteCurve(lmo,N0r)
gm.model.geo.mesh.setTransfiniteCurve(lmr,N0z)
gm.model.geo.mesh.setTransfiniteCurve(lml,N0z)
gm.model.geo.mesh.setTransfiniteCurve(lp2,N1+1,meshType='Progression',coef=rb) #===> Pilot
gm.model.geo.mesh.setTransfiniteCurve(lp1,N1+1,meshType='Progression',coef=rb)
gm.model.geo.mesh.setTransfiniteCurve(lpa,Nt+1,meshType='Progression',coef=  rj2)
gm.model.geo.mesh.setTransfiniteCurve(lpb,Nt+1,meshType='Progression',coef=1/rj2)
gm.model.geo.mesh.setTransfiniteCurve(lpr,N1z)
gm.model.geo.mesh.setTransfiniteCurve(lpl,N1z)
gm.model.geo.mesh.setTransfiniteCurve(lp4,N1+1,meshType='Progression',coef=  rb)
gm.model.geo.mesh.setTransfiniteCurve(lp3,N1+1,meshType='Progression',coef=1/rb)
gm.model.geo.mesh.setTransfiniteCurve(lpc,N1r)
gm.model.geo.mesh.setTransfiniteCurve(lpo,N1r)
gm.model.geo.mesh.setTransfiniteCurve(lc1,N2+1,meshType='Progression',coef=rb) #===> Coflow
gm.model.geo.mesh.setTransfiniteCurve(lc2,N2+1,meshType='Progression',coef=rb)
gm.model.geo.mesh.setTransfiniteCurve(lca,Nt+1,meshType='Progression',coef=  rj2)
gm.model.geo.mesh.setTransfiniteCurve(lcb,Nt+1,meshType='Progression',coef=1/rj2)
gm.model.geo.mesh.setTransfiniteCurve(lcc,N2r)
gm.model.geo.mesh.setTransfiniteCurve(lco,N2r)
gm.model.geo.mesh.setTransfiniteCurve(lcr,N2z)
gm.model.geo.mesh.setTransfiniteCurve(lcl,N2z)
gm.model.geo.synchronize()

#=====> Surfaces
lm=gm.model.geo.add_curve_loop([lma,lml,lmc,lm1,lmr,lmb,lm2,lmo]) #===> Main Jet
sm=gm.model.geo.add_plane_surface([lm])
lp=gm.model.geo.add_curve_loop([lpa,lpl,lp1,lpc,lp3,lpr,lpb,lp4,lpo,-lp2]) #===> Pilot Jet
sp=gm.model.geo.add_plane_surface([lp])
lc=gm.model.geo.add_curve_loop([lca,lcl,lc1,lcc,lcr,lcb,lco,-lc2]) #===> Coflow
sc=gm.model.geo.add_plane_surface([lc])
lo=gm.model.geo.add_curve_loop([lol,-lmo,-lm2,lme,lp2,-lpo,-lp4,lpe,lc2,-lco,lor,loo]) #===> Outlet
so=gm.model.geo.add_plane_surface([lo])
gm.model.geo.synchronize()

#=====> transfinite
gm.model.geo.mesh.setTransfiniteSurface(sm,cornerTags=[Points[ 0],Points[ 2],Points[ 4],Points[ 6]])
gm.model.geo.mesh.setTransfiniteSurface(sp,cornerTags=[Points[ 8],Points[11],Points[14],Points[16]])
gm.model.geo.mesh.setTransfiniteSurface(sc,cornerTags=[Points[18],Points[21],Points[23],Points[25]])
gm.model.geo.synchronize()

#=====> Recombine
gm.model.geo.mesh.setRecombine(gdim,sm)
gm.model.geo.mesh.setRecombine(gdim,sp)
gm.model.geo.mesh.setRecombine(gdim,sc)
gm.model.geo.synchronize()

#=====> Field
gm.model.mesh.field.add("Box", 1) # Flame center
gm.model.mesh.field.setNumber(1, "VIn" , hfc )
gm.model.mesh.field.setNumber(1, "VOut", hf  )
gm.model.mesh.field.setNumber(1, "YMin", 0  )
gm.model.mesh.field.setNumber(1, "YMax", r0 )
gm.model.mesh.field.setNumber(1, "XMin", y2   )
gm.model.mesh.field.setNumber(1, "XMax", y2+Lr)
gm.model.mesh.field.setNumber(1, "Thickness",  1e-1)
gm.model.mesh.field.add("Box", 2) # Flame side
gm.model.mesh.field.setNumber(2, "VIn" , hfs )
gm.model.mesh.field.setNumber(2, "VOut", hf  )
gm.model.mesh.field.setNumber(2, "YMin", r0   )
gm.model.mesh.field.setNumber(2, "YMax", r0+Dr)
gm.model.mesh.field.setNumber(2, "XMin", y2   )
gm.model.mesh.field.setNumber(2, "XMax", y2+Lr)
gm.model.mesh.field.setNumber(2, "Thickness",  1e-1)
gm.model.mesh.field.add("Box", 3) # Flame side
gm.model.mesh.field.setNumber(3, "VIn" , h00 )
gm.model.mesh.field.setNumber(3, "VOut", hf  )
gm.model.mesh.field.setNumber(3, "YMin", r0-rp)
gm.model.mesh.field.setNumber(3, "YMax", r0+rp)
gm.model.mesh.field.setNumber(3, "XMin", y2   )
gm.model.mesh.field.setNumber(3, "XMax", y2+Lf)
gm.model.mesh.field.setNumber(3, "Thickness",0.2)
# gm.model.mesh.field.add("BoundaryLayer", 3)
# gm.model.mesh.field.setNumber(3, "AnisoMax"  , 1e3 )
# gm.model.mesh.field.setNumber(3, "Quads"     , 1   )
# gm.model.mesh.field.setNumber(3, "Thickness" , 1   )
# gm.model.mesh.field.setNumber(3, "CurvesList", lmr )
# gm.model.mesh.field.setNumber(3, "NbLayers"  , 5   )
# gm.model.mesh.field.setNumber(3, "Ratio"     , 1.1 )
# gm.model.mesh.field.setNumber(3, "Size"      , h0*rb0 )
# gm.model.mesh.field.setNumber(3, "SizeFar"   , h0  )
gm.model.mesh.field.add("Min", 4)
gm.model.mesh.field.setNumbers(4, "FieldsList", [1,2,3])
gm.model.mesh.field.setAsBackgroundMesh(4)

#===================================================================================
util.Section( 'Physical groups : {:.3f} s'.format(time.time()-t0),1,5,'r' )
#===================================================================================

gm.model.addPhysicalGroup(1,[lma,lml,lol]    ,tag=1000,name='Axis'        )
gm.model.addPhysicalGroup(1,[lmc,lm1]        ,tag=2001,name='Inlet-Main'  )
gm.model.addPhysicalGroup(1,[lpc,lp1,lp3]    ,tag=2002,name='Inlet-Pilot' )
gm.model.addPhysicalGroup(1,[lcc,lc1]        ,tag=2003,name='Inlet-Coflow')
gm.model.addPhysicalGroup(1,[lmr,lmb]        ,tag=3001,name='Wall-Main'   )
gm.model.addPhysicalGroup(1,[lpa,lpl,lpr,lpb],tag=3002,name='Wall-Pilot'  )
gm.model.addPhysicalGroup(1,[lca,lcl]        ,tag=3003,name='Wall-Coflow' )
gm.model.addPhysicalGroup(1,[lme]            ,tag=4001,name='Leap-Main'   )
gm.model.addPhysicalGroup(1,[lpe]            ,tag=4002,name='Leap-Pilot'  )
gm.model.addPhysicalGroup(1,[lcr,lcb,lor]    ,tag=5001,name='External'    )
gm.model.addPhysicalGroup(1,[loo]            ,tag=5002,name='Outlet'      )
gm.model.addPhysicalGroup(2,[sm,sp,sc,so]    ,tag=6000,name='Fluid'       )

#===================================================================================
util.Section( 'Meshing : {:.3f} s'.format(time.time()-t0),1,5,'r' )
#===================================================================================

gm.option.setNumber("Mesh.Algorithm", 7)
gm.option.setNumber("Mesh.SmoothRatio",1)
gm.option.setNumber('Mesh.OptimizeThreshold',0.6)
gm.option.setNumber('Mesh.OptimizeNetgen',1)
gm.option.setNumber("Mesh.SaveElementTagType",2)
gm.option.setNumber("Mesh.MeshSizeExtendFromBoundary",0)
gm.option.setNumber("Mesh.MeshSizeFromPoints",0)
gm.option.setNumber("Mesh.MeshSizeFromCurvature",0)
# gm.model.mesh.setSmoothing(gdim,so,100)
gm.model.mesh.generate(gdim)
gm.model.mesh.optimize("Netgen")

#===================================================================================
util.Section( 'Writing : {:.3f} s'.format(time.time()-t0),1,5,'r' )
#===================================================================================

gm.model.geo.synchronize()
gm.write(d_mesh+name+'-Geo.geo_unrolled')
gm.write(d_mesh+name+'-MESH.msh')
gm.write(d_mesh+name+'-MESH.bdf')

gm.finalize()

if CONV :
	#===================================================================================
	util.Section( 'Conversion : {:.3f} s'.format(time.time()-t0),1,5,'r' )
	#===================================================================================
	mesh=mo.gmsh.read(d_mesh+name+'-MESH.msh')
	mo.ansys.write(d_mesh+name+'-ANSYS.msh',mesh)
	# mo.cgns.write(d_mesh+name+'-MESH.cgns',mesh)