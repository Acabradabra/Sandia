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
(Sysa,NSysa,Arg)=util.Parseur(['BD'],0,'Arg : ')
(                             [ BD ])=Arg

from numpy import *
import sys
import time
import pygmsh as pg
import gmsh   as gm
import Meshing as me
import meshio as mo

import matplotlib.pyplot as plt
import matplotlib

from ParamsGarnier import *

t0=time.time()
#%%=================================================================================
#                     Parameters
#===================================================================================

gdim = 2  # Geometric dimension of the mesh

#=====> files
d_mesh='MESH-GARNIER/'
d_plot='PLOT/'
name='Sandia-Garnier'

#=====> Mesh Size
h0=1e-4
h1=2e-4
he=1e-3
hf=2e-3

hfc=2e-4
hfs=5e-4

#=====> Flame refinement params
Lr=0.2  # Refinement length
Lf=4*ep # Refinement foot
Dr=1e-2

#=====> Boundary layer params
rb=1.1 # Ratio slices boundary layer
rt=1   # Ratio transition
rs=1   # Ratio surface caré/triangle
ra=0.7 # Aspect ratio last slice
rj0=0.45
N0=10
N1=10

#%%=================================================================================
util.Section( 'Geometrie : {:.3f} s'.format(time.time()-t0),1,5,'r' )
#===================================================================================

#=====> X dim
r0=0.5*D0
r1=0.5*D1
r2=0.5*D2

#=====> Boundary layer
rat=0.25*sqrt(3)*rs ; print( '=> Aspect Triangle=Square : {:.3f}'.format(rat) )
(h00,Lt0)=me.h_smooth(rb,N0,rt*h0) ; print('=> Fuel   : h0 {:.2e} , Lt {:.2e} , ra*h/h0 {:.2e}'.format(h00,Lt0,ra*h0/h00))
(h10,Lt1)=me.h_smooth(rb,N1,rt*h1) ; print('=> Coflow : h0 {:.2e} , Lt {:.2e} , ra*h/h0 {:.2e}'.format(h10,Lt1,ra*h1/h10))
N0r=int(round( (r0-    Lt0 )/h0 ,0))+1
N1r=int(round( (r2-(r1+Lt1))/h1 ,0))+1
(Nt,rj)=me.Nbl2(rj0,ra,rb) ; print('=> Transition : N {}  ,  r2 {:.3f}'.format(Nt,rj))
Ltj0=me.ht(h0*rj0,rj,Nt) ; print('=> Fuel   : L transition jet {:.3f} [mm]'.format(Ltj0*1e3))
Ltj1=me.ht(h1*rj0,rj,Nt) ; print('=> Coflow : L transition jet {:.3f} [mm]'.format(Ltj1*1e3))
N0z=int(round((Lc-Ltj0)/(h0*ra),0))+1
N1z=int(round((Lc-Ltj1)/(h1*ra),0))+1

BD=False
if BD : sys.exit('=> Boundary layer')

gm.initialize()
#=====> Points
Points=[
gm.model.geo.add_point( Lc     ,0      , 0 , meshSize=h0  ), #===> Fuel
gm.model.geo.add_point( Lc-Ltj0,0      , 0 , meshSize=h0  ),
gm.model.geo.add_point(  0     ,0      , 0 , meshSize=h0  ),
gm.model.geo.add_point(  0     ,r0-Lt0 , 0 , meshSize=h0  ),
gm.model.geo.add_point(  0     ,r0     , 0 , meshSize=h00 ),
gm.model.geo.add_point( Lc-Ltj0,r0     , 0 , meshSize=h00 ),
gm.model.geo.add_point( Lc     ,r0     , 0 , meshSize=h00 ),
gm.model.geo.add_point( Lc     ,r0-Lt0 , 0 , meshSize=h0  ),
gm.model.geo.add_point( Lc     ,r1     , 0 , meshSize=h10 ), #===> Coflow
gm.model.geo.add_point( Lc     ,r1+Lt1 , 0 , meshSize=h1  ),
gm.model.geo.add_point( Lc-Ltj1,r1+Lt1 , 0 , meshSize=h1  ),
gm.model.geo.add_point( Lc-Ltj1,r1     , 0 , meshSize=h10 ),
gm.model.geo.add_point(  0     ,r1     , 0 , meshSize=h10 ),
gm.model.geo.add_point(  0     ,r1+Lt1 , 0 , meshSize=h1  ),
gm.model.geo.add_point(  0     ,r2     , 0 , meshSize=he  ),
gm.model.geo.add_point( Lc-Ltj1,r2     , 0 , meshSize=he  ),
gm.model.geo.add_point( Lt     ,r2     , 0 , meshSize=hf  ), #===> Outlet
gm.model.geo.add_point( Lt     ,0      , 0 , meshSize=hf  ) ]
Np=len(Points)
gm.model.geo.synchronize()

#=====> Lines Fuel
lma=gm.model.geo.add_line(Points[ 0],Points[ 1]) #===> Fuel left transition
lml=gm.model.geo.add_line(Points[ 1],Points[ 2]) #===> Fuel left
lmc=gm.model.geo.add_line(Points[ 2],Points[ 3]) #===> Fuel Inlet flow
lm1=gm.model.geo.add_line(Points[ 3],Points[ 4]) #===> Fuel Inlet BL
lmr=gm.model.geo.add_line(Points[ 4],Points[ 5]) #===> Fuel right
lmb=gm.model.geo.add_line(Points[ 5],Points[ 6]) #===> Fuel right transition
lm2=gm.model.geo.add_line(Points[ 6],Points[ 7]) #===> Fuel Outlet BL
lmo=gm.model.geo.add_line(Points[ 7],Points[ 0]) #===> Fuel Outlet flow
lme=gm.model.geo.add_line(Points[ 6],Points[ 8]) #===> Fuel leap
#=====> Lines Coflow
lc2=gm.model.geo.add_line(Points[ 8],Points[ 9]) #===> Coflow Outlet BL 
lca=gm.model.geo.add_line(Points[ 8],Points[11]) #===> Coflow Left transition
lcl=gm.model.geo.add_line(Points[11],Points[12]) #===> Coflow Left
lc1=gm.model.geo.add_line(Points[12],Points[13]) #===> Coflow Inlet BL
lcc=gm.model.geo.add_line(Points[13],Points[14]) #===> Coflow Inlet
lcr=gm.model.geo.add_line(Points[14],Points[15]) #===> Coflow right
lci=gm.model.geo.add_line(Points[13],Points[10]) #===> Coflow outlet of BL
lcb=gm.model.geo.add_line(Points[10],Points[ 9]) #===> Coflow outlet of BL Transition
#=====> Lines Outlet
lor=gm.model.geo.add_line(Points[15],Points[16]) #===> Outlet right
loo=gm.model.geo.add_line(Points[16],Points[17]) #===> Outlet
lol=gm.model.geo.add_line(Points[17],Points[ 0]) #===> Outlet left
gm.model.geo.synchronize()

#=====> transfinite
gm.model.geo.mesh.setTransfiniteCurve(lm1,N0+1,meshType='Progression',coef=1/rb) #===> Main Jet
gm.model.geo.mesh.setTransfiniteCurve(lm2,N0+1,meshType='Progression',coef=  rb)
gm.model.geo.mesh.setTransfiniteCurve(lma,Nt+1,meshType='Progression',coef=  rj)
gm.model.geo.mesh.setTransfiniteCurve(lmb,Nt+1,meshType='Progression',coef=1/rj)
gm.model.geo.mesh.setTransfiniteCurve(lmc,N0r)
gm.model.geo.mesh.setTransfiniteCurve(lmo,N0r)
gm.model.geo.mesh.setTransfiniteCurve(lmr,N0z)
gm.model.geo.mesh.setTransfiniteCurve(lml,N0z)
gm.model.geo.mesh.setTransfiniteCurve(lc1,N1+1,meshType='Progression',coef=rb) #===> Coflow
gm.model.geo.mesh.setTransfiniteCurve(lc2,N1+1,meshType='Progression',coef=rb)
gm.model.geo.mesh.setTransfiniteCurve(lca,Nt+1,meshType='Progression',coef=  rj)
gm.model.geo.mesh.setTransfiniteCurve(lcb,Nt+1,meshType='Progression',coef=1/rj)
# gm.model.geo.mesh.setTransfiniteCurve(lcc,N1r)
gm.model.geo.mesh.setTransfiniteCurve(lci,N1z)
gm.model.geo.mesh.setTransfiniteCurve(lcl,N1z)
gm.model.geo.synchronize()

#=====> Surfaces
lm=gm.model.geo.add_curve_loop([lma,lml,lmc,lm1,lmr,lmb,lm2,lmo]) #===> Main Jet
sm=gm.model.geo.add_plane_surface([lm])
lc=gm.model.geo.add_curve_loop([lca,lcl,lc1,lci,lcb,-lc2]) #===> Coflow BL
sc=gm.model.geo.add_plane_surface([lc])
lo=gm.model.geo.add_curve_loop([lol,-lmo,-lm2,lme,lc2,-lcb,-lci,lcc,lcr,lor,loo]) #===> Outlet + Coflow
so=gm.model.geo.add_plane_surface([lo])
gm.model.geo.synchronize()

#=====> transfinite
gm.model.geo.mesh.setTransfiniteSurface(sm,cornerTags=[Points[ 0],Points[ 2],Points[ 4],Points[ 6]]) # Fuel
gm.model.geo.mesh.setTransfiniteSurface(sc,cornerTags=[Points[ 8],Points[12],Points[13],Points[ 9]]) # Coflow BL
gm.model.geo.synchronize()

#=====> Recombine
gm.model.geo.mesh.setRecombine(gdim,sm)
gm.model.geo.mesh.setRecombine(gdim,sc)
gm.model.geo.synchronize()

#=====> Fields
gm.model.mesh.field.add("Box", 1) # Flame center
gm.model.mesh.field.setNumber(1, "VIn" , hfc )
gm.model.mesh.field.setNumber(1, "VOut", hf  )
gm.model.mesh.field.setNumber(1, "YMin", 0  )
gm.model.mesh.field.setNumber(1, "YMax", r0 )
gm.model.mesh.field.setNumber(1, "XMin", Lc   )
gm.model.mesh.field.setNumber(1, "XMax", Lc+Lr)
gm.model.mesh.field.setNumber(1, "Thickness",  2e-2)
gm.model.mesh.field.add("Box", 2) # Flame side
gm.model.mesh.field.setNumber(2, "VIn" , hfs )
gm.model.mesh.field.setNumber(2, "VOut", hf  )
gm.model.mesh.field.setNumber(2, "YMin", r0   )
gm.model.mesh.field.setNumber(2, "YMax", r0+Dr)
gm.model.mesh.field.setNumber(2, "XMin", Lc   )
gm.model.mesh.field.setNumber(2, "XMax", Lc+Lr)
gm.model.mesh.field.setNumber(2, "Thickness",  2e-2)
gm.model.mesh.field.add("Box", 3) # Flame side
gm.model.mesh.field.setNumber(3, "VIn" , h10 )
gm.model.mesh.field.setNumber(3, "VOut", hf  )
gm.model.mesh.field.setNumber(3, "YMin", r0  )
gm.model.mesh.field.setNumber(3, "YMax", r1  )
gm.model.mesh.field.setNumber(3, "XMin", Lc   )
gm.model.mesh.field.setNumber(3, "XMax", Lc+Lf)
gm.model.mesh.field.setNumber(3, "Thickness",1e-2)
gm.model.mesh.field.add("Box", 4) # Coflow
gm.model.mesh.field.setNumber(4, "VIn" , h1 )
gm.model.mesh.field.setNumber(4, "VOut", hf  )
gm.model.mesh.field.setNumber(4, "YMin", r1  )
gm.model.mesh.field.setNumber(4, "YMax", r2*0.1)
gm.model.mesh.field.setNumber(4, "XMin", 0   )
gm.model.mesh.field.setNumber(4, "XMax", Lc  )
gm.model.mesh.field.setNumber(4, "Thickness",2e-2)
gm.model.mesh.field.add("Min", 5)
gm.model.mesh.field.setNumbers(5, "FieldsList", [1,2,3,4])
gm.model.mesh.field.setAsBackgroundMesh(5)
gm.model.geo.synchronize()

#===================================================================================
util.Section( 'Physical groups : {:.3f} s'.format(time.time()-t0),1,5,'r' )
#===================================================================================

gm.model.addPhysicalGroup(1,[lma,lml,lol],tag=1000,name='Axis'        )
gm.model.addPhysicalGroup(1,[lmc,lm1]    ,tag=2001,name='Inlet-Main'  )
gm.model.addPhysicalGroup(1,[lcc,lc1]    ,tag=2003,name='Inlet-Coflow')
gm.model.addPhysicalGroup(1,[lmr,lmb]    ,tag=3001,name='Wall-Main'   )
gm.model.addPhysicalGroup(1,[lca,lcl]    ,tag=3003,name='Wall-Coflow' )
gm.model.addPhysicalGroup(1,[lme]        ,tag=4001,name='Leap-Main'   )
gm.model.addPhysicalGroup(1,[lcr,lor]    ,tag=5001,name='External'    )
gm.model.addPhysicalGroup(1,[loo]        ,tag=5002,name='Outlet'      )
gm.model.addPhysicalGroup(2,[sm,sc,so]   ,tag=6000,name='Fluid'       )

#===================================================================================
util.Section( 'Meshing : {:.3f} s'.format(time.time()-t0),1,5,'r' )
#===================================================================================

# gm.option.setNumber("Mesh.Algorithm", 7)
gm.option.setNumber("Mesh.Algorithm", 6)
# gm.option.setNumber("Mesh.Smoothing",10)
gm.option.setNumber("Mesh.SmoothRatio",1)
gm.option.setNumber('Mesh.OptimizeThreshold',0.6)
gm.option.setNumber('Mesh.OptimizeNetgen',1)
gm.option.setNumber("Mesh.SaveElementTagType",2)
gm.option.setNumber("Mesh.MeshSizeExtendFromBoundary",1)
gm.option.setNumber("Mesh.MeshSizeFromParametricPoints",1 )
gm.option.setNumber("Mesh.MeshSizeFromPoints",0 )
gm.option.setNumber("Mesh.MeshSizeFromCurvature",0)
gm.model.mesh.setSmoothing(gdim,so,100)
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
# %%
