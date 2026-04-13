#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#=====> Dimensions
D0=  5e-3
D1=  6e-3
D2=350e-3
Lt=1
Lc=1e-2

ep=D1-D0

#====================> Case
flame='d'

if   flame=='a' : re=15 ; h2=55  ; Umoy=98.2  ; Uox=0.778
elif flame=='b' : re=15 ; h2=45  ; Umoy=84.5  ; Uox=0.755
elif flame=='c' : re=15 ; h2=37  ; Umoy=75.8  ; Uox=0.739
elif flame=='d' : re=12 ; h2=55  ; Umoy=78.6  ; Uox=0.622
elif flame=='e' : re=18 ; h2=55  ; Umoy=117.8 ; Uox=0.933
ch4=100-h2

#=====> Boundary conditions
BC_f={'T':300 , 'V':Umoy , 'H2':h2/100 , 'O2':0    ,'N2':0 , 'H2O':0 , 'CH4':ch4/100 , 'CO2':0    , 'unit':'X'}
BC_o={'T':300 , 'V':Uox  , 'H2':0      , 'O2':0.32 ,'N2':0 , 'H2O':0 , 'CH4':0       , 'CO2':0.68 , 'unit':'X'}

#====================> Directories
dir0='/mnt/d/Python/Sandia/'
dirs=dir0+'DATA-Oxy/Mean_RMS/'
# dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Sevault/RUN-01-CADFEM/DUMP/'
dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Sevault/RUN-01-CADFEM/DUMP-00-D-EDC/'
dird=dirc+'DATA/'
dirp=dirc+'PLOT/'
slice='Data-all.dat'
# slice='Data-all-Reac.dat'
par=''
# par='DO-4433'
if par :
    slice='Data-%s.dat'%(par)
    dirp=dirc+'PLOT-%s/'%(par)

# D_compa=[
# ]

#====================> Variables
Vars_TNF=[
'CO2'  ,
'O2'   ,
'CO'   ,
'N2'   ,
'CH4'  ,
'H2O'  ,
'H2'   ,
'COLIF',
'F560' ,
'b3'   ,
'Tray' ,
'Tpg'  ,
'phi'  ,
'xmm'  ,
'shot' ,
'dens' ,
'FH'   ,
'FC'
]

#====================> Correspondences TNF variables
Cor={
# 'r'   :'xmm'  ,
'r'   :'shot'  ,
'T'   :'Tray' ,
'Tpg' :'Tpg'  ,
'Rho' :'dens' ,
'mix' :'phi'  ,
'mixH':'FH'   ,
'mixC':'FC'   ,
'o2'  :'O2'   ,
'h2'  :'H2'   ,
'n2'  :'N2'   ,
'co'  :'CO'   ,
'ch4' :'CH4'  ,
'h2o' :'H2O'  ,
'co2' :'CO2'
}

#====================> Uncertainties
Err={
'r'   :0 ,
'T'   :0.02 ,
'Tpg' :0.02 ,
'Rho' :0 ,
'mix' :0 ,
'mixH':0 ,
'mixC':0 ,
'o2'  :0 ,
'h2'  :0.08 ,
'n2'  :0.02 ,
'co'  :0.08 ,
'ch4' :0 ,
'h2o' :0 ,
'co2' :0
}