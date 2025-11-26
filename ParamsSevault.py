#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#=====> Boudary conditions
BC_f={'T':300 , 'V':98.2  , 'H2':0.55 , 'O2':0    ,'N2':0 , 'H2O':0 , 'CH4':0.45 , 'CO2':0    , 'unit':'X'} # 55p H2 , Re 15k
BC_o={'T':300 , 'V':0.778 , 'H2':0    , 'O2':0.32 ,'N2':0 , 'H2O':0 , 'CH4':0    , 'CO2':0.68 , 'unit':'X'}

#=====> Dimensions
D0=  5e-3
D1=  6e-3
D2=350e-3
Lt=1
Lc=1e-2

ep=D1-D0

#====================> Helium content 
h2=55
re=15

if   re==15 and h2==55 : flame='a' ; Umoy=98.2
elif re==15 and h2==45 : flame='b' ; Umoy=84.5
elif re==15 and h2==37 : flame='c' ; Umoy=75.8
elif re==12 and h2==55 : flame='d' ; Umoy=78.6
elif re==18 and h2==55 : flame='e' ; Umoy=117.8

#====================> Directories
dir0='/mnt/d/Python/Sandia/'
dirs=dir0+'DATA-Oxy/Mean_RMS/'
# dirc='/mnt/d/FLUENT/Sandia-Sevault/RUN-00-COLD/DUMP-05-Side/'
# dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Sevault/RUN-00-55p-15k/DUMP-03-Tw1400/'
# dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Sevault/RUN-00-55p-15k/DUMP-04-DF3/'
# dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Sevault/RUN-00-55p-15k/DUMP-04-Tprof/'
# dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Sevault/RUN-00-55p-15k/DUMP-05-ke/'
# dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Sevault/RUN-00-55p-15k/DUMP-05-DF4/'
# dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Sevault/RUN-00-55p-15k/DUMP-05-Tw400-Cmix/'
dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Sevault/RUN-00-55p-15k/DUMP/'
dird=dirc+'DATA/'
dirp=dirc+'PLOT/'
slice='Data-all.dat'
par=''
# par='1k'
# par='ko'
# par='C019'
# par='DF30'
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