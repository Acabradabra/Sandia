#!/usr/bin/env python3
# -*- coding: utf-8 -*-

flame='D'

#=====> Mixture fraction
BC_m={'T':294 ,'CH4':0.156,'O2':0.197,'CO2':0    ,'H2O':0     ,'CO':0,'NO':0,'N2':0.647 }
BC_p={'T':1880,'CH4':0    ,'O2':0.054,'CO2':0.111,'H2O':0.0942,'CO':0,'NO':0,'N2':0.7408}
BC_o={'T':291 ,'CH4':0    ,'O2':0.230,'CO2':0    ,'H2O':0     ,'CO':0,'NO':0,'N2':0.77  }

#=====> Dimensions
D0= 7.2e-3
D1=18.2e-3
# D2=50.0e-3
# D2=0.1
D2=0.3
# Lt= 0.5
Lt=1
Lc=2e-2
Ls=1e-2
ep= 1.0e-3

#=====> Flow
Um=49.6
Up=11.4
Uc=0.9

Umoy=Um

#====================> Directories IN
# dir0='/mnt/scratch/ZEUS/'
dir0='/mnt/beegfs/ZEUS/'
dirs=dir0+'Python/Sandia/DATA-Pilote/'
dirc=dir0+'/FLUENT/Sandia-Jaravel/RUN-D300-00/COMPA-BFER-Laera/'
# dirc=dir0+'/FLUENT/Sandia-Jaravel/RUN-D300-00/DUMP/'
# dirc=dir0+'/FLUENT/Sandia-Jaravel/RUN-D300-00/DUMP-01-EDM/'
# dirc=dir0+'/FLUENT/Sandia-Jaravel/RUN-D300-00/DUMP-02-FRED/'
# dirc=dir0+'/FLUENT/Sandia-Jaravel/RUN-D300-00/DUMP-03-EDC-BFER/'
# dirc=dir0+'/FLUENT/Sandia-Jaravel/RUN-D300-00/DUMP-04-EDC-Laera/'

#====================> Directories OUT
dird=dirc+'DATA/'
dirp=dirc+'PLOT/'
slice='Data-all.dat'
par=''
if par :
    slice='Data-%s.dat'%(par)
    dirp=dirc+'PLOT-%s/'%(par)

#====================> Correspondences TNF variables
Cor={
'r'   :'r/d'  ,
'T'   :'T(K)' ,
'mix' :'F'    ,
'o2'  :'YO2'  ,
'oh'  :'YOH'  ,
'no'  :'YNO'  ,
'h2'  :'YH2'  ,
'n2'  :'YN2'  ,
'co'  :'YCO'  ,
'ch4' :'YCH4' ,
'h2o' :'YH2O' ,
'co2' :'YCO2'
}

#====================> Uncertainties
Err={
'r'   : 0 ,
'T'   : 0.03 ,
'k'   : 0 ,
'Vel' : 0 ,
'mix' : 0 ,
'o2'  : 0 ,
'oh'  : 0.10 ,
'no'  : 0.15 ,
'h2'  : 0.12 ,
'n2'  : 0.03 ,
'co'  : 0.20 ,
'ch4' : 0 ,
'h2o' : 0.04 ,
'co2' : 0.04
}