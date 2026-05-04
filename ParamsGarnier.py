#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#=====> Boudary conditions
# BC_f={'T':298 , 'V':296 , 'H2':1    , 'O2':0     ,'N2':0      , 'H2O':0      , 'CH4':0 , 'He':0.57 } #  0p He
BC_f={'T':298 , 'V':296 , 'H2':0.43 , 'O2':0     ,'N2':0      , 'H2O':0      , 'CH4':0 , 'He':0.57 } # 40p He
BC_o={'T':296 , 'V':1   , 'H2':0    , 'O2':0.2303,'N2':0.7625 , 'H2O':0.0072 , 'CH4':0 , 'He':0    }

#=====> Dimensions
D0=  3.75e-3
D1=  5e-3
# D2=0.1
# Lt=0.5
D2=336e-3
Lt=1
Lc= 1e-2

ep=D1-D0

#=====> Lines positions
Pos=[1/8,1/4,3/8,1/2,5/8,3/4,1]

#====================> Helium content 
# he='he4'
he=0.4
if   he==0   : Ld=180 ; Umoy=296
elif he==0.2 : Ld=150 ; Umoy=294
elif he==0.4 : Ld=100 ; Umoy=256

#====================> Directories
dir0='/mnt/d/Python/Sandia/'
dirv=dir0+'DATA-H2/ETHZ_H2/%.0fhe/'%(he*100)
dirs=dir0+'DATA-H2/SANDH2_A/he%.0f/he%.0fstatY/'%(he*10,he*10)
# dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-02-Big-40p/DUMP/'
# dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-02-Big-40p/DUMP-00-He-FD3/'
dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-02-Big-40p/DUMP-01-He-Uprof2/'
# dirc='/mnt/scratch/ZEUS/FLUENT/Sandia-Garnier/RUN-02-Big-40p/DUMP-04-He-FD39/'
dird=dirc+'DATA/'
dirp=dirc+'PLOT/'

#====================> Velocity plots
Rlims=[3,8,10,20,20,20,40,40,40]
Pv=[ '0','1/16','1/8','1/4','3/8','1/2','5/8','3/4','1' ]

#====================> Correspondences TNF variables
Cor={
'r'   :'r(mm)',
'T'   :'T(K)' ,
'mixH':'Fblgr',
'o2'  :'YO2'  ,
'h2'  :'YH2'  ,
'n2'  :'YN2'  ,
'h2o' :'YH2O'
}

#====================> Uncertainties
Err={
'T'   :0.03,
'mixH':0.04,
'o2'  :0.04,
'h2'  :0.04,
'n2'  :0.04,
'h2o' :0.04
}
