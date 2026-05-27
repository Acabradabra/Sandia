#!/usr/bin/env python3
# -*- coding: utf-8 -*-

flame='a'

#=====> Dimensions
D0=  5e-3
D1=  6e-3
D2=350e-3
Lt=1
Lc=1e-2

ep=D1-D0

#====================> Directories IN
dir0='/mnt/beegfs/ZEUS/'
dirs=dir0+'Python/Sandia/DATA-Oxy/Mean_RMS/'
dirc=dir0+'/FLUENT/Sandia-Sevault/COMPA-A-EDC-EDM/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-01-EDC/DUMP/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-01-EDC/DUMP-00-A-EDC/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-01-EDC/DUMP-01-A-Dyn/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-01-EDC/DUMP-01-A-Cv3-Ct03/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-01-EDC/DUMP-02-B-EDC/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-01-EDC/DUMP-00-D-EDC/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-01-EDC/DUMP-02-C-EDC/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-01-EDC/DUMP-02-F-Vpol/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-01-EDC/DUMP-03-F-Dyn/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-02-EDM/DUMP/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-02-EDM/DUMP-03-A-CO/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-02-EDM/DUMP-03-B-CO/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-02-EDM/DUMP-03-C-CO/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-02-EDM/DUMP-03-F-CO/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-03-PDF/DUMP/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-03-PDF/DUMP-01-A-SLF/'
# dirc=dir0+'/FLUENT/Sandia-Sevault/RUN-03-PDF/DUMP-03-A-FGM/'
if   '-A-' in dirc : flame='a'
elif '-B-' in dirc : flame='b'
elif '-C-' in dirc : flame='c'
elif '-D-' in dirc : flame='d'
elif '-E-' in dirc : flame='e'
elif '-F-' in dirc : flame='f'

#====================> Case
if   flame=='a' : re=15 ; h2=55  ; Umoy=98.2  ; Uox=0.778
elif flame=='b' : re=15 ; h2=45  ; Umoy=84.5  ; Uox=0.755
elif flame=='c' : re=15 ; h2=37  ; Umoy=75.8  ; Uox=0.739
elif flame=='d' : re=12 ; h2=55  ; Umoy=78.6  ; Uox=0.622
elif flame=='e' : re=18 ; h2=55  ; Umoy=117.8 ; Uox=0.933
elif flame=='f' : re=15 ; h2=100 ; Umoy=198.7 ; Uox=0.877
ch4=100-h2

#=====> Boundary conditions
BC_f={'T':300 , 'V':Umoy , 'H2':h2/100 , 'O2':0    ,'N2':0 , 'H2O':0 , 'CH4':ch4/100 , 'CO2':0    , 'unit':'X'}
BC_o={'T':300 , 'V':Uox  , 'H2':0      , 'O2':0.32 ,'N2':0 , 'H2O':0 , 'CH4':0       , 'CO2':0.68 , 'unit':'X'}


#====================> Directories OUT
dird=dirc+'DATA/'
dirp=dirc+'PLOT/'
slice='Data-all.dat'
par=''
if par :
    slice='Data-%s.dat'%(par)
    dirp=dirc+'PLOT-%s/'%(par)

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
# 'T'   :'Tpg' ,
'Tpg' :'Tpg'  ,
'Rho' :'dens' ,
'mix' :'phi'  ,
'mixH':'FH' ,
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