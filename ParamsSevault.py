#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#=====> Boudary conditions
BC_f={'T':300 , 'V':98.2  , 'H2':0.55 , 'O2':0    ,'N2':0 , 'H2O':0 , 'CH4':0.45 , 'CO2':0    } # 55p H2 , Re 15k
BC_o={'T':300 , 'V':0.778 , 'H2':0    , 'O2':0.32 ,'N2':0 , 'H2O':0 , 'CH4':0    , 'CO2':0.68 }

#=====> Dimensions
D0=  5e-3
D1=  6e-3
D2=350e-3
Lt=1
Lc=1e-2

ep=D1-D0

#=====> Lines positions
Pos=[1/8,1/4,3/8,1/2,5/8,3/4,1]