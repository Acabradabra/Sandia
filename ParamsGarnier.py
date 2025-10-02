#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#=====> Boudary conditions
BC_f={'T':298 , 'V':296 , 'H2':1 , 'O2':0     ,'N2':0      , 'H2O':0      }
BC_o={'T':296 , 'V':1   , 'H2':0 , 'O2':0.2303,'N2':0.7625 , 'H2O':0.0072 }

#=====> Dimensions
D0=  3.75e-3
D1=  5e-3
# D2=336e-3
# Lt= 1
D2=100e-3
Lt=0.5
Lc= 1e-2

ep=D1-D0