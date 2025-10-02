#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#=====> Mixture fraction
Y_m={'T':294 ,'CH4':0.156,'O2':0.197,'CO2':0    ,'H2O':0     ,'CO':0,'NO':0,'N2':0.647 }
Y_p={'T':1880,'CH4':0    ,'O2':0.054,'CO2':0.111,'H2O':0.0942,'CO':0,'NO':0,'N2':0.7408}
Y_o={'T':291 ,'CH4':0    ,'O2':0.230,'CO2':0    ,'H2O':0     ,'CO':0,'NO':0,'N2':0.77  }

#=====> Dimensions
D0= 7.2e-3
D1=18.2e-3
# D2=50.0e-3
D2=0.1
Lt= 0.5
Lc=2e-2
Ls=1e-2
ep= 1.0e-3

#=====> Flow
Um=49.6
Up=11.4
Uc=0.9

if False :
    import Fluent as fl
    
    Yc_m=fl.Yc(Ych4[0],Yco2[0],Yco[0],Mol_m)
    Yc_p=fl.Yc(Ych4[1],Yco2[1],Yco[1],Mol_m)
    Yc_o=fl.Yc(Ych4[2],Yco2[2],Yco[2],Mol_m)

    Zp_f=(Yc_p-Yc_o)/(Yc_m-Yc_o) ; print('=> {:.5f}'.format(Zp_f))
    Zo_s=(Yc_o-Yc_p)/(Yc_m-Yc_p) ; print('=> {:.5f}'.format(Zo_s))