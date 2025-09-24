#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Utilities as util
(Sysa,NSysa,Arg)=util.Parseur([],0,'Arg : ')
(                             [])=Arg

from numpy import *
import sys
import time
import Fluent as fl

t0=time.time()
(plt,mtp)=util.Plot0()
#===================================================================================
#                     Parameters
#===================================================================================

file_in='/mnt/scratch/PRECIZE/Sandia-Jaravel/RUN-D100-02-Laera/INIT-Fluent/Interp-00-Fluent.ip'
file_ou=file_in[:-3]+'-New.ip'

#===================================================================================
util.Section('=> Reading : {:.3f}'.format(time.time()-t0),2,5,'r')
#===================================================================================

fin=open(file_in,'r')
Lfin=fin.readlines()
fin.close()

Nd=int(float(Lfin[1]))
Np=int(float(Lfin[2]))
Nv=int(float(Lfin[3]))
Nl=len(Lfin)

Fields=[ f[:-1] for f in Lfin[4:4+Nv]]

print( '=> Nd={:d}, Np={:d}, Nv={:d}, Nl={:d}'.format(Nd,Np,Nv,Nl) )
print( '=> Fields : {}'.format(Fields) )

DATA=zeros((Np,Nv+Nd))
for n in range(Nv+Nd) :
    Sv=
    DATA[:,n]=array( [ float(Lfin[4+Nv+i].split()[n]) for i in range(Np) ] )