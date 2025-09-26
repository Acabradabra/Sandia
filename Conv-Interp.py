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

# file_in='/mnt/scratch/PRECIZE/Sandia-Jaravel/RUN-D100-02-Laera/INIT-Fluent/Interp-00-Fluent.ip'
# file_in='/mnt/scratch/PRECIZE/Sandia-Jaravel/RUN-D100-02-Laera/DUMP-02-EDC-PB-OD2-FD3/DATA/Interp-EDC-FD3.ip'
file_ou=file_in[:-3]+'-New.ip'

#===================================================================================
util.Section('Reading : {:.3f}'.format(time.time()-t0),1,5,'r')
#===================================================================================
fin=open(file_in,'r')
Lfin=fin.readlines()
fin.close()
util.Section('Reading 1 : {:.3f}'.format(time.time()-t0),0,1,'b')

Nd=int(Lfin[1])
Np=int(Lfin[2])
Nv=int(Lfin[3])
Nl=len(Lfin)

Fields_in=[ f[:-1] for f in Lfin[4:4+Nv]]

print( '=> Nd={:d}, Np={:d}, Nv={:d}, Nl={:d}'.format(Nd,Np,Nv,Nl) )
print( '=> Fields_in : {}'.format(Fields_in) )

DATA_in=zeros((Np,Nv+Nd))
for n in range(Nv+Nd) :
    n0=4+Nv+n*(Np+1)
    Sv=Lfin[n0:n0+Np] ; Sv[0]=Sv[0][1:]
    DATA_in[:,n]=array( [float(s) for s in Sv ] )

#===================================================================================
util.Section('Processing : {:.3f}'.format(time.time()-t0),1,5,'r')
#===================================================================================

DATA_ou=DATA_in.copy()

#===================================================================================
util.Section('Writing : {:.3f}'.format(time.time()-t0),1,5,'r')
#===================================================================================

Fields_ou=Fields_in.copy()+['fmean']
Nv=len(Fields_ou)

fou=open(file_ou,'w')
fou.write( '3\n' )
fou.write( '{}\n'.format(Nd) )
fou.write( '{}\n'.format(Np) )
fou.write( '{}\n'.format(Nv) )
for f in Fields_ou :
    fou.write( '{}\n'.format(f) )
for v in range(Nv) :
    fou.write('(')
    for p in range(Np) : fou.write( ' {:.12e}\n'.format(DATA_ou[p,v]) )
    fou.write(')\n')
fou.close()

#===================================================================================
util.Section('Programme completed : {:.3f}'.format(time.time()-t0),1,5,'r')
#===================================================================================