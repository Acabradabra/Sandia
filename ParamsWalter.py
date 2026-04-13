#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#====================> Directories
dir0='/mnt/d/Python/Sandia/'
dirs=dir0+'DATA-Oxy/Mean_RMS/'
dirc='/mnt/scratch/ZEUS/FLUENT/PRECIZE/Walter-50pH2-60pJet/DUMP/'
# dirc='/mnt/scratch/ZEUS/FLUENT/PRECIZE/Walter-50pH2-60pJet/DUMP-06-NFR/'
dird=dirc+'DATA/'
dirp=dirc+'PLOT/'
# slice='Data-all.dat'
slice_m='Data-FS.dat'
slice_f='Data-F.dat'
slice_o='Data-OUT.dat'

BC_f={'unit':'Y'}
BC_o={'unit':'Y'}

ray_p=2e-2
Pos_p=[
    [3.098,1.712],
    [2.17 ,1.712],
    [0.51 ,1.712],
    [2.296,0.926],
    [1.618,0.869],
    [1.417,1.435],
    [1.1575,2.19]
]
Txt_p=[ '1','2','3' , 'a','c','d','e' ]