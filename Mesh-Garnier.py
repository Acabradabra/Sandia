#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Utilities as util
(Sysa,NSysa,Arg)=util.Parseur(['BD','Conv','Mix','VelProf','Pos','Flow'],0,'Arg : ')
(                             [ BD , CONV , MIX , VELPROF , POS , FLOW ])=Arg

from numpy import *
import sys
import time
import pygmsh as pg
import gmsh   as gm
import Meshing as me
import meshio as mo

import matplotlib.pyplot as plt
import matplotlib

from Params import *

t0=time.time()
#===================================================================================
#                     Parameters
#===================================================================================

