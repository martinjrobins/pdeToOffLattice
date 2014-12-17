# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 00:18:14 2014

@author: mrobins
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Dec 15 14:05:32 2014

@author: robinsonm
"""

from fipy import *
from scipy.special import erf
from math import sqrt
import pyTyche as tyche


nx = 50
D = 1.
L = 1.
dx = L/nx
conversion_rate = 100.

interface = L/2.0

valueLeft = 1
valueRight = 0

timeStepDuration = 9 * dx**2 / (2 * D)
steps = 200
mol_dt = timeStepDuration/100.

A = tyche.new_species(D)
A.set_grid(tyche.new_structured_grid([0,0,0],[L,0,0],[dx,0,0]))

xminboundary = tyche.new_reflective_boundary(tyche.new_xplane(0,1))
xmaxboundary = tyche.new_reflective_boundary(tyche.new_xplane(L,-1))

#set off-lattice sink and source
sink = tyche.new_uni_reaction(conversion_rate,[[A],[A.pde()]])
sink.set_geometry(0,tyche.new_xplane(interface,-1))
source = tyche.new_zero_reaction_lattice(conversion_rate,[[A.pde()],[A]])
source.set_geometry(tyche.new_xplane(interface,1))

algorithm = tyche.group([tyche.new_diffusion(), sink, xminboundary, xmaxboundary])
algorithm.add_species(A)

mesh = Grid1D(nx=nx, dx=dx)
phi = CellVariable(name="solution variable",  mesh=mesh, value=1.)

print 'begin loop'
for step in range(steps):
    #set off-lattice generators 
    #(set B compartments equal to pde)
    A.set_pde(phi.value)
    print 'middle loop'

    #integrate off-lattice model    
    time = algorithm.integrate_for_time(timeStepDuration,mol_dt)
    
    

