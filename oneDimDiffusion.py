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


mesh = Grid1D(nx=nx, dx=dx)

phi = CellVariable(name="solution variable",  mesh=mesh, value=0.)
phi.constrain(valueRight, mesh.facesRight)
phi.constrain(valueLeft, mesh.facesLeft)

eq = TransientTerm() == DiffusionTerm(coeff=D) - conversion_rate.*ImplicitSourceTerm((mesh.x > interface))

A = tyche.new_species(D,new_structured_grid([0,0,0],[L,0,0],[dx,0,0])

xminboundary = tyche.new_reflective_boundary(tyche.new_xplane(0,1),[L,0,0])
xmaxboundary = tyche.new_reflective_boundary(tyche.new_xplane(L,-1),[-L,0,0])

#set off-lattice sink and source
sink = tyche.new_uni_reaction(conversion_rate,[[A],[A.lattice]],tyche.new_xplane(interface,-1))
source = tyche.new_uni_reaction(conversion_rate,[[A.lattice],[A]],tyche.new_xplane(interface,1))

algorithm = tyche.group([tyche.new_diffusion(), xminboundary, xmaxboundary, sink, source])
algorithm.add_species(A)


for step in range(steps):
    #set off-lattice generators 
    #(set B compartments equal to pde)
    A.set_compartments(np.array(phi.value))

    #integrate pde model    
    eq.solve(var=phi, dt=timeStepDuration)    
    
    #integrate off-lattice model    
    time = algorithm.integrate_for_time(timeStepDuration,mol_dt)
    
    #transfer sink particles to pde  
    #(transfer change in B to phi)
    phi.value += A.get_compartments - np.array(phi.oldValue)
    
    

phiAnalytical = CellVariable(name="analytical value", mesh=mesh)
x = mesh.cellCenters[0]
t = timeStepDuration * steps
phiAnalytical.setValue(1 - erf(x / (2 * sqrt(D * t))))
viewer = Viewer(vars=(phi, phiAnalytical), datamin=0., datamax=1.)
viewer.plot()
raw_input("Press <return> to proceed...")

