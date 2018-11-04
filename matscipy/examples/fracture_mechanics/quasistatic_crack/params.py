#! /usr/bin/env python

from math import log

import numpy as np

import ase.io

from matscipy.fracture_mechanics.clusters import diamond, set_groups

import atomistica

###

# Interaction potential
calc = atomistica.KumagaiScr()

# Fundamental material properties
el              = 'Si'
a0              = 5.429
C11             = 166.   # GPa
C12             = 65.    # GPa
C44             = 77.    # GPa
surface_energy  = 1.08  * 10    # GPa*A = 0.1 J/m^2

# Crack system
n               = [ 4, 6, 1 ]
crack_surface   = [ 1,-1, 0 ]
crack_front     = [ 1, 1, 0 ]
crack_tip       = [ 41, 56 ]
skin_x, skin_y = 1, 1

vacuum          = 6.0

# Simulation control
nsteps          = 31
# Increase stress intensity factor
k1              = np.linspace(0.8, 1.5, nsteps)
# Don't move crack tip
tip_dx          = np.zeros_like(k1)
tip_dz          = np.zeros_like(k1)

fmax            = 0.05

# Setup crack system
cryst = diamond(el, a0, n, crack_surface, crack_front)
set_groups(cryst, n, skin_x, skin_y)

ase.io.write('cryst.cfg', cryst)

# Compute crack tip position
r0 = np.sum(cryst.get_positions()[crack_tip,:], axis=0)/len(crack_tip)
tip_x0, tip_y0, tip_z0 = r0
