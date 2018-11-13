#! /usr/bin/env pytho

# ======================================================================
# matscipy - Python materials science tools
# https://github.com/libAtoms/matscipy
#
# Copyright (2014) James Kermode, King's College London
#                  Lars Pastewka, Karlsruhe Institute of Technology
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ======================================================================

from __future__ import print_function

import random
import unittest

import numpy as np
from numpy.linalg import norm

import ase.io as io
from ase.constraints import StrainFilter, UnitCellFilter
from ase.lattice.compounds import B1, B2, L1_0, L1_2
from ase.lattice.cubic import FaceCenteredCubic
from ase.lattice.hexagonal import HexagonalClosedPacked
from ase.optimize import FIRE
from ase.units import GPa

import matscipytest
from matscipy.calculators.eam import EAM
from matscipy.elasticity import fit_elastic_constants, Voigt_6x6_to_cubic

###

class TestEAMCalculator(matscipytest.MatSciPyTestCase):

    disp = 1e-6
    tol = 2e-6

    def test_forces(self):
        for calc in [EAM('Au-Grochola-JCP05.eam.alloy')]:
            a = io.read('Au_923.xyz')
            a.center(vacuum=10.0)
            a.set_calculator(calc)
            f = a.get_forces()
            fn = calc.calculate_numerical_forces(a)
            self.assertArrayAlmostEqual(f, fn, tol=self.tol)

    def test_stress(self):
        a = FaceCenteredCubic('Au', size=[2,2,2])
        calc = EAM('Au-Grochola-JCP05.eam.alloy')
        a.set_calculator(calc)
        self.assertArrayAlmostEqual(a.get_stress(), calc.calculate_numerical_stress(a), tol=self.tol)

        sx, sy, sz = a.cell.diagonal()
        a.set_cell([sx, 0.9*sy, 1.2*sz], scale_atoms=True)
        self.assertArrayAlmostEqual(a.get_stress(), calc.calculate_numerical_stress(a), tol=self.tol)

        a.set_cell([[sx, 0.1*sx, 0], [0, 0.9*sy, 0], [0, -0.1*sy, 1.2*sz]], scale_atoms=True)
        self.assertArrayAlmostEqual(a.get_stress(), calc.calculate_numerical_stress(a), tol=self.tol)

    def test_Grochola(self):
        a = FaceCenteredCubic('Au', size=[2,2,2])
        calc = EAM('Au-Grochola-JCP05.eam.alloy')
        a.set_calculator(calc)
        FIRE(StrainFilter(a, mask=[1,1,1,0,0,0]), logfile=None).run(fmax=0.001)
        a0 = a.cell.diagonal().mean()/2
        self.assertTrue(abs(a0-4.0701)<2e-5)
        self.assertTrue(abs(a.get_potential_energy()/len(a)+3.924)<0.0003)
        C, C_err = fit_elastic_constants(a, symmetry='cubic', verbose=False)
        C11, C12, C44 = Voigt_6x6_to_cubic(C)
        self.assertTrue(abs((C11-C12)/GPa-32.07)<0.7)
        self.assertTrue(abs(C44/GPa-45.94)<0.5)

    def test_CuAg(self):
        a = FaceCenteredCubic('Cu', size=[2,2,2])
        calc = EAM('CuAg.eam.alloy')
        a.set_calculator(calc)
        FIRE(StrainFilter(a, mask=[1,1,1,0,0,0]), logfile=None).run(fmax=0.001)
        e_Cu = a.get_potential_energy()/len(a)

        a = FaceCenteredCubic('Ag', size=[2,2,2])
        a.set_calculator(calc)
        FIRE(StrainFilter(a, mask=[1,1,1,0,0,0]), logfile=None).run(fmax=0.001)
        e_Ag = a.get_potential_energy()/len(a)
        self.assertTrue(abs(e_Ag+2.85)<1e-6)

        a = L1_2(['Ag', 'Cu'], size=[2,2,2], latticeconstant=4.0)
        a.set_calculator(calc)
        FIRE(UnitCellFilter(a, mask=[1,1,1,0,0,0]), logfile=None).run(fmax=0.001)
        e = a.get_potential_energy()
        syms = np.array(a.get_chemical_symbols())
        self.assertTrue(abs((e-(syms=='Cu').sum()*e_Cu-
                               (syms=='Ag').sum()*e_Ag)/len(a)-0.096)<0.0005)

        a = B1(['Ag', 'Cu'], size=[2,2,2], latticeconstant=4.0)
        a.set_calculator(calc)
        FIRE(UnitCellFilter(a, mask=[1,1,1,0,0,0]), logfile=None).run(fmax=0.001)
        e = a.get_potential_energy()
        syms = np.array(a.get_chemical_symbols())
        self.assertTrue(abs((e-(syms=='Cu').sum()*e_Cu-
                               (syms=='Ag').sum()*e_Ag)/len(a)-0.516)<0.0005)

        a = B2(['Ag', 'Cu'], size=[2,2,2], latticeconstant=4.0)
        a.set_calculator(calc)
        FIRE(UnitCellFilter(a, mask=[1,1,1,0,0,0]), logfile=None).run(fmax=0.001)
        e = a.get_potential_energy()
        syms = np.array(a.get_chemical_symbols())
        self.assertTrue(abs((e-(syms=='Cu').sum()*e_Cu-
                               (syms=='Ag').sum()*e_Ag)/len(a)-0.177)<0.0003)

        a = L1_2(['Cu', 'Ag'], size=[2,2,2], latticeconstant=4.0)
        a.set_calculator(calc)
        FIRE(UnitCellFilter(a, mask=[1,1,1,0,0,0]), logfile=None).run(fmax=0.001)
        e = a.get_potential_energy()
        syms = np.array(a.get_chemical_symbols())
        self.assertTrue(abs((e-(syms=='Cu').sum()*e_Cu-
                               (syms=='Ag').sum()*e_Ag)/len(a)-0.083)<0.0005)

    def test_CuZr(self):
        # This is a test for the potential published in:
        # Mendelev, Sordelet, Kramer, J. Appl. Phys. 102, 043501 (2007)
        a = FaceCenteredCubic('Cu', size=[2,2,2])
        calc = EAM('CuZr_mm.eam.fs', kind='eam/fs')
        a.set_calculator(calc)
        FIRE(StrainFilter(a, mask=[1,1,1,0,0,0]), logfile=None).run(fmax=0.001)
        a_Cu = a.cell.diagonal().mean()/2
        #print('a_Cu (3.639) = ', a_Cu)
        self.assertAlmostEqual(a_Cu, 3.639, 3)

        a = HexagonalClosedPacked('Zr', size=[2,2,2])
        a.set_calculator(calc)
        FIRE(StrainFilter(a, mask=[1,1,1,0,0,0]), logfile=None).run(fmax=0.001)
        a, b, c = a.cell/2
        #print('a_Zr (3.220) = ', norm(a), norm(b))
        #print('c_Zr (5.215) = ', norm(c))
        self.assertAlmostEqual(norm(a), 3.220, 3)
        self.assertAlmostEqual(norm(b), 3.220, 3)
        self.assertAlmostEqual(norm(c), 5.215, 3)

        # CuZr3
        a = L1_2(['Cu', 'Zr'], size=[2,2,2], latticeconstant=4.0)
        a.set_calculator(calc)
        FIRE(UnitCellFilter(a, mask=[1,1,1,0,0,0]), logfile=None).run(fmax=0.001)
        self.assertAlmostEqual(a.cell.diagonal().mean()/2, 4.324, 3)

        # Cu3Zr
        a = L1_2(['Zr', 'Cu'], size=[2,2,2], latticeconstant=4.0)
        a.set_calculator(calc)
        FIRE(UnitCellFilter(a, mask=[1,1,1,0,0,0]), logfile=None).run(fmax=0.001)
        self.assertAlmostEqual(a.cell.diagonal().mean()/2, 3.936, 3)

        # CuZr
        a = B2(['Zr', 'Cu'], size=[2,2,2], latticeconstant=3.3)
        a.set_calculator(calc)
        FIRE(UnitCellFilter(a, mask=[1,1,1,0,0,0]), logfile=None).run(fmax=0.001)
        self.assertAlmostEqual(a.cell.diagonal().mean()/2, 3.237, 3)

###

if __name__ == '__main__':
    unittest.main()
