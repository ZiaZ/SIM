#coding=utf-8
#!/usr/bin/env python

import sys

import numpy as np
from scipy.interpolate import interp1d
import ase.io
from ase import units
from ase.md.langevin import Langevin
from ase.md.nvtberendsen import NVTBerendsen
from ase.io.netcdftrajectory import NetCDFTrajectory
from ase.atoms import Atoms
from ase.md import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as mbd
from ase.optimize.fire import FIRE

# from matscipy.fracture_mechanics.idealbrittlesolid import (IdealBrittleSolid,
#                                                            triangular_lattice_slab,
#                                                            find_crack_tip,
#                                                            set_initial_velocities,
#                                                            set_constraints,
#                                                            extend_strip)

from matscipy_local.matscipy.fracture_mechanics.idealbrittlesolid import (IdealBrittleSolid,
                                                           triangular_lattice_slab,
                                                           find_crack_tip,
                                                           set_initial_velocities,
                                                           set_constraints,
                                                           extend_strip)


from matscipy.fracture_mechanics.crack import (thin_strip_displacement_y,
                                               ConstantStrainRate)

sys.path.insert(0, '.')

import modifiers
import create_folder as cf

def ribs(params, frame_count = 1000):
       
       calc = IdealBrittleSolid(rc=params.rc, k=params.k, a=params.a, beta=params.beta)

       x_dimer = np.linspace(params.a-(params.rc-params.a),
                            params.a+1.1*(params.rc-params.a),51)
       dimers = [Atoms('Si2', [(0, 0, 0), (x, 0, 0)],
                     cell=[10., 10., 10.], pbc=True) for x in x_dimer]
       calc.set_reference_crystal(dimers[0])
       e_dimer = []
       f_dimer = []
       f_num = []
       for d in dimers:
              d.set_calculator(calc)
              e_dimer.append(d.get_potential_energy())
              f_dimer.append(d.get_forces())
              f_num.append(calc.calculate_numerical_forces(d))
       e_dimer = np.array(e_dimer)
       f_dimer = np.array(f_dimer)
       f_num = np.array(f_num)
       assert abs(f_dimer - f_num).max() < 0.1

       #! crystal is created here, the length and height can be modified here as well
       #! edit 3N changed to different values to test
       crystal = triangular_lattice_slab(params.a, params.lm*params.N, params.N)
       calc.set_reference_crystal(crystal)
       crystal.set_calculator(calc)

       e0 = crystal.get_potential_energy()
       l = crystal.cell[0,0]
       h = crystal.cell[1,1]
       print('l=', l, 'h=', h)

       # compute surface (Griffith) energy
       b = crystal.copy()
       b.set_calculator(calc)
       shift = calc.parameters['rc']*2
       y = crystal.positions[:, 1]
       b.positions[y > h/2, 1] += shift
       b.cell[1, 1] += shift
       e1 = b.get_potential_energy()
       E_G = (e1 - e0)/l
       print ('Griffith energy', E_G)

       # compute Griffith strain
       eps = 0.0   # initial strain is zero
       eps_max = 2/np.sqrt(3)*(params.rc-params.a)*np.sqrt(params.N-1)/h # Griffith strain assuming harmonic energy
       deps = eps_max/100. # strain increment
       e_over_l = 0.0     # initial energy per unit length is zero
       energy = []
       strain = []
       while e_over_l < E_G:
              c = crystal.copy()
              c.set_calculator(calc)
              c.positions[:, 1] *= (1.0 + eps)
              c.cell[1,1] *= (1.0 + eps)
              e_over_l = c.get_potential_energy()/l
              energy.append(e_over_l)
              strain.append(eps)
              eps += deps

       energy = np.array(energy)
       eps_of_e = interp1d(energy, strain, kind='linear')
       eps_G = eps_of_e(E_G)

       print ('Griffith strain', eps_G)

       c = crystal.copy()
       c.info['E_G'] = E_G
       c.info['eps_G'] = eps_G

       # open up the cell along x and y by introducing some vaccum
       orig_cell_width = c.cell[0, 0]
       orig_cell_height = c.cell[1, 1]
       c.center(params.vacuum, axis=0)
       c.center(params.vacuum, axis=1)

       # centre the slab on the origin
       c.positions[:, 0] -= c.positions[:, 0].mean()
       c.positions[:, 1] -= c.positions[:, 1].mean()

       c.info['cell_origin'] = [-c.cell[0,0]/2, -c.cell[1,1]/2, 0.0]
       ase.io.write('crack_1.xyz', c, format='extxyz')

       width = (c.positions[:, 0].max() -
              c.positions[:, 0].min())
       height = (c.positions[:, 1].max() -
              c.positions[:, 1].min())

       c.info['OrigHeight'] = height

       print(('Made slab with %d atoms, original width and height: %.1f x %.1f A^2' %
              (len(c), width, height)))

       top = c.positions[:, 1].max()
       bottom = c.positions[:, 1].min()
       left = c.positions[:, 0].min()
       right = c.positions[:, 0].max()

       crack_seed_length = 0.2*width
       strain_ramp_length = 8.0*params.a # make this bigger until crack looks nicer
       delta_strain = params.strain_rate*params.dt

       # fix top and bottom rows, and setup Stokes damping mask
       # initial use constant strain
       set_constraints(c, params.a)

       # apply initial displacment field
       c.positions[:, 1] += thin_strip_displacement_y(
                                   c.positions[:, 0],
                                   c.positions[:, 1],
                                   params.delta*eps_G,
                                   left + crack_seed_length,
                                   left + crack_seed_length +
                                          strain_ramp_length)

       print('Applied initial load: delta=%.2f strain=%.4f' %
       (params.delta, params.delta*eps_G))

       ase.io.write('crack_2.xyz', c, format='extxyz')

       c.set_calculator(calc)

       cl, cs, cr = calc.get_wave_speeds(c)

       print("rayleigh speed = %f" %cr)

       # relax initial structure
       # opt = FIRE(c)
       # opt.run(fmax=1e-3)

       ase.io.write('crack_3.xyz', c, format='extxyz')

       #Atomic Void Simulations ""
       void = "😵"
       #------------------------------------------------------------------------------
       #the following lines of code were written to convert 1D positions into a 2D array
       #so as to make manipulation of slab easier
       if True:
              L = params.N*params.lm
              H = params.N*2
              #this reference code is for 160x40 slab
              #in steps of 40 i.e. the height, create a list upto the length of slab
              #row0 = range(0,6400,40)
              row0 = range(0,len(c),H)
              #the 2D array will be held in slab
              slab = [] 
              for col in range(H):
                     row = []
                     for r in row0:
                            i = col + r
                            row.append(c.positions[i])
                     slab.append(row)

              slab = np.array(slab)
              #all items in the array reversed, needed because the salb is built from bottom left up
              #in other words, reflected in x axis
              slab = slab[::-1]
              # print(slab[0])
              
              #around a max of [0.3--1.7] recommended
              mid_offset = 1.0
              #y offset is a fraction of distance from the end of the slab
              y_offset = 0.5
              rad = 3
              slab[modifiers.mask(h=H,w=L, center=[int(L*y_offset),int((H/2 - 1)*mid_offset)], radius=rad)] = 0
              #reversed the slab back again here
              slab = slab[::-1]
              # # this is a useful text-array representation of the slab, for debugging purposes
              # mtext = open('masktest.txt','w')
              # for row in slab:
              #        mtext.write(str(row).replace("\n",",")+"\n")
              # mtext.close
              
              slab_1d = []
              for col in range(L):
                     for row in range(H):
                            slab_1d.append(slab[row,col])
              slab_1d = np.array(slab_1d)
              
              todel = []
              for i in range(len(c)):
                     if slab_1d[i][2] == 0:
                            todel.append(i)
              print(todel)
              del c[todel]
              # return
       #-------------------------------------------------------------------------------
       #! replaced velcityVerlet with Lagevin to add temperature parameter
       if params.v_verlet:
              dyn = VelocityVerlet(c, params.dt * units.fs, logfile=None)
              # set_initial_velocities(dyn.atoms)
       else: 
              print("Using NVT!")
              # mbd(c, 20 * units.kB, force_temp = True)
              # dyn = Langevin(c,params.dt*units.fs,params.T*units.kB, 5)
              dyn = NVTBerendsen(c, params.dt * units.fs, params.T, taut=0.5*1000*units.fs)

       #dyn.atoms.rattle(1e-3) # non-deterministic simulations - adjust to suit

       #!simulation outputs numbered, avoids deleting exisiting results
       iterFile = open("simIteration.txt",'r+')
       iteration = int(iterFile.readlines()[0]) + 1
       if params.overwrite_output:
              iteration -= 1
       iterFile.seek(0,0)
       iterFile.write(str(iteration))
       iterFile.close()

       dir = "./.simout/sim_"+str(iteration)+"_"+str(params.keep_test).lower()+"_"+params.desc+"/"
       cf.createFolder(dir)

       #!Saving parameter values for each iteration of the simulation
       logFile = open(dir+"params_log.txt",'w')
       for p, value in params.compose_params().iteritems():
              logFile.write(p+" ==> "+str(value)+"\n")
       logFile.close
       crack_pos = []
       if params.keep_test:
              traj = NetCDFTrajectory(dir+'traj'+str(iteration)+'.nc', 'w', c)
              dyn.attach(traj.write, 10, dyn.atoms, arrays=['stokes', 'momenta'])

       # #! isolating crack tip_x for saving
       # crack_tip_file2 = open(dir+'tip_x.txt','w')
       # crack_tip_file2.close()
       tip_x_file = open(dir+'tip_x.txt','a')
       console_output = open(dir+'console_output.txt','a')
       coord_file = open(dir+'coordinates.csv','a')

       coordinates = []
       distances = []

       dyn.attach(find_crack_tip, 10, dyn.atoms, tipxfile = tip_x_file, cout=console_output, coord=coordinates, d=distances,
              dt=params.dt*10, store=True, results=crack_pos)

       # run for 2000 time steps to reach steady state at initial load
       # for i in range(10):
       #     dyn.run(250)
       #     if extend_strip(dyn.atoms, params.a, params.N, params.M, params.vacuum):
       #         set_constraints(dyn.atoms, params.a)

       # start decreasing strain
       #set_constraints(dyn.atoms, params.a, delta_strain=delta_strain)

       # strain_atoms = ConstantStrainRate(dyn.atoms.info['OrigHeight'],
       #                                   delta_strain)
       # dyn.attach(strain_atoms.apply_strain, 1, dyn. atoms   )

       # for i in range(50):
       #     dyn.run(100)
       #     if extend_strip(dyn.atoms, params.a, params.N, params.M, params.vacuum):
       #         set_constraints(dyn.atoms, params.a)

       # #cleardel dyn.observers[-1] # stop increasing the strain

       # for i in range(1000):
       #     dyn.run(100)
       #     if extend_strip(dyn.atoms, params.a, params.N, params.M, params.vacuum):
       #         set_constraints(dyn.atoms, params.a)
       
       
       dyn.run(int(1*frame_count)*10+10)

       # print("\n\n\n\n\n -----Adding Temperature------\n\n\n\n\n")
       # # mbd(c, 2*params.T * units.kB, force_temp = True)
       # dyn.set_temperature(params.T*units.kB)
       # dyn.run(int(0.5*frame_count)*10+10)
       for c in coordinates:
              coord_file.write(str(c[0])+','+str(c[1])+'\n')
       coord_file.close()

       if params.keep_test:
              traj.close()
       tip_x_file.close()
       console_output.close()
       # time = 10.0*dyn.dt*np.arange(dyn.get_number_of_steps()/10)
       # np.savetxt(dir+'crackpos.txt', np.c_[time, crack_pos])

if __name__ == '__main__':
       import params
       params.keep_test = True
       params.delta = 2
       params.k = 0.5
       params.v_verlet = True
       params.desc = "hole_test"
       ribs(params, frame_count = 2000)