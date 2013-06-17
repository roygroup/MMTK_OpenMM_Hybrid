#This script simulates 2 spcfw-q water molecules using OpenMM's GPU libraries
#The 2 spcfw-q's are held together by a centre of mass restraint
#Last modified: June 17, 2013

"""
---------------------------------------------------------------------------------------
Copyright (c) 2013 Kevin Bishop
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
---------------------------------------------------------------------------------------
"""

from MMTK import *
from MMTK.ForceFields import LennardJonesForceField, Amber94ForceField
from MMTK.Trajectory import Trajectory, TrajectoryOutput, LogOutput, StandardLogOutput
from MMTK.Dynamics import VelocityVerletIntegrator, VelocityScaler, \
                          TranslationRemover, BarostatReset                   
from LangevinDynamics import LangevinIntegrator
from time import time
from math import sqrt
import sys,os,string

# Define thermodynamic parameters from input
temperature = 5.*Units.K
dt = 2.0*Units.fs
EquilSteps = 1000
ProdSteps = 10000
skipSteps = 100
nb = int(sys.argv[1])
k_restraint = float(sys.argv[2])
r0_restraint = float(sys.argv[3])
temperature = float(sys.argv[4])
dt = float(sys.argv[5])*Units.fs
#The friction and steps need to be chosen more systematically
friction=0.001/dt
EquilSteps = int(sys.argv[6])
ProdSteps = int(sys.argv[7])
skipSteps = int(sys.argv[8])

print "\n--------Parameters-------------------------------"
print "Number of beads:\t" + str(nb)
print "Restraint k:\t\t" + str(k_restraint) + " kJ/mol/nm^2"
print "Restraint r0:\t\t" + str(r0_restraint) + " nm"
print "Temperature:\t\t" + str(temperature) + " K"
print "Timestep:\t\t" + str(dt) + " ps"
print "Friction:\t\t" + str(friction)
print "Equilibration Steps:\t" + str(EquilSteps)
print "Production Steps:\t" + str(ProdSteps)
print "Skip Steps:\t\t" + str(skipSteps)
print "-----------End-----------------------------------\n"

#Create Universe
universe = InfiniteUniverse(Amber94ForceField(15.*Units.Ang,None))
universe.addObject(Environment.PathIntegrals(temperature, True))

#Add 2 spcfw-q's to the universe with a separation of 3 Angstroms
ar_dist = 3.0*Units.Ang
pos1 = Vector(0.0,0.0,0.0)
pos2 = Vector(ar_dist,0.0,0.0)
universe.addObject(Molecule('spcfw-q', position=pos1))
universe.addObject(Molecule('spcfw-q', position=pos2))

#Set the number of beads for each atom in the universe
for atom in universe.atomList():
        atom.setNumberOfBeads(nb)
universe.environmentObjectList(Environment.PathIntegrals)[0].include_spring_terms = False
universe._changed(True)

# Initialize velocities to temperature
universe.initializeVelocitiesToTemperature(temperature)

#This creates the integrator object by creating the class defined in LangevinDynamics.py
#Note that the restraint parameter is a list where the first element is a list of lists. The first list
#contains the particle indices for which the first centre of mass should be calculated and the second list
#contains the particle indices for the second centre of mass, the 2nd parameter is the k value for the restraint,
#third is the equilibrium value, and finally a string that is 'com' for centre of mass or 'cog' for centre of geometry (and more could be added)
integrator = LangevinIntegrator(universe, delta_t=dt,
                                friction=friction, temperature=temperature, restraint=[[0,1,2],[3,4,5],k_restraint,r0_restraint,'com'])

#Create the trajectory files for both equilbration and production runs
dir = "spcfw-q_2_"+str(nb)+"_"+str(k_restraint)+"_"+str(r0_restraint)
try:
  os.mkdir(dir)
except (OSError):
  print 'Do not need to create directory, '  + dir +  ', it is already present'                     

trajectory_eq = Trajectory(universe,   dir + "//" + dir + "_eq.nc", "w")
trajectory_prod = Trajectory(universe, dir + "//" + dir + "_prod.nc", "w")

# Periodical actions for trajectory output and text log output.
eq_output_actions = [TrajectoryOutput(trajectory_eq,
                                   ('configuration', 'energy', 'thermodynamic',
                                    'time', 'auxiliary','velocities'), 0, None, 100)]
#                                    StandardLogOutput(100)]
prod_output_actions = [TrajectoryOutput(trajectory_prod,
                                   ('configuration', 'energy', 'thermodynamic',
                                    'time', 'auxiliary','velocities'), 0, None, skipSteps)]

#Perform the equilibration portion
print '.....Beginning Equlibration.....................................................'

integrator(steps = EquilSteps, actions =  eq_output_actions)
trajectory_eq.close()

print '.....Equilibrating Done.........................................................'

#Perform the production portion
start = time()       
print '\n.....Beginning Production.......................................................'     
integrator = LangevinIntegrator(universe, delta_t=dt,
                                friction=friction, temperature=temperature, restraint=[[0,1,2],[3,4,5],k_restraint,r0_restraint,'com'])      
integrator(steps = ProdSteps, actions = prod_output_actions)

print '.....Production Done............................................................'
end = time()

print '\n.....Production Statistics......................................................'
print 'Completed', ProdSteps, 'production steps in ', (end-start), 'seconds or ', (end-start)/ProdSteps, 'seconds/step'
num_ns_prod = int(ProdSteps)*float(dt) / 1000.0
ns_per_day = num_ns_prod * 86400.0 / (end-start)
print 'Production was done at a rate of ' + str(ns_per_day) + ' ns/day'
trajectory_prod.close()
print '.....End of Production Statistics...............................................\n'

####################################################################################################################################################################
#This is a function that returns the centre of mass between 2 lists of atoms for a particular bead number
#It is used in the analysis portion of this script only.
def BeadToBeadCOM(atoms1, atoms2, bead):
  x = 0       #These three values are just the x,y,z indices of a particle
  y = 1
  z = 2

  m1 = 0.       #Sum of the masses for list of atoms 1
  m2 = 0.       #Sum of the masses for list of atoms 2

  x1 = 0.
  y1 = 0.
  z1 = 0.
  x2 = 0.
  y2 = 0.
  z2 = 0.

  for a1 in atoms1:
    m1 += a1.mass()

    x1 += a1.beadPositions()[bead][x] * a1.mass()

    y1 += a1.beadPositions()[bead][y] * a1.mass()

    z1 += a1.beadPositions()[bead][z] * a1.mass()

  for a2 in atoms2:
    m2 += a2.mass()

    x2 += a2.beadPositions()[bead][x] * a2.mass()

    y2 += a2.beadPositions()[bead][y] * a2.mass()

    z2 += a2.beadPositions()[bead][z] * a2.mass()


  return sqrt(((x1/m1) - (x2/m2))**2 + ((y1/m1) - (y2/m2))**2 + ((z1/m1) - (z2/m2))**2)

####################################################################################################################################################################

print '.....Performing Analysis........................................................'
#Read in trajectory from file
trajectory = Trajectory(None, dir + '//' + dir + '_prod.nc')
#Output the variables stored in the trajectory
print 'The available variables in the trajectory are:', trajectory.variables()
universe = trajectory.universe
natoms = universe.numberOfAtoms()
P = universe[0].atomList()[0].numberOfPoints()
print 'There are', natoms, 'atoms and', P, 'beads per atom present in the trajectory'

#Output files to write centre of masses to.
com_out = open(dir + '//' + dir + '_centroid_com_distance.txt','w')
bead_com_out = open(dir + '//' + dir + '_bead_com_distance.txt','w')

for step in trajectory:
    #set the configuration to each new step
    universe.setConfiguration(step['configuration'])

    for object1 in universe.objectList():
        for object2 in universe.objectList():
            if object1 < object2:
                #plain centroid centre of mass
                com_out.write(str(universe.distance(object1.centerOfMass(),object2.centerOfMass())) + '\n')
                atoms1 = object1.atomList()
                atoms2 = object2.atomList()
                #Bead to bead centre of mass
                for bead in range(P):
                  bead_com_out.write(str(BeadToBeadCOM(atoms1,atoms2,bead)) + "\n")
com_out.close()
bead_com_out.close()

print ".....Analysis Complete.........................................................."

print '\nScript Complete!\n'