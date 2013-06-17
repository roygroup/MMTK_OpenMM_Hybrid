# This module implements a Langevin integrator.
# Last modified: June 17, 2013

"""
Copyright (c) 2013 Stephen Constable, Nabil Faruk, Kevin Bishop, Pierre-Nicholas Roy
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
"""

from MMTK import Dynamics, Environment, Features, Trajectory, \
                 Units, ParticleProperties, Universe, Molecule
import numpy as N
from MMTK.Features import PathIntegralsFeature
import MMTK_langevin
import MMTK_forcefield
import operator, copy
from RestraintString import CreateRestraintString

#
# Langevin integrator
#
class LangevinIntegrator(Dynamics.Integrator):

    def __init__(self, universe, **options):
        Dynamics.Integrator.__init__(self, universe, options)
        # Supported features: path integrals, to do PIMD
        self.features = [PathIntegralsFeature]
        assert PathIntegralsFeature.isInUniverse(universe)
        self.nbeads = universe.maxNumberOfBeads()
	
    def initOmmSystem(self):
    	atoms = copy.copy(self.universe.atomList())
    	atoms.sort(key=operator.attrgetter('index'))
    	masses = [atom.mass() for atom in atoms]
    	MMTK_langevin.initOpenMM(N.array(masses, dtype=N.float64), len(atoms))
		
    def destroyOmmSystem(self):
		MMTK_langevin.destroyOpenMM()
	
    def addOmmCMRemover(self, skip):
        MMTK_langevin.addCMRemover(skip)
    
    def createOmmForces(self):
        params = self.universe.energyEvaluatorParameters()
        natoms = len(self.universe.atomList())
        nbeads = self.nbeads
        
        #here we convert MMTK's internal representation of the forcefield parameters
        #to OpenMM.  For a new MMTK forcefield to work, it must implement the 
        #energyEvaluatorParameters method, and this code must be modified to make use
        #of those parameters
        
        #first, lets start with the bonded forcefield
        bonds = params['harmonic_distance_term']

        atom_index_1 = []
        atom_index_2 = []
        eq_distance = []
        spring_k = []
        
        for bond in bonds:
        	if ((bond[0] % nbeads == 0) and (bond[1] % nbeads == 0)):
        		atom_index_1.append(bond[0]/nbeads)
        		atom_index_2.append(bond[1]/nbeads)
        		eq_distance.append(bond[2])
        		spring_k.append(2.*bond[3]) # OpenMM divides by 2 so must counteract here by multiplying by 2

        MMTK_langevin.makeOmmBondedForce(N.array(atom_index_1, dtype=N.int32), 
        								 N.array(atom_index_2, dtype=N.int32), 
        								 N.array(eq_distance, dtype=N.float64),
        								 N.array(spring_k, dtype=N.float64),
        								 len(atom_index_1))
        
        
        #angles
        angles = params['harmonic_angle_term']
        
        atom_index_1 = []
        atom_index_2 = []
        atom_index_3 = []
        eq_angle = []
        spring_k = []
        
        for angle in angles:
          if ((angle[0] % nbeads == 0) and (angle[1] % nbeads == 0) and (angle[2] % nbeads == 0)):
            atom_index_1.append(angle[0]/nbeads)
            atom_index_2.append(angle[1]/nbeads)
            atom_index_3.append(angle[2]/nbeads)
            eq_angle.append(angle[3])
            spring_k.append(2.*angle[4]) # OpenMM divides by 2 so must counteract here by multiplying by 2
        
        MMTK_langevin.makeOmmAngleForce(N.array(atom_index_1, dtype=N.int32),
        								N.array(atom_index_2, dtype=N.int32),
        								N.array(atom_index_3, dtype=N.int32),
        								N.array(eq_angle, dtype=N.float64),
        								N.array(spring_k, dtype=N.float64),
        								len(atom_index_1))
        
        
        #dihedrals
        dihedrals = params['cosine_dihedral_term']
        
        atom_index_1 = []
        atom_index_2 = []
        atom_index_3 = []
        atom_index_4 = []
        periodicity = []
        eq_dihedral = []
        spring_k = []
        
        for dihedral in dihedrals:
        	if ((dihedral[0] % nbeads == 0) and (dihedral[1] % nbeads == 0) and (dihedral[2] % nbeads == 0) and (dihedral[3] % nbeads == 0)):
        		atom_index_1.append(dihedral[0]/nbeads)
        		atom_index_2.append(dihedral[1]/nbeads)
        		atom_index_3.append(dihedral[2]/nbeads)
        		atom_index_4.append(dihedral[3]/nbeads)
        		periodicity.append(dihedral[4])
        		eq_dihedral.append(dihedral[5])
        		spring_k.append(dihedral[6])
        
        MMTK_langevin.makeOmmDihedralForce(N.array(atom_index_1, dtype=N.int32),
        								   N.array(atom_index_2, dtype=N.int32),
        								   N.array(atom_index_3, dtype=N.int32),
        								   N.array(atom_index_4, dtype=N.int32),
        								   N.array(periodicity, dtype=N.int32),
        								   N.array(eq_dihedral, dtype=N.float64),
        								   N.array(spring_k, dtype=N.float64),
        								   len(atom_index_1))
        
        #now, nonbonded forces
        lj = params['lennard_jones']
        lj_14_factor = lj['one_four_factor']
        lj_cutoff = lj['cutoff']
        
        e_s_matrix = lj['epsilon_sigma']
        epsilon_list = []
        sigma_list = []
        
        for i in range(0, len(e_s_matrix)):
        	epsilon, sigma = e_s_matrix[i][i]
        	epsilon_list.append(epsilon)
        	sigma_list.append(sigma)
        
        #contains per-particle information, but we want per-atom
        e_s_types = lj['type']
        epsilon_p = [] # per-particle epsilons
        sigma_p = [] # per-particle sigmas
        
        for type in e_s_types:
        	epsilon_p.append(epsilon_list[type])
        	sigma_p.append(sigma_list[type])
        
        epsilon = [] # per-atom epsilons
        sigma = [] # per-atom sigmas
        
        for i in range(0,natoms*self.nbeads,self.nbeads):
        	epsilon.append(epsilon_p[i])
        	sigma.append(sigma_p[i])
        
        #now epsilon and sigma contain the correct values in the correct order
        
        elec = params['electrostatic']
        elec_14_factor = elec['one_four_factor']

        #check if periodic universe so that Ewald and proper cutoff is set
        is_periodic = (1,0)[isinstance(self.universe, Universe.InfiniteUniverse)]
        #print "periodicity: ", is_periodic

        if is_periodic:
            elec_cutoff = elec['real_cutoff'] #for step_0_
        else:
            elec_cutoff = elec['cutoff'] #for example
        
        charge_p = elec['charge'] # per-particle charges
        elec_method = elec['algorithm']
        
        charge = [] # per-atom charges
        for i in range(0,natoms*self.nbeads,self.nbeads):
        	charge.append(charge_p[i])
                
        MMTK_langevin.makeOmmEsAndLjForce(is_periodic,
                                          N.array(sigma, dtype=N.float64),
        								  N.array(epsilon, dtype=N.float64),
        								  N.array(charge, dtype=N.float64),
        								  elec_14_factor,
        								  lj_14_factor,
        								  elec_cutoff,
        								  len(sigma))

    #This will add a restraint to the forcefield    
    def addRestraint(self,restraint):
        
        g1 = restraint[0]
        g2 = restraint[1]
        k_compound = N.float64(restraint[2])
        r0_compound = N.float64(restraint[3])
        option = restraint[4]

        energy_exp = CreateRestraintString(g1,g2,k_compound,r0_compound,option)
        print energy_exp

        num_of_particles = len(g1) + len(g2)

        MMTK_langevin.addRestraint(N.array(restraint[0]+restraint[1], dtype=N.int32),
                                   energy_exp,num_of_particles)   

    def __call__(self, **options):
        # Process the keyword arguments
        self.setCallOptions(options)
        # Check if the universe has features not supported by the integrator
        Features.checkFeatures(self, self.universe)
        
        #the following two lines are required due to MMTK mixing up atom indices.  This sorts them.
        #atoms = self.universe.atomList()
        #atoms.sort(key=operator.attrgetter('index'))
        
        # Get the universe variables needed by the integrator
        configuration = self.universe.configuration()
        velocities = self.universe.velocities()
        if velocities is None:
            raise ValueError("no velocities")

        # Get the friction coefficients. First check if a keyword argument
        # 'friction' was given to the integrator. Its value can be a
        # ParticleScalar or a plain number (used for all atoms). If no
        # such argument is given, collect the values of the attribute
        # 'friction' from all atoms (default is zero).
        friction = self.getOption('friction')  
        
        #call this method to create the OpenMM System.  Without this, creating
        #the forces will fail!!
        self.initOmmSystem()
        
        
        #call this method to ensure that the forcefield parameters get
        #passed to OpenMM.  Without this, there will be no forces to integrate!
        self.createOmmForces()

        #The restraint option must be included in the creation of the integrator for this part to be executed
        try:
            restraint = self.getOption('restraint')
        except ValueError:
            restraint = None
            print 'No restraints (Value Error)'
        except:
            print sys.exc_info()[0]
        if restraint != None:
            self.addRestraint(restraint) 
        
        #ugly hack to get timestep skip information
        skip = -1
        if ('actions' in self.call_options):
            actions = self.call_options['actions']
            for action in actions:
                if isinstance(action, Dynamics.TranslationRemover):
                    self.addOmmCMRemover(action.skip)
                if isinstance(action, Trajectory.LogOutput):
                    break # to avoid Logoutput skip setting problem PNR NFF SJC
                if isinstance(action, Trajectory.TrajectoryOutput):
                    skip = action.skip
                    action.skip = 1
        
            if skip < 0 :
                for action in actions:
                    if issubclass(action.__class__, Trajectory.TrajectoryOutput):
                        skip = action.skip
                        action.skip = 1
                        break
        
        #we are not logging anything, so we don't have to report any intermediate values
        if skip < 0:
            skip = 100
            steps = 1

        # Run the C integrator
        atoms = self.universe.atomList()
        masses = [atom.mass() for atom in atoms]
        natoms = len(atoms)
        nmolecules = len(self.universe.objectList(Molecule))
        molecule_lengths = N.zeros(nmolecules, dtype=N.int32, order='C')
        for i in range(0,len(self.universe.objectList(Molecule))):
            molecule_lengths[i] = len(self.universe.objectList(Molecule)[i].atomList())

        MMTK_langevin.integrateLD(self.universe, configuration.array, 
                                  velocities.array, N.array(masses, dtype=N.float64), 
                                  friction,
                                  self.getOption('temperature'),
                                  self.getOption('delta_t'),
                                  self.getOption('steps'), skip, natoms, self.nbeads,
                                  nmolecules, molecule_lengths,
                                  self.getActions())
        
        #call this to clean up after ourselves
        self.destroyOmmSystem()

	
        
        
        

        
        
