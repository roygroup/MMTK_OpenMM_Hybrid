/* Low-level Langevin dynamics integrators
 *
 * Written by Konrad Hinsen
 * Modified by Stephen Constable, Nabil Faruk and Kevin Bishop
 * Roy Group
 * 2013
 */

#include "MMTK/universe.h"
#include "MMTK/forcefield.h"
#include "MMTK/trajectory.h"
#include "OpenMMCWrapper.h"
#include "RPMDOpenMMCWrapper.h"

static void PI_fold_into_box(vector3 *x, int nmolecules, int *molecule_lengths, int natoms, int nbeads, double *data);

/* Global variables */
OpenMM_System* omm_system;
OpenMM_BondArray* bonded_list;
double kB;                      /* Boltzman constant */

/* Allocate and initialize Output variable descriptors */

static PyTrajectoryVariable * get_data_descriptors(int n, 
                                                  PyUniverseSpecObject *universe_spec,
                                            		  PyArrayObject *configuration, PyArrayObject *velocities,
                                            		  PyArrayObject *gradients, PyArrayObject *masses,
                                            		  int *ndf, double *time,
                                            		  double *p_energy, double *k_energy,
                                            		  double *temperature, double *pressure,
                                            		  double *box_size)
{
  PyTrajectoryVariable *vars = (PyTrajectoryVariable *)
                               malloc((n+1)*sizeof(PyTrajectoryVariable));
  int i = 0;
  if (vars != NULL) {
    if (time != NULL && i < n) {
      vars[i].name = "time";
      vars[i].text = "Time: %lf\n";
      vars[i].unit = time_unit_name;
      vars[i].type = PyTrajectory_Scalar;
      vars[i].class = PyTrajectory_Time;
      vars[i].value.dp = time;
      i++;
    }
    if (p_energy != NULL && i < n) {
      vars[i].name = "potential_energy";
      vars[i].text = "Potential energy: %lf, ";
      vars[i].unit = energy_unit_name;
      vars[i].type = PyTrajectory_Scalar;
      vars[i].class = PyTrajectory_Energy;
      vars[i].value.dp = p_energy;
      i++;
    }
    if (k_energy != NULL && i < n) {
      vars[i].name = "kinetic_energy";
      vars[i].text = "Kinetic energy: %lf\n";
      vars[i].unit = energy_unit_name;
      vars[i].type = PyTrajectory_Scalar;
      vars[i].class = PyTrajectory_Energy;
      vars[i].value.dp = k_energy;
      i++;
    }
    if (temperature != NULL && i < n) {
      vars[i].name = "temperature";
      vars[i].text = "Temperature: %lf\n";
      vars[i].unit = temperature_unit_name;
      vars[i].type = PyTrajectory_Scalar;
      vars[i].class = PyTrajectory_Thermodynamic;
      vars[i].value.dp = temperature;
      i++;
    }
    if (pressure != NULL && i < n) {
      vars[i].name = "pressure";
      vars[i].text = "Pressure: %lf\n";
      vars[i].unit = pressure_unit_name;
      vars[i].type = PyTrajectory_Scalar;
      vars[i].class = PyTrajectory_Thermodynamic;
      vars[i].value.dp = pressure;
      i++;
    }
    if (configuration != NULL && i < n) {
      vars[i].name = "configuration";
      vars[i].text = "Configuration:\n";
      vars[i].unit = length_unit_name;
      vars[i].type = PyTrajectory_ParticleVector;
      vars[i].class = PyTrajectory_Configuration;
      vars[i].value.array = configuration;
      i++;
    }
    if (box_size != NULL && i < n) {
      vars[i].name = "box_size";
      vars[i].text = "Box size:";
      vars[i].unit = length_unit_name;
      vars[i].type = PyTrajectory_BoxSize;
      vars[i].class = PyTrajectory_Configuration;
      vars[i].value.dp = box_size;
      vars[i].length = universe_spec->geometry_data_length;
      i++;
  }
    if (velocities != NULL && i < n) {
      vars[i].name = "velocities";
      vars[i].text = "Velocities:\n";
      vars[i].unit = velocity_unit_name;
      vars[i].type = PyTrajectory_ParticleVector;
      vars[i].class = PyTrajectory_Velocities;
      vars[i].value.array = velocities;
      i++;
    }
    if (gradients != NULL && i < n) {
      vars[i].name = "gradients";
      vars[i].text = "Energy gradients:\n";
      vars[i].unit = energy_gradient_unit_name;
      vars[i].type = PyTrajectory_ParticleVector;
      vars[i].class = PyTrajectory_Gradients;
      vars[i].value.array = gradients;
      i++;
    }
    if (masses != NULL && i < n) {
      vars[i].name = "masses";
      vars[i].text = "Masses:\n";
      vars[i].unit = mass_unit_name;
      vars[i].type = PyTrajectory_ParticleScalar;
      vars[i].class = PyTrajectory_Internal;
      vars[i].value.array = masses;
      i++;
    }
    if (ndf != NULL && i < n) {
      vars[i].name = "degrees_of_freedom";
      vars[i].text = "Degrees of freedom: %d\n";
      vars[i].unit = "";
      vars[i].type = PyTrajectory_IntScalar;
      vars[i].class = PyTrajectory_Internal;
      vars[i].value.ip = ndf;
      i++;
    }
    vars[i].name = NULL;
  }
  return vars;
}

/* Langevin dynamics integrator */

static PyObject * integrateLD(PyObject *dummy, PyObject *args)
{
  /* The parameters passed from the Python code */
  PyObject *universe;                 /* universe */
  PyArrayObject *configuration;       /* array of positions */
  PyArrayObject *velocities;          /* array of velocities */
  PyArrayObject *masses;              /* array of masses */
  PyListObject *spec_list;            /* list of periodic actions */
  double friction;                    /* friction coefficient */
  double ext_temp;                    /* temperature of the heat bath */
  double delta_t;                     /* time step */
  int steps;                          /* number of steps */
  int skip;							              /* number of steps to skip */
  int natoms;                         /* number of atoms */
  int nbeads;						              /* number of beads per atom */
  int nmolecules;                     /* number of molecules */
  PyArrayObject *molecule_lengths;    /* array of molecule lengths */

  /* Other variables, see below for explanations */
  PyThreadState *this_thread;
  PyUniverseSpecObject *universe_spec;
  PyArrayObject *gradients, *random1, *random2;
  PyTrajectoryOutputSpec *output;
  PyTrajectoryVariable *data_descriptors = NULL;
  vector3 *x, *v, *g;
  double *f;
  energy_data p_energy;
  double time, k_energy, temperature;
  int atoms, df;
  int pressure_available;
  int i, j, k;
  
  /* OpenMM required variables */
  OpenMM_Context* context;
  OpenMM_Platform* platform;
  OpenMM_RPMDIntegrator* integrator;
  OpenMM_State* state;
  OpenMM_Vec3Array* omm_pos;
  OpenMM_Vec3Array* omm_vel;
  OpenMM_Vec3* pos_i;
  OpenMM_Vec3* vel_i;

  /* Get arguments passed from Python code */
  if (!PyArg_ParseTuple(args, "OO!O!O!dddiiiiiO!O!", &universe,
			&PyArray_Type, &configuration,
			&PyArray_Type, &velocities,
			&PyArray_Type, &masses,
			&friction, &ext_temp, &delta_t, &steps, 
			&skip, &natoms, &nbeads, &nmolecules,
      &PyArray_Type, &molecule_lengths,
      &PyList_Type, &spec_list))
    return NULL;

  /* Obtain the universe specification */
  universe_spec = (PyUniverseSpecObject *)
                   PyObject_GetAttrString(universe, "_spec");
  if (universe_spec == NULL)
    return NULL;

  /* Create the array for energy gradients */
#if defined(NUMPY)
  gradients = (PyArrayObject *)PyArray_Copy(configuration);
#else
  gradients = (PyArrayObject *)PyArray_FromDims(configuration->nd,
						configuration->dimensions,
						PyArray_DOUBLE);
#endif
  if (gradients == NULL)
    return NULL;
    
    

  /* Set some convenient variables */
  atoms = natoms;                          /* number of atoms */
  x = (vector3 *)configuration->data;      /* pointer to positions */
  v = (vector3 *)velocities->data;         /* pointer to velocities */
  g = (vector3 *)gradients->data;          /* pointer to energy gradients */
  df = 3*atoms; 						   /* number of degrees of freedom */
  
    /* Initialize trajectory output and periodic actions */
  data_descriptors =
    get_data_descriptors(10 + pressure_available,
			 universe_spec,
			 configuration, velocities,
			 gradients, masses, &df,
			 &time, &p_energy.energy, &k_energy,
			 &temperature, NULL,
			 (universe_spec->geometry_data_length > 0) ?
			           universe_spec->geometry_data : NULL);
  if (data_descriptors == NULL){
    goto error2;}
  output = PyTrajectory_OutputSpecification(universe, spec_list,
					    "Langevin Dynamics",
					    data_descriptors);
  if (output == NULL){
    goto error2;}

  /* Initial coordinate correction (for periodic universes etc.) */
  if (universe_spec->is_periodic == 1) {
    PI_fold_into_box(x, nmolecules, molecule_lengths, natoms, nbeads, universe_spec->geometry_data);
  }
  
  if (universe_spec->is_periodic == 1) {
    if (universe_spec->is_orthogonal == 1) {
    	/*double min =  universe_spec->geometry_data[0];
		int geometry_counter = 0;
		for (geometry_counter = 0; geometry_counter < universe_spec->geometry_data_length; geometry_counter++) {
			if (universe_spec->geometry_data[geometry_counter] < min) {
		  		min = universe_spec->geometry_data[geometry_counter];
			}
		}*/

      OpenMM_Vec3 a = {universe_spec->geometry_data[0],0.0,0.0};
      OpenMM_Vec3 b = {0.0,universe_spec->geometry_data[1],0.0};
      OpenMM_Vec3 c = {0.0,0.0,universe_spec->geometry_data[2]};
      
      //printf("after creation of box %lf \n",b.y);
  
  		OpenMM_System_setDefaultPeriodicBoxVectors(omm_system, &a, &b, &c);
  		
  		printf("REMARK Using OpenMM periodic box with dimensions (%f, %f, %f)\n", universe_spec->geometry_data[0], universe_spec->geometry_data[1], universe_spec->geometry_data[2]);
  		
  	} else {
  		/* abort! OpenMM does not support non-orthorhombic periodic universes */
  		printf("ERROR OpenMM does not support non-orthorhombic universes\n");
  		goto error2;
   	}
  }

  integrator = (OpenMM_Integrator*) OpenMM_RPMDIntegrator_create(nbeads, ext_temp, friction, delta_t);

  /* Let OpenMM Context choose best platform. */
/*  context = OpenMM_Context_create(omm_system, integrator);
  platform = OpenMM_Context_getPlatform(context);*/
  platform = OpenMM_Platform_getPlatformByName("Reference"); 
  context = OpenMM_Context_create_2(omm_system, integrator, platform); 
  printf( "REMARK  Using OpenMM platform %s\n", 
      OpenMM_Platform_getName(platform));
  /*printf("here is the default precision %s\n",
      OpenMM_Platform_getPropertyDefaultValue(platform,"ReferencePrecision"));*/

  /* Set starting positions and velocities of the atoms. Leave time zero. */
  
  
  /*printf("atoms: %d\n", atoms);*/

  for (j = 0; j < nbeads; j++) {
    omm_pos = OpenMM_Vec3Array_create(atoms);
    omm_vel = OpenMM_Vec3Array_create(atoms);
    for (i =  0; i < atoms; i++) { 
      pos_i = x[j+i*nbeads];
      vel_i = v[j+i*nbeads];

      if (universe_spec->is_periodic == 1) {
        pos_i->x = pos_i->x + universe_spec->geometry_data[0]/2.;
        pos_i->y = pos_i->y + universe_spec->geometry_data[1]/2.;
        pos_i->z = pos_i->z + universe_spec->geometry_data[2]/2.;
      }

      OpenMM_Vec3Array_set(omm_pos, i, *pos_i);
      OpenMM_Vec3Array_set(omm_vel, i, *vel_i);
    } 
  	OpenMM_RPMDIntegrator_setPositions(integrator, j, omm_pos);
  	OpenMM_RPMDIntegrator_setVelocities(integrator, j, omm_vel);
  }

  /* Collect properties of interest for which to query OpenMM */
  int infoMask = OpenMM_State_Positions + OpenMM_State_Velocities + OpenMM_State_Energy;
  
  /*
   * Main integration loop
   */
  time  = 0.;
  
  printf("completing %d steps while skipping %d per log output\n", steps, skip);

/*  OpenMM_Vec3 ap;
  OpenMM_Vec3 bp;
  OpenMM_Vec3 cp;
  state = OpenMM_RPMDIntegrator_getState(integrator, 0, infoMask, 0);
  OpenMM_State_getPeriodicBoxVectors(state, &ap, &bp, &cp);
  printf("The y dimension of the periodic box is %f \n", bp.y);
  OpenMM_State_destroy(state);*/
  
  for (k = 0; k < steps; k+= skip) {
    /* Calculation of thermodynamic properties */
    p_energy.energy = 0.;
    k_energy = 0.;
    for (j = 0; j < nbeads; j++) {
    	state = OpenMM_RPMDIntegrator_getState(integrator, j, infoMask, 0);
    	omm_pos = OpenMM_State_getPositions(state);
      omm_vel = OpenMM_State_getVelocities(state);
    	int omm_index = 0;
    	for (i = j; i < atoms*nbeads; i+= nbeads) {
    		pos_i = OpenMM_Vec3Array_get(omm_pos, omm_index);
        vel_i = OpenMM_Vec3Array_get(omm_vel, omm_index);

        if (universe_spec->is_periodic == 1) {
          pos_i->x = pos_i->x - universe_spec->geometry_data[0]/2.;
          pos_i->y = pos_i->y - universe_spec->geometry_data[1]/2.;
          pos_i->z = pos_i->z - universe_spec->geometry_data[2]/2.;
        }

        x[i][0] = pos_i->x;
    		x[i][1] = pos_i->y;
    		x[i][2] = pos_i->z;

        v[i][0] = vel_i->x;
        v[i][1] = vel_i->y;
        v[i][2] = vel_i->z;

    		omm_index++;
    	}
    	k_energy += OpenMM_State_getKineticEnergy(state);
    	p_energy.energy += OpenMM_State_getPotentialEnergy(state);
    	time = OpenMM_State_getTime(state);
    	OpenMM_State_destroy(state);
    }
    k_energy/=(double)nbeads;
    p_energy.energy /=(double)nbeads;
    temperature = 2.*k_energy/(df*kB);


    /* Coordinate correction (for periodic universes etc.) */
    if (universe_spec->is_periodic == 1) {
      PI_fold_into_box(x, nmolecules, molecule_lengths, natoms, nbeads, universe_spec->geometry_data);
    } 

    /* Trajectory and log output */
    PyTrajectory_Output(output, k, data_descriptors, NULL);
    
    OpenMM_RPMDIntegrator_step(integrator, skip);

  }
  /** End of main integration loop **/

  /* Final thermodynamic property evaluation */
  /* Calculation of thermodynamic properties */
  p_energy.energy = 0.;
  k_energy = 0.;
  for (j = 0; j < nbeads; j++) {
    state = OpenMM_RPMDIntegrator_getState(integrator, j, infoMask, 0);
    omm_pos = OpenMM_State_getPositions(state);
    omm_vel = OpenMM_State_getVelocities(state);
    int omm_index = 0;
    for (i = j; i < atoms*nbeads; i+= nbeads) {
      pos_i = OpenMM_Vec3Array_get(omm_pos, omm_index);
      vel_i = OpenMM_Vec3Array_get(omm_vel, omm_index);

      if (universe_spec->is_periodic == 1) {
        pos_i->x = pos_i->x - universe_spec->geometry_data[0]/2.;
        pos_i->y = pos_i->y - universe_spec->geometry_data[1]/2.;
        pos_i->z = pos_i->z - universe_spec->geometry_data[2]/2.;
      }

      x[i][0] = pos_i->x;
      x[i][1] = pos_i->y;
      x[i][2] = pos_i->z;

      v[i][0] = vel_i->x;
      v[i][1] = vel_i->y;
      v[i][2] = vel_i->z;

      omm_index++;
    }
    k_energy += OpenMM_State_getKineticEnergy(state);
    p_energy.energy += OpenMM_State_getPotentialEnergy(state);
    time = OpenMM_State_getTime(state);
    OpenMM_State_destroy(state);
  }
  k_energy/=(double)nbeads;
  p_energy.energy /=(double)nbeads;
  temperature = 2.*k_energy/(df*kB);

  /* Coordinate correction (for periodic universes etc.) */
  if (universe_spec->is_periodic == 1) {
    PI_fold_into_box(x, nmolecules, molecule_lengths, natoms, nbeads, universe_spec->geometry_data);
  }

  /* Final trajectory and log output */
  PyTrajectory_Output(output, k, data_descriptors, NULL);
  
  /* Cleanup */
  PyUniverseSpec_StateLock(universe_spec, -2);
  /*PyEval_RestoreThread(this_thread);*/
  PyTrajectory_OutputFinish(output, k, 0, 1, data_descriptors);
  free(data_descriptors);
  Py_DECREF(gradients);
  Py_INCREF(Py_None);
  return Py_None;

  /* Error return */
error:
  PyTrajectory_OutputFinish(output, k, 1, 1, data_descriptors);
error2:
  free(data_descriptors);
  Py_DECREF(gradients);
  return NULL;
}


static PyObject* initOpenMM(PyObject* dummy, PyObject* args) {
	PyArrayObject* masses;
	int natoms;
		
	if (!PyArg_ParseTuple(args, "O!i",
			&PyArray_Type, &masses,
			&natoms))
	return NULL;
		
	double* m = (double*)masses->data;
	omm_system = OpenMM_System_create();
		
	int iter;
	for (iter = 0; iter < natoms; iter++) {
		OpenMM_System_addParticle(omm_system, m[iter]);
		/* printf("Added particle %d with mass %f\n", iter, m[iter]); */
	}
	return Py_None;
}

static PyObject* destroyOpenMM(PyObject* dummy, PyObject* args) {
	OpenMM_System_destroy(omm_system);
	OpenMM_BondArray_destroy(bonded_list);
	return Py_None;
}

static PyObject* addCMRemover(PyObject *dummy, PyObject *args) {
	int skip;
	
	if (!PyArg_ParseTuple(args, "i",
		&skip))
	return NULL;
	
	OpenMM_CMMotionRemover* remover = OpenMM_CMMotionRemover_create(skip);
	OpenMM_System_addForce(omm_system, (OpenMM_Force*)remover);
	return Py_None;
}

static PyObject*  makeOmmBondedForce(PyObject *dummy, PyObject *args) {
	PyArrayObject *atom_index_1;
	PyArrayObject *atom_index_2;
	PyArrayObject *eq_distance;
	PyArrayObject *spring_k;
	int nbonds;
	
	if (!PyArg_ParseTuple(args, "O!O!O!O!i",
			&PyArray_Type, &atom_index_1,
			&PyArray_Type, &atom_index_2,
			&PyArray_Type, &eq_distance,
			&PyArray_Type, &spring_k,
			&nbonds))
    return NULL;
    
    int* i = (int*)atom_index_1->data;
    int* j = (int*)atom_index_2->data;
    double* r = (double*)eq_distance->data;
    double* k_eq = (double*)spring_k->data;
    
    OpenMM_HarmonicBondForce* bonded = OpenMM_HarmonicBondForce_create();
    bonded_list = OpenMM_BondArray_create(nbonds);
    
    printf("Adding bonds\n");
    
    int iter;
    for (iter = 0; iter< nbonds; iter++) {
    	OpenMM_HarmonicBondForce_addBond(bonded, i[iter], j[iter], r[iter], k_eq[iter]);
    	/*printf("Added bond %d between atoms %d and %d with length %f and strength %f\n", iter, i[iter], j[iter], r[iter], k_eq[iter]);*/
    	OpenMM_BondArray_append(bonded_list, i[iter], j[iter]);
    }
    
    OpenMM_System_addForce(omm_system, (OpenMM_Force*)bonded);
    return Py_None;
}

static PyObject* makeOmmAngleForce(PyObject *dummy, PyObject *args) {
	PyArrayObject *atom_index_1;
	PyArrayObject *atom_index_2;
	PyArrayObject *atom_index_3;
	PyArrayObject *eq_angle;
	PyArrayObject *spring_k;
	int nangles;
	
	if (!PyArg_ParseTuple(args, "O!O!O!O!O!i",
			&PyArray_Type, &atom_index_1,
			&PyArray_Type, &atom_index_2,
			&PyArray_Type, &atom_index_3,
			&PyArray_Type, &eq_angle,
			&PyArray_Type, &spring_k,
			&nangles))
    return NULL;
    
    int* i = (int*)atom_index_1->data;
    int* j = (int*)atom_index_2->data;
    int* k = (int*)atom_index_3->data;
    double* theta = (double*)eq_angle->data;
    double* k_eq = (double*)spring_k->data;
    
    OpenMM_HarmonicAngleForce* angles = OpenMM_HarmonicAngleForce_create();
    
    printf("Adding angles\n");
    
    int iter;
    for (iter = 0; iter < nangles; iter++) {
    	OpenMM_HarmonicAngleForce_addAngle(angles, i[iter], j[iter], k[iter], theta[iter], k_eq[iter]);
    }
    
    OpenMM_System_addForce(omm_system, (OpenMM_Force*)angles);
    return Py_None;
}    
    
static PyObject* makeOmmDihedralForce(PyObject* dummy, PyObject* args) {
	PyArrayObject* atom_index_1;
	PyArrayObject* atom_index_2;
	PyArrayObject* atom_index_3;
	PyArrayObject* atom_index_4;
	PyArrayObject* periodicity;
	PyArrayObject* eq_dihedral;
	PyArrayObject* spring_k;
	int ndihedrals;
	
	if (!PyArg_ParseTuple(args, "O!O!O!O!O!O!O!i",
			&PyArray_Type, &atom_index_1,
			&PyArray_Type, &atom_index_2,
			&PyArray_Type, &atom_index_3,
			&PyArray_Type, &atom_index_4,
			&PyArray_Type, &periodicity,
			&PyArray_Type, &eq_dihedral,
			&PyArray_Type, &spring_k,
			&ndihedrals))
    return NULL;
    
    int* i = (int*)atom_index_1->data;
    int* j = (int*)atom_index_2->data;
    int* k = (int*)atom_index_3->data;
    int* l = (int*)atom_index_4->data;
    int* n = (int*)periodicity->data;
    double* phi = (double*)eq_dihedral->data;
    double* k_eq = (double*)spring_k->data;
    
    OpenMM_PeriodicTorsionForce* dihedrals = OpenMM_PeriodicTorsionForce_create();
    
    printf("Adding dihedrals\n");
    
    int iter;
    for (iter = 0; iter < ndihedrals; iter++) {
    	OpenMM_PeriodicTorsionForce_addTorsion(dihedrals, i[iter], j[iter], k[iter], l[iter], n[iter], phi[iter], k_eq[iter]);
    }
    
    OpenMM_System_addForce(omm_system, (OpenMM_Force*)dihedrals);
    return Py_None;
}

static PyObject* makeOmmEsAndLjForce(PyObject *dummy, PyObject *args) {
  int is_periodic;
	PyArrayObject* sigma;
	PyArrayObject* epsilon;
	PyArrayObject* charge;
	int natoms;
	double es_14_factor, lj_14_factor, cutoff;

  PyUniverseSpecObject *universe_spec;
	
	if (!PyArg_ParseTuple(args, "iO!O!O!dddi",
      &is_periodic,
			&PyArray_Type, &sigma,
			&PyArray_Type, &epsilon,
			&PyArray_Type, &charge,
			&es_14_factor, &lj_14_factor,
			&cutoff, &natoms))
	return NULL;
	
	double* s = (double*)sigma->data;
	double* e = (double*)epsilon->data;
	double* z = (double*)charge->data;
	
	OpenMM_NonbondedForce* nonbond = OpenMM_NonbondedForce_create();
	
	printf("Adding nonbonded\n");
	
	int iter;
	for (iter = 0; iter < natoms; iter++) {
		OpenMM_NonbondedForce_addParticle(nonbond,  z[iter], s[iter], e[iter]);
	}
	
	OpenMM_NonbondedForce_createExceptionsFromBonds(nonbond, bonded_list, es_14_factor, lj_14_factor);
	
  if (is_periodic == 1) {
    OpenMM_NonbondedForce_setNonbondedMethod(nonbond, OpenMM_NonbondedForce_Ewald); /*for step_0*/
    OpenMM_NonbondedForce_setCutoffDistance(nonbond, cutoff);
  }
  else {
    OpenMM_NonbondedForce_setNonbondedMethod(nonbond, OpenMM_NonbondedForce_NoCutoff); /*for example*/
  }
	
	OpenMM_System_addForce(omm_system, (OpenMM_Force*)nonbond);
	return Py_None;
}

static PyMethodDef langevin_methods[] = {
  {"integrateLD", integrateLD, 1},
  {"initOpenMM", initOpenMM, 1},
  {"destroyOpenMM", destroyOpenMM, 1},
  {"addCMRemover", addCMRemover, 1},
  {"makeOmmBondedForce", makeOmmBondedForce, 1},
  {"makeOmmAngleForce", makeOmmAngleForce, 1},
  {"makeOmmDihedralForce", makeOmmDihedralForce, 1},
  {"makeOmmEsAndLjForce", makeOmmEsAndLjForce, 1},
  {NULL, NULL}		/* sentinel */
};

/* Initialization function for the module */

void
initMMTK_langevin()
{
  PyObject *m, *dict;
  PyObject *universe, *trajectory, *forcefield, *units;

  /* Create the module and add the functions */
  m = Py_InitModule("MMTK_langevin", langevin_methods);
  dict = PyModule_GetDict(m);

  /* Import the array module */
  import_array();

  /* Import MMTK modules */
  import_MMTK_universe();
  import_MMTK_forcefield();
  import_MMTK_trajectory();
  
  /* OpenMM library initialization */
  OpenMM_StringArray* pluginList = OpenMM_Platform_loadPluginsFromDirectory(
        OpenMM_Platform_getDefaultPluginsDirectory());
  int num_plugins = OpenMM_StringArray_getSize(pluginList);
  int i;
  for (i = 0; i < num_plugins; i++) {
  	printf("loaded plugin %s\n", OpenMM_StringArray_get(pluginList, i));
  }
  OpenMM_StringArray_destroy(pluginList);

  /* Get the Boltzman constant factor from MMTK.Units */
  units = PyImport_ImportModule("MMTK.Units");
  if (units != NULL) {
    PyObject *module_dict = PyModule_GetDict(units);
    PyObject *factor = PyDict_GetItemString(module_dict, "k_B");
    kB = PyFloat_AsDouble(factor);
  }

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module MMTK_langevin");
}

static void
PI_fold_into_box(vector3 *x, int nmolecules, int *molecule_lengths, int natoms, int nbeads, double *data)
{
  double a = data[0];
  double b = data[1];
  double c = data[2];
  double ah = 0.5*a, bh = 0.5*b, ch = 0.5*c;
  double center[3];
  int i, j, k, *l, m;

  k = 0;
  m = 0;
  for (i = 0; i < nmolecules; i++) {
    /* Find the molecule centroid. */
    memset(center, 0, 3);
    l = (int *)PyArray_GETPTR1(molecule_lengths,i);

    for (j = k; j < k+(*l)*nbeads; j++) {
      center[0] += x[j][0];
      center[1] += x[j][1];
      center[2] += x[j][2];
    }
    k = j;

    center[0] *= 1.0/((*l)*nbeads);
    center[1] *= 1.0/((*l)*nbeads);
    center[2] *= 1.0/((*l)*nbeads);

    /* reposition centroid in terms of a box with origin {0,0,0} */
    center[0] += ah;
    center[1] += bh;
    center[2] += ch;

    /* Find the centroid displacement to move it into the first periodic box. */

    int xcell = (int) floor(center[0]/a);
    int ycell = (int) floor(center[1]/b);
    int zcell = (int) floor(center[2]/c);
    double dx = xcell*a;
    double dy = ycell*b;
    double dz = zcell*c;

    /* Translate all the beads in the molecule. */
    for (j = m; j < m+(*l)*nbeads; j++) {
      x[j][0] -= dx;
      x[j][1] -= dy;
      x[j][2] -= dz;
    }
    m = j;
  }
}

