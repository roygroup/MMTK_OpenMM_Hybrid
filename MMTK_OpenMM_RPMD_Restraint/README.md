Before running the example scripts here, you will need to have installed MMTK and OpenMM as described in the upper directory. The installation guide is available at MMTK_OpenMM_Hybrid/InstallationGuide.pdf

After installation, navigate to this directory (MMTK_OpenMM_RPMD_Restraint)

THe first step is to compile the C code. Do this by entering:
    pygpu setup.py build_ext --inplace

Note: pygpu will be whatever alias you have set for your development python

Now, you can run the example script by entering:

    pygpu spcfw-q_restraint.py 1 10.0 1.0 300 1.0 1000 10000 100

Note: the order of parameters is as follows:
    filename                        (spcfw-q_restraint.py)
    number of beads
    k constant for restraint        (in kJ/mol/nm^2)
    r0 value for restraint          (in nm)
    temperature                     (in K)
    timestep                        (in fs)
    number of equilibration steps
    number of production steps
    number of steps to skip 