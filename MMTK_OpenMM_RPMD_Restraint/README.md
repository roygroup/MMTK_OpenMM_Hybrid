...........................................................................................

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

...........................................................................................


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