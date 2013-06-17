#Setup file to compile C executables
#Last Modified: June 17, 2013

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

from distutils.core import setup, Extension
import os, sys

compile_args = []
library_dirs = ['/home/kpbishop/Dev_GPU/openmm5/lib', '/home/kpbishop/Dev_GPU/openmm5/lib/plugins']
libraries = ['OpenMM', 'OpenMMRPMD','OpenMMOpenCL', 'OpenMMRPMDOpenCL', 'OpenMMCUDA', 'OpenMMRPMDCUDA']
include_dirs = [
'/home/kpbishop/Dev_GPU/openmm5/include',
'/home/kpbishop/Dev_GPU/include']

from Scientific import N
try:
    num_package = N.package
except AttributeError:
    num_package = "Numeric"
if num_package == "NumPy":
    compile_args.append("-DNUMPY=1")
    if sys.platform == 'win32':
        include_dirs.append(os.path.join(sys.prefix,
                            "Lib/site-packages/numpy/core/include"))
    else:
        include_dirs.append(os.path.join(sys.prefix,
                            "lib/python%s.%s/site-packages/numpy/core/include"
                             % sys.version_info [:2]))

setup (name = "MMTK-LangevinDynamics_OpenMM",
       version = "0.1",
       description = "Langevin dynamics module for the Molecular Modelling Toolkit using OpenMM",
       author = "Stephen Constable and Konrad Hinsen",
       author_email = "sjconstable@uwaterloo.ca",
       url = "http://dirac.cnrs-orleans.fr/MMTK/",
       license = "CeCILL-C",

       py_modules = ['LangevinDynamics'],
       ext_modules = [Extension('MMTK_langevin',
                                ['./MMTK_langevin.c'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs,
                                library_dirs=library_dirs,
                                libraries=libraries)]
       )
