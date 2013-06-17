%% LyX 1.6.9 created this file.  For more info, see http://www.lyx.org/.
%% Do not edit unless you really know what you are doing.
\documentclass[11pt,english]{article}
\usepackage{geometry}
\geometry{verbose,tmargin=1in,bmargin=1in,lmargin=1in,rmargin=1in}
\usepackage{xcolor}
\newcommand{\tab}[1]{\hspace{.05\textwidth}\rlap{#1}}
\begin{document}

\section{Installation Guide for MMTK and OpenMM}

This installation guide has been developed to work for Mac OS X Snow Leopard for MMTK 2.7.8 and OpenMM 4.1.1

Created by: Chris Ing and Stephen Constable

Updated by: Nabil Faruk and Kevin Bishop

\subsection{Pre-installation Requirements}

\begin{enumerate}
 \item Xcode
 
 You will require XCode to be installed on your Mac. It can be found in the additional installs directory of the discs that come with your Mac. X11 is another 	requirement if you need to use GnuPlot.
Note: If you do not have the install discs, you can download Xcode from apple at: \color{blue} https://developer.apple.com/xcode/ \color{black}
 
 \item Development Directory
 
 Create a new directory for Python binaries. A common one is \color{blue} \$HOME/Dev \color{black}. You should not use the default python installation since files in /Library/Python/2.*/site-packages can be modified by Apple patches.
 
\end{enumerate}

\subsection{MMTK Installation}

\begin{enumerate}
 \item Python
 	\begin{itemize}
		\item Download the latest compressed source tarball from \color{blue}http://www.python.org/download of Python 2.* \color{black}
		\item Extract this tarball to a temporary directory with:
 			
		\color{blue} tar -xf Python2.*.tar \color{black}
		
		and cd to that directory
		\item Bootstrap your installation with the following command:
		
		\color{blue} ./configure -prefix=\$HOME/Dev \color{black}
		
		\item Build and install with the command:
		
		\color{blue} make install \color{black}
	\end{itemize}
	Note: You must make your new python accessible from the command line:
	
	\begin{itemize}
		\item Modify your \color{blue} \~{}/.bash\_profile \color{black} script by adding the following line:
		
		\color{blue} alias pydev='"\$HOME/Dev/bin/python" \$*' \color{black}
		\item Then reload your bash\_profile by entering:
		
		\color{blue} source \~{}/.bash\_profile \color{black}
	\end{itemize}
	
	Now, when you type \color{blue}pydev \color{black} at the command line, your new python will open.
	
 \item Cython
 \begin{itemize}
 	\item Download the latest tar file and extract by entering:
	
	\color{blue} tar -xf Cython-0.17.2.tar.gz  \color{black}
	\item Install Cython by entering (inside the Cython directory):
	
	\color{blue} pydev setup.py install  \color{black}
	
	Ensure that you use \color{blue}pydev  \color{black}as we do not want to modify the original python environment. This will be critical for the rest of the future installations as well.

 \end{itemize}
 
 \item zlib
  \begin{itemize}
 	\item Download zlib 1.2.5 from \color{blue} ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/ \color{black}
	\item Extract by entering:
	
	\color{blue} tar -xf zlib-1.2.5.tar.gz  \color{black}
	\item Build and install by entering the following 2 lines:
	
	\color{blue} ./configure -prefix=/Users/kpbishop/Dev 
	
	make check install \color{black}
	
	Note: \color{blue} kpbishop \color{black} should be your own username  

 \end{itemize}
 
 \item HDF5
  \begin{itemize}
 	\item Download HDF5 1.8.6 from \color{blue} ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4/ \color{black}
	\item Extract by entering:
	
	\color{blue} tar -xf hdf5-1.8.6.tar.gz  \color{black}
	
	\item Build and install by entering the following 2 lines:
	
	\color{blue} ./configure -with-zlib=/Users/kpbishop/Dev -prefix=/Users/kpbishop/Dev
	
	make check install \color{black}
 \end{itemize}
 
 \item netCDF
  \begin{itemize}
 	\item Download the latest C netCDF source from:
	
	\color{blue} http://www.unidata.ucar.edu/downloads/netcdf/index.jsp \color{black}
	
	netCDF used to come for C, C++ and Fortran but we only need the C portion. The most recent release at time of writing is the netCDF C library and utilities, version 4.2.1.1
	\item Extract this file by entering:
	
	\color{blue} tar -xf netcdf4.*.tar.gz \color{black}
 	
	\item Build and install netcdf and specify the locations of zlib and HDF5 with the following 4 commands:
	
	\color{blue} CPPFLAGS=-I/Users/kpbishop/Dev/include 
	
	LDFLAGS=-L/Users/kpbishop/Dev/lib 
	
	 ./configure -prefix=/Users/kpbishop/Dev
	
	 make check install  \color{black}
 \end{itemize}
 
 \item NumPy
  \begin{itemize}
 	\item Download the latest version of NumPy as a tar.gz (currently 1.6.2)
	\item Extract it using the following line:
	
	\color{blue} tar -xf numpy-1.6.2.tar.gz		\color{black}
	\item To build and install NumPy, enter:
	
	\color{blue} pydev setup.py install		\color{black}

 \end{itemize}
 
 \item Scientific Python
  \begin{itemize}
 	\item Download the latest version of Scientific Python as a tar.gz (currently 2.9.1) from \color{blue}http://sourcesup.cru.fr/projects/scientific-py/ \color{black}
	\item Extract is using the following line:
	
	\color{blue} tar -xf ScientificPython-2.9.1.tar		\color{black}
	\item To build and install SciPy, enter:
	
\color{blue}	pydev setup.py install				\color{black}
	\item Test your installation using any of the scripts in the Examples folder, try protein.py in /Examples/MolecularDynamics/protein.py with:
	
	\color{blue} pydev protein.py					\color{black}

 \end{itemize}
 
 \item MMTK
  \begin{itemize}
 	\item Download the latest development release of MMTK from:
	
\color{blue} http://sourcesup.cru.fr/projects/mmtk		\color{black}
	\item Extract MMTK using:
	
	\color{blue} tar -xf MMTK-2.*.tar				\color{black}
 	\item Build and install MMTK with: (This setup file may need to be changed depending on your install path)
	
	\color{blue} pydev setup.py install				\color{black}
  \end{itemize}
  
 \newpage 
  
 Notes:
 \begin{itemize}
	\item To run the Langevin Dynamics example script in version 2.7.1 of MMTK, you may need to manually add the netCDF include directory to the setup.py script by changing Line 20:
 
	\color{blue} include\_dirs=['./'])			\color{black}
 
changed to:
 
	\color{blue} include\_dirs=['./','/Users/kpbishop/Dev/include'])  \color{black}
	\item then the Langevin integrator can be compiled with:
 
 	\color{blue} pydev setup.py build\_ext -inplace			\color{black}
	\item and the example can be run with:
	
	\color{blue} pydev example.py			\color{black}

	\item In some cases fftw may be required to be installed. 
	
	\begin{itemize}
	
	\item You can obtain a tarball from \color{blue} http://www.fftw.org/download.html\color{black}  
	
	and extract using:
	
	\color{blue} tar -xf fftw-3.3.3.tar.gz	\color{black}
	\item Configure using:
	
	\color{blue} ./configure --prefix=\$HOME/Dev --enable-shared		\color{black}
	\item and build and install using:
	
	\color{blue} make install		\color{black}
	\end{itemize}
 \end{itemize}
\end{enumerate}

\subsection{OpenMM Installation}

\begin{enumerate}
 \item Copy over a working OpenMM directory:

	\color{blue} /home/nffaruk/Sugar\_GPU/OpenMM4.1.1-Source/ \color{black}
	
	or try from the OpenMM website 
	
	(\color{blue} https://simtk.org/project/xml/downloads.xml?group\_id=161\color{black})
 \item You will need CMake for the installation process. Download and install it from

	\color{blue} http://www.cmake.org/files/v2.8/cmake-2.8.9-Darwin64-universal.dmg \color{black}
 \item GCC-XML
 \begin{itemize}
 \item Copy to your machine by entering:
 
 \color{blue} git clone git://github.com/gccxml/gccxml.git \color{black}
 
 \item Install it with CMake by cding into the gccxml source directory and using the command:
 
 \color{blue} /usr/bin/cmake -i  \color{black}
 
Choose your Dev folder for the install path when prompted
 \end{itemize}
 
 \item Get and install CUDA toolkit, drivers, and SDK for Mac:

\color{blue} http://developer.nvidia.com/cuda/cuda-downloads \color{black}

May also need CUDA driver from:

\color{blue} http://www.nvidia.com/object/mac-driver-archive.html \color{black}

 \item Update your \color{blue} bash\_profile \color{black} so that it has the following lines at the bottom:

\color{blue}

export MMTK\_USE\_CYTHON=1

export DYLD\_LIBRARY\_PATH=\$HOME/Dev/openmm/lib

export DYLD\_LIBRARY\_PATH= 

\tab\tab \${DYLD\_LIBRARY\_PATH}:\$HOME/Dev/openmm/lib/plugins

export DYLD\_LIBRARY\_PATH=\${DYLD\_LIBRARY\_PATH}:\$HOME/Dev/lib

export DYLD\_LIBRARY\_PATH=\${DYLD\_LIBRARY\_PATH}:/usr/local/cuda/lib

export OPENCL\_DIR=/System/Library/Frameworks/OpenCL.framework

export OPENMM\_PLUGIN\_DIR=\$HOME/Dev/openmm/lib/plugins

\color{black}

Reload it using:

\color{blue} source \~{}/.bash\_profile \color{black}

\item OpenMM

\begin{itemize} 

\item You should be able to install OpenMM now by cding to the top level of its source directory and use CMake:

\color{blue} /usr/bin/cmake -i \color{black}

In the CMake prompts watch out for:

\begin{itemize}

\item the installation directory... should be \color{blue} \~{}/Dev/openmm \color{black}

\item where gcc-xml executable is located.. should be \color{blue}\~{}/Dev/bin/gccxml \color{black}

\item whether to build C and fortran wrappers... \color{blue}yes \color{black}
\item whether to build python wrapper.. \color{blue}no \color{black}

\end{itemize}

\item Then build/install with:

\color{blue} make OpenMM \color{black}

\color{blue} make install \color{black}

\item Then test with:

\color{blue} make test \color{black}

Remember to \color{blue} 'make clean' \color{black} if you experience any problems and want to reinstall after making modifications

\end{itemize}

\end{enumerate}

\end{document}
