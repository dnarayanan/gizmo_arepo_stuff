#-----------------------------------------------------------------
#
# You might be looking for the compile-time Makefile options of the code...
#
# They have moved to a separate file.
#
# To build the code, do the following:
#
#  (1) Copy the file "Template-Config.sh"  to  "Config.sh"
#
#        cp Template-Config.sh Config.sh 
#
#  (2) Edit "Config.sh" as needed for your application
#
#  (3) Run "make"
#
#
#  New compile-time options should be added to the 
#  file "Template-Config.sh" only. Usually, the should be added
#  there in the disabled/default version.
#
#  "Config.sh" should *not* be checked in to the repository
#
#  Note: It is possible to override the default name of the 
#  Config.sh file, if desired, as well as the name of the
#  executable. For example:
#
#   make  CONFIG=MyNewConf.sh  EXEC=GIZMO
# 
#-----------------------------------------------------------------
#
# You might also be looking for the target system SYSTYPE option
#
# It has also moved to a separate file.
#
# To build the code, do the following:
#
# (A) set the SYSTYPE variable in your .bashrc (or similar file):
#
#        e.g. export SYSTYPE=Magny
# or
#
# (B) set SYSTYPE in Makefile.systype 
#     This file has priority over your shell variable.:
#
#    (1) Copy the file "Template-Makefile.systype"  to  "Makefile.systype"
#
#        cp Template-Makefile.systype Makefile.systype 
#
#    (2) Uncomment your system in  "Makefile.systype".
#
# If you add an ifeq for a new system below, also add that systype to
# Template-Makefile.systype
#
###########
#
# This file was originally part of the GADGET3 code developed by
#   Volker Springel (volker.springel@h-its.org). The code has been modified
#   slighty by Phil Hopkins (phopkins@caltech.edu) for GIZMO (mostly 
#   dealing with new files and filename conventions)
#
#############

ifdef SYSTYPE
SYSTYPE := "$(SYSTYPE)"
-include Makefile.systype
else
include Makefile.systype
endif

ifeq ($(wildcard Makefile.systype), Makefile.systype)
INCL = Makefile.systype
else
INCL =
endif
FINCL =


CONFIG   =  Config.sh
PERL     =  /usr/bin/perl

RESULT     := $(shell CONFIG=$(CONFIG) PERL=$(PERL) make -f config-makefile)
CONFIGVARS := $(shell cat GIZMO_config.h)

#HG_COMMIT := $(shell hg id 2>/dev/null)
#HG_REPO := $(shell hg path default)
#HG_BRANCH := $(shell hg branch)
BUILDINFO = "Build on $(HOSTNAME) by $(USER) from $(HG_BRANCH):$(HG_COMMIT) at $(HG_REPO)"
OPT += -DBUILDINFO='$(BUILDINFO)'

ifeq (FIRE_PHYSICS_DEFAULTS,$(findstring FIRE_PHYSICS_DEFAULTS,$(CONFIGVARS)))  # using 'fire default' instead of all the above
    CONFIGVARS += COOLING COOL_LOW_TEMPERATURES COOL_METAL_LINES_BY_SPECIES
    CONFIGVARS += GALSF METALS GALSF_SFR_MOLECULAR_CRITERION GALSF_SFR_VIRIAL_SF_CRITERION=0
    CONFIGVARS += GALSF_FB_GASRETURN GALSF_FB_HII_HEATING GALSF_FB_SNE_HEATING=1 GALSF_FB_RT_PHOTONMOMENTUM
    CONFIGVARS += GALSF_FB_LOCAL_UV_HEATING GALSF_FB_RPWIND_LOCAL GALSF_FB_RPROCESS_ENRICHMENT=4 GALSF_SFR_IMF_VARIATION
endif


CC       = mpicc        # sets the C-compiler (default)
CXX       = mpiCC       # sets the C++-compiler (default)

FC 	 = mpif90

OPTIMIZE = -Wall  -g   # optimization and warning flags (default)

MPICHLIB = -lmpich

GRACKLEINCL =
GRACKLELIBS = -lgrackle


ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
    FFTW_LIBNAMES =  #-lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
    FFTW_LIBNAMES =  #-ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
    FFTW_LIBNAMES =  #-lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif

# we only need fftw if PMGRID is turned on
ifeq (PMGRID, $(findstring PMGRID, $(CONFIGVARS)))
ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
  FFTW_LIBNAMES = -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIBNAMES = -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIBNAMES = -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif
else
# or if POWERSPEC_GRID is activated
ifeq (POWERSPEC_GRID, $(findstring POWERSPEC_GRID, $(CONFIGVARS)))
ifeq (NOTYPEPREFIX_FFTW,$(findstring NOTYPEPREFIX_FFTW,$(CONFIGVARS)))  # fftw installed without type prefix?
  FFTW_LIBNAMES = -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
else
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(CONFIGVARS)))  # test for double precision libraries
  FFTW_LIBNAMES = -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIBNAMES = -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
endif
else
  FFTW_LIBNAMES = #
endif

endif



# RT/RD SYSTEMS ----------------------------------------

ifeq ($(SYSTYPE),"archer")
CC   =  cc
CXX  =  CC
MPICHLIB = 
FC   =  $(CC)
OPT += -DH5_USE_16_API -DCONFIG_BFLOAT_8
OPTIMIZE = -O2 #-hfp3
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -fopenmp  # openmp required compiler flags
endif
FFTW_INCL = #-I/home/d150/d150/rad/gizmo-extra/fftw-2.1.5/include
FFTW_LIBS = #-L/home/d150/d150/rad/gizmo-extra/fftw-2.1.5/lib
#MPICHLIB  = -L/home/rad/gadget-3.25utils/openmpi/open64/openmpi-1.6.5/lib
GRACKLEINCL = -I/home/d150/d150/rad/gizmo-extra/grackle/include
GRACKLELIBS = -L/home/d150/d150/rad/gizmo-extra/grackle/lib -lgrackle
HDF5INCL  = #-I/home/d150/d150/rad/gizmo-extra/hdf5-1.10.1/include 
HDF5LIB   = #-L/home/d150/d150/rad/gizmo-extra/hdf5-1.10.1/lib -lhdf5 -ldl -lz
#HDF5INCL  = -I/home/d150/d150/rad/gizmo-extra/hdf5-1.8.13/include
#HDF5LIB   = -L/home/d150/d150/rad/gizmo-extra/hdf5-1.8.13/lib -lhdf5
endif

# RT/RD SYSTEMS ----------------------------------------
ifeq ($(SYSTYPE),"hipergator-intel")
#CC   =  mpiicc
#CXX  =  mpiicpc
CC   =  mpicc
CXX  =  mpicxx
FC   =  $(CC)
OPT += -DH5_USE_16_API #-DCONFIG_BFLOAT_8
GSL_INCL    = -I$(HPC_GSL_INC)
GSL_LIBS    = -L$(HPC_GSL_LIB)
FFTW_INCL   = -I$(HPC_FFTW_INC)
FFTW_LIBS   = -L$(HPC_FFTW_LIB)
HDF5LIB     = -L$(HPC_HDF5_LIB) -lhdf5
HDF5INCL    = -I$(HPC_HDF5_INC)
BLAS_LIBS   = -L$(HPC_MKL_LIB) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
MPICHLIB    =
GRACKLEINCL = -I${HPC_GRACKLE_INC}
GRACKLELIBS = -L${HPC_GRACKLE_LIB} -lgrackle
#OPTIMIZE += -O3 -mavx -axcore-avx2 -DUSE_MPI_IN_PLACE  -DNO_ISEND_IRECV_IN_DOMAIN -funroll-loops -no-prec-div
#OPTIMIZE += -O0 -mavx -axcore-avx2 -DUSE_MPI_IN_PLACE -DNO_ISEND_IRECV_IN_DOMAIN 
#OPTIMIZE += -O2 -xcore-avx2 -DNO_ISEND_IRECV_IN_DOMAIN
OPTIMIZE += -O2 -xcore-avx2 

#---------------------------
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE   += -qopenmp
endif
endif

ifeq ($(SYSTYPE),"hipergator-gnu")
HOME = /home/desika.narayanan
CC   =  mpicc
CXX  =  mpicxx
FC   =  $(CC)
OPT += -DH5_USE_16_API #-DCONFIG_BFLOAT_8

GSL_INCL    = -I$(HPC_GSL_INC)
GSL_LIBS    = -L$(HPC_GSL_LIB)
FFTW2_HOME  =   $(HOME)/fftw2
FFTW_INCL   = -I$(FFTW2_HOME)/include
FFTW_LIBS   = -L$(FFTW2_HOME)/lib
HDF5LIB     = -L$(HPC_HDF5_LIB) -lhdf5
HDF5INCL    = -I$(HPC_HDF5_INC)
#GRACKLEINCL = -I/home/desika.narayanan/local_gnu_hpg2/
#GRACKLELIBS = -L/home/desika.narayanan/local_gnu_hpg2/ -lgrackle
#HDF5LIB     = -lhdf5
MPICHLIB    = #

OPTIMIZE   += -O3 -funroll-loops
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE   += -fopenmp  # openmp required compiler flags
endif
endif

#----------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Stampede")
CC       =  mpicc
CXX      =  mpic++
FC       =  mpif90 -nofor_main
OPTIMIZE = -O3 -xhost -ipo -funroll-loops -no-prec-div -fp-model fast=2  # speed
OPTIMIZE += -g -Wall # compiler warnings
#OPTIMIZE += -parallel -openmp  # openmp (comment out this line if OPENMP not used)
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -parallel -openmp  # openmp required compiler flags
endif
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I$(TACC_MKL_INC)
MKL_LIBS = -L$(TACC_MKL_LIB) -mkl=sequential
##MKL_LIBS = -L$(TACC_MKL_LIB) -lm -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC) -DH5_USE_16_API
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
#MPICHLIB =
OPT     += -DUSE_MPI_IN_PLACE
## modules to load: 
## module load intel mvapich2 gsl hdf5 fftw2
##  -- performance is very similar with impi (intel-mpi) instead of mpavich2, 
##   if preferred use that with MPICHLIB line uncommented
## newest version of code needed for compatibility with calls in MPI-2 libraries
##
endif



#----------------------------
ifeq ($(SYSTYPE),"MacBookPro")
CC       =  mpicc
CXX      =  mpiccxx
FC       =  $(CC) #mpifort
OPTIMIZE = -O1 -funroll-loops
OPTIMIZE += -g -Wall # compiler warnings
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = #
MKL_LIBS = #
GSL_INCL = -I/usr/local/include -I$(PORTINCLUDE)
GSL_LIBS = -L/usr/local/lib -L$(PORTLIB)
FFTW_INCL= -I/usr/local/include
FFTW_LIBS= -L/usr/local/lib
HDF5INCL = -I/usr/local/include -I$(PORTINCLUDE) -DH5_USE_16_API
HDF5LIB  = -L/usr/local/lib -L$(PORTLIB) -lhdf5 -lz
MPICHLIB = #
OPT     += #
##
## PFH: this is my own laptop installation (2013 MacBook Pro running Yosemite)
## --
## I have installed GSL and HDF5 through MacPorts (once you have it installed, just use:
## sudo port install gsl
## sudo port install hdf5
## then the shortcut PORTINCLUDE/PORTLIB are just my own links to the macports installation
##  directories. in my case they are the default:
## PORTLIB=/opt/local/lib
## PORTINCLUDE=/opt/local/include
## --
## Unfortunately, FFTW is more complicated, since macports, fink, and other repository systems
## do not support direct installation of the MPI version of FFTW2, which is what GIZMO needs
## if you want to run with PMGRID or POWERSPEC enabled (if not, it should just compile without
## FFTW just fine). Be sure to install FFTW 2.1.5: get it from http://www.fftw.org/
## then unpack it, go into the unpacked directory, and configure it with:
## ./configure --enable-mpi --enable-type-prefix --enable-float
## (this set of commands is important to install the correct version)
## then "make" and finally "sudo make install"
## that should install it to its default location, /usr/local/, which is where FFTW_INCL/FFW_LIBS
## are set to point (to the respective include and lib sub-directories). check to make sure you
## have the fftw libraries correctly installed.
## --
## With this done, and the code successfully compiled, you should be able to run it with
## mpirun -np X ./GIZMO 1>gizmo.out 2>gizmo.err &
## (here "X" is the number of processes you want to use, I'm assuming youre running from the
##  same directory with the code so ./GIZMO is just in the local directory, and GIZMO is the
##  compiled file, and the 1> and 2> commands route stdin and stderr to the desired files)
##--
## If you're having trouble, I recommend the excellent guides to installing GADGET-2 at:
## http://astrobites.org/2011/04/02/installing-and-running-gadget-2/
## and
## https://gauge.wordpress.com/2009/06/16/pitp-2009-installing-gadget2/
## (by Nathan Goldbaum and Javiera Guedes, respectively) -- the installation should be
## nearly identical here
##
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Quest")
CC       =  mpicc
CXX      =  mpic++
FC       =  $(CC)
OPTIMIZE = -O2 -xhost -ipo -funroll-loops -no-prec-div -fp-model fast=2
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = -I$(MKLROOT)/include
MKL_LIBS = -L$(MKLROOT)/lib/intel64 -lm -lmkl_core -lmkl_sequential -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64
GSL_INCL = -I/projects/b1026/pascal/software/gsl/1.16/include
GSL_LIBS = -L/projects/b1026/pascal/software/gsl/1.16/lib -lgsl -lgslcblas -lm
FFTW_INCL= -I/projects/b1026/pascal/software/fftw/2.1.5-mvp/include
FFTW_LIBS= -L/projects/b1026/pascal/software/fftw/2.1.5-mvp/lib
HDF5INCL = -I/projects/b1026/pascal/software/hdf5/1.8.12/include -DH5_USE_16_API
HDF5LIB  = -L/projects/b1026/pascal/software/hdf5/1.8.12/lib -lhdf5 -lz
MPICHLIB = -lmpich
OPT     += -DUSE_MPI_IN_PLACE
#### modules to load:
#module load mpi/mvapich2-intel2013.2
#module use /projects/b1026/pascal/software/modules
#module load hdf5/1.8.12.1 gsl/1.16 fftw/2.1.5-mvp
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Titan")
CC       =  cc
CXX      =  CC
FC       =  $(CC) #ftn
OPTIMIZE = -O3 -ipo -funroll-loops -no-prec-div -fp-model fast=2 -static
OPTIMIZE += -g
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -openmp # (intel) openmp required compiler flags
FC       = $(CC)
endif
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = #
MKL_LIBS = #
GSL_INCL = -I$(GSL_DIR)/include
GSL_LIBS = -L$(GSL_DIR)/lib -lgsl -lgslcblas -lm
FFTW_INCL= -I/opt/cray/fftw/2.1.5.8/include
FFTW_LIBS= -L/opt/cray/fftw/2.1.5.8/lib
HDF5INCL = -I$(HDF5_DIR)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5_DIR)/lib -lhdf5 -lz
MPICHLIB =
OPT     += -DUSE_MPI_IN_PLACE
endif
## in your .bashrc file, include
## module swap PrgEnv-pgi PrgEnv-intel
## module load cray-hdf5-parallel fftw/2.1.5.8 gsl mercurial



#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Mira")
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
CC       = mpixlc_r # xl compilers appear to give significant speedup vs gcc
CXX      = mpixlcxx_r # _r for thread-safe versions, desired with openmp
OPTIMIZE = -openmp -qsmp=omp:noauto # -fopenmp for gcc or bgclang
else
CC       = mpixlc
CXX      = mpixlcxx
OPTIMIZE =
endif
FC       =  $(CC)
OPTIMIZE += -O3 -static -qnostrict -lm -ldl #-lpthread
OPTIMIZE += -g #-qlist -qsource -qreport -qlistopt # compiler warnings: qlist, etc produce list of opts
GMP_INCL = #
GMP_LIBS = #
MKL_INCL = #
MKL_LIBS = #
GSL_INCL = -I$(MIRA_GSL_INC)
GSL_LIBS = -lm -ldl -lpthread -L$(MIRA_GSL_LIB) -lgsl
FFTW_INCL= -I$(MIRA_FFTW2_INC)
FFTW_LIBS= -L$(MIRA_FFTW2_LIB)
HDF5INCL = -I$(MIRA_HDF5_INC) -DH5_USE_16_API -I$(MIRA_SZIP_INC) -I$(MIRA_LZIP_INC)
HDF5LIB  = -L$(MIRA_SZIP_LIB) -lszip -L$(MIRA_LZIP_LIB) -lz -L$(MIRA_HDF5_LIB) -lhdf5 -lz -lszip
MPICHLIB = #
OPT     += -DUSE_MPI_IN_PLACE -DREDUCE_TREEWALK_BRANCHING
##
## in .bashrc, need to define environmental variables:
##   export MIRA_HDF5_INC=/soft/libraries/hdf5/current/cnk-xl/current/include
##   export MIRA_HDF5_LIB=/soft/libraries/hdf5/current/cnk-xl/current/lib
##   export MIRA_GSL_INC=/soft/libraries/3rdparty/gsl/1.9/xl/include
##   export MIRA_GSL_LIB=/soft/libraries/3rdparty/gsl/1.9/xl/lib
##   export MIRA_SZIP_INC=/soft/libraries/alcf/current/xl/SZIP/include
##   export MIRA_SZIP_LIB=/soft/libraries/alcf/current/xl/SZIP/lib
##   export MIRA_LZIP_INC=/soft/libraries/alcf/current/xl/ZLIB/include
##   export MIRA_LZIP_LIB=/soft/libraries/alcf/current/xl/ZLIB/lib
##   export MIRA_FFTW2_INC=/home/phopkins/fftw/include
##   export MIRA_FFTW2_LIB=/home/phopkins/fftw/lib
##   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MIRA_LZIP_LIB:$MIRA_SZIP_LIB:$MIRA_FFTW2_LIB:$MIRA_GSL_LIB:$MIRA_HDF5_LIB:/bgsys/drivers/ppcfloor/comm/lib
##
## for HDF5,GSL,LZIP,SZIP these link to the current general-use versions of these libraries. the last command (adding these to the LD_LIBRARY_PATH)
##   is also critical, since they are not in the paths by default and it will be unable to find them even with the links above.
## for FFTW2, the pre-compiled libraries do not contain the correct mpi libraries, so you will have to compile your own. FFTW2 is installed and
##   compiled in the directory shown for me: you have to install it and link it yourself since the directory cannot be shared. for the gcc
##   compilers above, FFTW should be compiled with the following settings:
##     for xl compilers:
##       ./configure --prefix=$HOME/fftw --enable-mpi --enable-type-prefix --enable-float LDFLAGS=-L$HOME/lib CFLAGS=-I$HOME/include CC=mpixlc
##
## also in your .soft file, you want to enable:
##   for XL compilers:
##     +mpiwrapper-xl
##     +python
##     @default
## to load the mpi compilers and mpi wrappers, and MPICH libraries (python there is optional)
## xl appears to provide some improvement over gcc; xl-ndebug provides no noticeable further improvement, despite being more unsafe
endif



#----------------------------------------------------------------------------------------------
ifeq (Pleiades,$(findstring Pleiades,$(SYSTYPE)))
CC       =  icc -lmpi
CXX      =  icc -lmpi -lmpi++
FC       =  ifort -nofor_main -lmpi
ifeq ($(SYSTYPE),"Pleiades-Haswell")
OPTIMIZE = -O3 -ip -funroll-loops -no-prec-div -fp-model fast=2 -xCORE-AVX2 # Haswell cores
endif
ifeq ($(SYSTYPE),"Pleiades-SIBridge")
OPTIMIZE = -O3 -ip -funroll-loops -no-prec-div -fp-model fast=2 -xAVX # Sandy or Ivy-Bridge cores
endif
OPTIMIZE += -Wall # compiler warnings
ifeq (OPENMP,$(findstring OPENMP,$(CONFIGVARS)))
OPTIMIZE += -parallel -qopenmp
endif
GMP_INCL =
GMP_LIBS =
GSL_INCL =
GSL_LIBS =
FFTW_INCL= -I$(FFTW2_HOME)/include
FFTW_LIBS= -L$(FFTW2_HOME)/lib
HDF5INCL = -I$(HDF5)/include -DH5_USE_16_API
HDF5LIB  = -L$(HDF5)/lib -lhdf5 -lz -L/nasa/szip/2.1/lib -lsz
MPICHLIB =
OPT     += -DUSE_MPI_IN_PLACE
endif
##
## Notes:
##   1. modules to load:
##          module load comp-intel mpi-sgi/mpt hdf5/1.8.3/intel/mpt gsl python/2.7.9 szip
##   2. make sure you set the correct core-type: runs submitted to the wrong cores will not run
##   3. FFTW2: the pre-existing installation on Pleiades is incomplete and problematic.
##      you will need to install your own in your home directory. when building the library, use
##          ./configure --prefix=$HOME/fftw --enable-mpi --enable-type-prefix --enable-float
##      where "$HOME/fftw" can be renamed but is the install director (should be your home directory);
##      then you need to define the variable (here or in your bashrc file)
##          FFTW2_HOME=$HOME/fftw
##      (matching what you used for the installation) so that the code can find fftw2
##   4. in your job submission file, it is recommended for certain core types that additional settings
##      are used. for Sandy Bridge, they recommend:
##          setenv MPI_DSM_DISTRIBUTE 0
##          setenv KMP_AFFINITY disabled
##      before the actual lines submitting your job
##


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Ranger_intel")
CC       =  mpicc
CXX      =  mpiCC
FC       =  $(CC)
OPTIMIZE = -O3 -xO -ipo -funroll-loops -no-prec-div -fp-model fast=2  # speed
OPTIMIZE += -parallel -openmp  # openmp
OPTIMIZE += -g -Wall -debug parallel # compiler warnings
GMP_INCL = -I$(TACC_GMP_INC)
GMP_LIBS = -L$(TACC_GMP_LIB)
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC)
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
MPICHLIB =      # must be empty if using openmpi
OPT     += -DFIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
##
## Notes:
## 
## include the following in your .bashrc file (there is no default fftw2 module):
## module load intel/10.1 openmpi/1.2.4 gmp gsl hdf5 #now have to add fftw2 manually
## export TACC_FFTW2_INC=/opt/apps/intel10_1/openmpi_1_2_4/fftw2/2.1.5/include
## export TACC_FFTW2_LIB=/opt/apps/intel10_1/openmpi_1_2_4/fftw2/2.1.5/lib
## export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/binutils-amd/070220/lib64
##
## Options
## OPT += -DNOCALLSOFSYSTEM -DNO_ISEND_IRECV_IN_DOMAIN -DMPICH_IGNORE_CXX_SEEK
##   are not necessary, but may improve stability in some cases
##
endif


#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"Ranger_pgi")
CC       =  mpicc 
CXX      =  mpiCC
FC       =  $(CC)
OPTIMIZE = -tp barcelona-64 -fast -Mipa=fast,inline -Munroll -Mvect -O4
OPTIMIZE += -mp -Mconcur  # openmp
OPTIMIZE += -Wall  # compiler warnings
GMP_INCL = -I$(TACC_GMP_INC)
GMP_LIBS = -L$(TACC_GMP_LIB)
GSL_INCL = -I$(TACC_GSL_INC)
GSL_LIBS = -L$(TACC_GSL_LIB)
FFTW_INCL= -I$(TACC_FFTW2_INC)
FFTW_LIBS= -L$(TACC_FFTW2_LIB)
HDF5INCL = -I$(TACC_HDF5_INC)
HDF5LIB  = -L$(TACC_HDF5_LIB) -lhdf5 -lz
OPT     += -DFIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
OPT     += -DNOCALLSOFSYSTEM -DNO_ISEND_IRECV_IN_DOMAIN -DMPICH_IGNORE_CXX_SEEK
## 
## Notes:
##
## include the following in your .bashrc file:
##   module load pgi mvapich gmp gsl fftw2 hdf5
##   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/binutils-amd/070220/lib64
## 
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"odyssey")
CC       =  mpicc     # sets the C-compiler
OPT      +=  -DMPICH_IGNORE_CXX_SEEK 
FC       =  $(CC)
OPTIMIZE = -g -O2 -Wall -Wno-unused-but-set-variable
GSL_INCL =
GSL_LIBS =
FFTW_INCL=
FFTW_LIBS=
MPICHLIB =
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5 -lz
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"antares")
CC       =  mpicc     # sets the C-compiler
OPT      +=  -DMPICH_IGNORE_CXX_SEEK
FC       =  $(CC)
OPTIMIZE = -g -O2 -Wall -Wno-unused-but-set-variable
GSL_INCL = -I/home/ptorrey/local/gsl-2.1/include
GSL_LIBS = -L/home/ptorrey/local/gsl-2.1/lib -lgsl -lm
FFTW_INCL= -I/home/ptorrey/local/fftw-2.1.5/include
FFTW_LIBS= -L/home/ptorrey/local/fftw-2.1.5/libGSL_INCL
MPICHLIB =
HDF5INCL =  -DH5_USE_16_API
HDF5LIB  =  -lhdf5 -lz
endif

#----------------------------------------------------------------------------------------------
ifeq ($(SYSTYPE),"CITA")
CC       =  mpicc
CXX      =  mpicxx
OPTIMIZE =  -O3 -Wall
GSL_INCL =  -I/usr/include/gsl
GSL_LIBS =  -L/usr/lib/libgsl
FFTW_INCL=  -I/opt/fftw-2.1.5/include
FFTW_LIBS=  -L/opt/fftw-2.1.5/lib
MPICHLIB =  -L/usr/lib/libmpi
HDF5INCL =  -I/usr/include
HDF5LIB  =  -L/usr/lib/libhdf5 -static -lhdf5 -lz
endif 

#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------



ifneq (HAVE_HDF5,$(findstring HAVE_HDF5,$(CONFIGVARS)))
HDF5INCL =
HDF5LIB  =
endif


ifeq (GRACKLE,$(findstring GRACKLE,$(CONFIGVARS)))
OPT += -DCONFIG_BFLOAT_8
else
GRACKLEINCL =
GRACKLELIBS =
endif



SYSTEM_OBJS =   system/system.o system/allocate.o system/mymalloc.o system/parallel_sort.o \
                system/peano.o system/parallel_sort_special.o system/mpi_util.o

GRAVITY_OBJS  = gravity/forcetree.o gravity/cosmology.o gravity/pm_periodic.o gravity/potential.o \
                gravity/gravtree.o gravity/forcetree_update.o gravity/pm_nonperiodic.o gravity/longrange.o \
                gravity/ags_hsml.o

HYDRO_OBJS = hydro/hydra_master.o hydro/density.o hydro/gradients.o eos/eos.o


L3_OBJS =


OPTIONS = $(OPTIMIZE) $(OPT) 

FOPTIONS = $(OPTIMIZE) $(FOPT)

EXEC   = GIZMO

OBJS  =  main.o accel.o  timestep.o init.o restart.o io.o \
         predict.o global.o begrun.o run.o allvars.o read_ic.o \
         domain.o driftfac.o kicks.o ngb.o compile_time_info.o merge_split.o
FOBJS =

OBJS	+= $(GRAVITY_OBJS) $(HYDRO_OBJS) $(SYSTEM_OBJS)
OBJS	+= $(L3_OBJS)

INCL    += allvars.h proto.h gravity/forcetree.h domain.h system/myqsort.h kernel.h eos/eos.h Makefile \


ifeq (GALSF_SUBGRID_VARIABLEVELOCITY_DM_DISPERSION,$(findstring GALSF_SUBGRID_VARIABLEVELOCITY_DM_DISPERSION,$(CONFIGVARS)))
OBJS    += galaxy_sf/dm_dispersion_hsml.o
endif

ifeq (GRAIN_FLUID,$(findstring GRAIN_FLUID,$(CONFIGVARS)))
OBJS    += solids/grain_physics.o
endif

ifeq (GALSF,$(findstring GALSF,$(CONFIGVARS)))
OBJS    += galaxy_sf/sfr_eff.o
endif

ifeq (GALSF_FB_HII_HEATING,$(findstring GALSF_FB_HII_HEATING,$(CONFIGVARS)))
OBJS    += galaxy_sf/hII_heating.o
endif

ifeq (RT_CHEM_PHOTOION,$(findstring RT_CHEM_PHOTOION,$(CONFIGVARS)))
OBJS    += galaxy_sf/hII_heating.o
endif



ifeq (TWOPOINT_FUNCTION_COMPUTATION_ENABLED,$(findstring TWOPOINT_FUNCTION_COMPUTATION_ENABLED,$(CONFIGVARS)))
OBJS    += structure/twopoint.o
endif

ifeq (GALSF_FB_SNE_HEATING,$(findstring GALSF_FB_SNE_HEATING,$(CONFIGVARS)))
OBJS    += galaxy_sf/mechanical_fb.o
endif

ifeq (GALSF_FB_RPWIND_LOCAL,$(findstring GALSF_FB_RPWIND_LOCAL,$(CONFIGVARS)))
OBJS    += galaxy_sf/rp_localwinds.o
endif

ifeq (BLACK_HOLES,$(findstring BLACK_HOLES,$(CONFIGVARS)))
OBJS    += galaxy_sf/blackholes/blackhole.o
OBJS    += galaxy_sf/blackholes/blackhole_util.o
OBJS    += galaxy_sf/blackholes/blackhole_environment.o
OBJS    += galaxy_sf/blackholes/blackhole_feed.o
OBJS    += galaxy_sf/blackholes/blackhole_swallow_and_kick.o
INCL    += galaxy_sf/blackholes/blackhole.h
endif


ifeq (SINGLE_STAR,$(findstring SINGLE_STAR,$(CONFIGVARS)))
OBJS	+= radiation/rt_utilities.o radiation/rt_CGmethod.o radiation/rt_source_injection.o radiation/rt_chem.o radiation/rt_cooling.o
OBJS    += galaxy_sf/sfr_eff.o galaxy_sf/hII_heating.o galaxy_sf/mechanical_fb.o galaxy_sf/rp_localwinds.o
OBJS    += galaxy_sf/blackholes/blackhole.o galaxy_sf/blackholes/blackhole_util.o galaxy_sf/blackholes/blackhole_environment.o galaxy_sf/blackholes/blackhole_feed.o galaxy_sf/blackholes/blackhole_swallow_and_kick.o
INCL    += galaxy_sf/blackholes/blackhole.h
endif



ifeq (SCFPOTENTIAL,$(findstring SCFPOTENTIAL,$(CONFIGVARS)))
OBJS    += modules/potentials/scf.o modules/potentials/scf_util.o
endif

ifeq (FOF,$(findstring FOF,$(CONFIGVARS)))
OBJS    += structure/fof.o
INCL	+= structure/fof.h
endif

ifeq (OUTPUTLINEOFSIGHT,$(findstring OUTPUTLINEOFSIGHT,$(CONFIGVARS)))
OBJS    += structure/lineofsight.o
endif

ifeq (COOLING,$(findstring COOLING,$(CONFIGVARS)))
OBJS    += cooling/cooling.o
INCL	+= cooling/cooling.h
endif

ifeq (GRACKLE,$(findstring GRACKLE,$(CONFIGVARS)))
OBJS    += cooling/grackle3.o
endif

ifeq (BUBBLES,$(findstring BUBBLES,$(CONFIGVARS)))
OBJS    += modules/bubbles/bubbles.o
endif

ifeq (EOS_HELMHOLTZ,$(findstring EOS_HELMHOLTZ,$(CONFIGVARS)))
OBJS    += eos/eos_interface.o
INCL    += eos/helmholtz/helm_wrap.h
FOBJS   += eos/helmholtz/helm_impl.o eos/helmholtz/helm_wrap.o
FINCL   += eos/helmholtz/helm_const.dek eos/helmholtz/helm_implno.dek eos/helmholtz/helm_table_storage.dek eos/helmholtz/helm_vector_eos.dek
endif

ifeq (IMPOSE_PINNING,$(findstring IMPOSE_PINNING,$(CONFIGVARS)))
OBJS	+= system/pinning.o
endif

ifeq (DISTORTIONTENSORPS,$(findstring DISTORTIONTENSORPS,$(CONFIGVARS)))
OBJS	+= modules/phasespace/phasespace.o modules/phasespace/phasespace_math.o
endif

ifeq (RT_,$(findstring RT_,$(CONFIGVARS)))
OBJS	+= radiation/rt_utilities.o radiation/rt_CGmethod.o radiation/rt_source_injection.o radiation/rt_chem.o radiation/rt_cooling.o
endif

ifeq (SUBFIND,$(findstring SUBFIND,$(CONFIGVARS)))
OBJS	+= subfind/subfind.o subfind/subfind_vars.o subfind/subfind_collective.o subfind/subfind_serial.o subfind/subfind_so.o subfind/subfind_cont.o \
	subfind/subfind_distribute.o subfind/subfind_findlinkngb.o subfind/subfind_nearesttwo.o subfind/subfind_loctree.o subfind/subfind_alternative_collective.o subfind/subfind_reshuffle.o \
	subfind/subfind_potential.o subfind/subfind_density.o
INCL	+= subfind/subfind.h
endif

ifeq (SIDM,$(findstring SIDM,$(CONFIGVARS)))
OBJS    +=  sidm/sidm_core.o sidm/sidm_allvars.o
INCL    +=  sidm/sidm_proto.h
endif

ifeq (NUCLEAR_NETWORK,$(findstring NUCLEAR_NETWORK,$(CONFIGVARS)))
OBJS	+=  nuclear/nuclear_network_solver.o nuclear/nuclear_network.o
INCL	+=  nuclear/nuclear_network.h
endif

ifeq (TURB_DRIVING,$(findstring TURB_DRIVING,$(CONFIGVARS)))
OBJS	+= turb/turb_driving.o turb/turb_powerspectra.o
endif

######################
## ADDED BY ROMEEL ###
######################
ifeq (FOFRAD,$(findstring FOFRAD,$(CONFIGVARS)))
OBJS    += structure/fofrad.o
INCL    += structure/fofrad.h
endif

ifeq (GALSF_INSTANTANEOUS_METALS,$(findstring GALSF_INSTANTANEOUS_METALS,$(CONFIGVARS)))
OBJS    += galaxy_sf/metals.o
endif

ifeq (HOTGAS_QUENCH,$(findstring HOTGAS_QUENCH,$(CONFIGVARS)))
OBJS    += galaxy_sf/hotgas_quench.o
endif

######################
## ADDED BY QI LI ####
######################
ifeq (GALSF_DUST,$(findstring GALSF_DUST,$(CONFIGVARS)))
OBJS    += dust/dust.o
INCL    += dust/dust.h
endif
######################

CFLAGS = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL) $(GMP_INCL) $(GRACKLEINCL)

ifeq (VIP,$(findstring VIP,$(CONFIGVARS)))
FFLAGS = $(FOPTIONS)
else
FFLAGS = $(OPTIONS)
endif


ifeq (ALTERNATIVE_PSORT,$(findstring ALTERNATIVE_PSORT,$(CONFIGVARS)))
OBJS  += fof_alt_psort.o modules/psort-1.0/error_handling.o
CXXFLAGS = $(CFLAGS)
FC    = $(CXX)
endif

FFTW = $(FFTW_LIBS)  $(FFTW_LIBNAMES) 


LIBS   = -lm $(HDF5LIB) -g $(GSL_LIBS) -lgsl $(FFTW) $(GRACKLELIBS) $(BLAS_LIBS)

ifeq (PTHREADS_NUM_THREADS,$(findstring PTHREADS_NUM_THREADS,$(CONFIGVARS))) 
LIBS   +=  -lpthread
endif

$(EXEC): $(OBJS) $(FOBJS)  
	$(FC) $(OPTIMIZE) $(OBJS) $(FOBJS) $(LIBS) $(RLIBS) -o $(EXEC)

$(OBJS): $(INCL)  $(CONFIG)  compile_time_info.c

$(FOBJS): %.o: %.f90
	$(FC) $(OPTIMIZE) -c $< -o $@

compile_time_info.c: $(CONFIG)
	$(PERL) prepare-config.perl $(CONFIG)

clean:
	rm -f $(OBJS) $(FOBJS) $(EXEC) *.oo *.c~ compile_time_info.c GIZMO_config.h


