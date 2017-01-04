
# Important: Specify the full pathname here; can't do ~/Desktop ...
# If you get the error message that a module *.mod can't be found, it's because you didn't put the precise
# path name in the line below

# *** Mac ***
export AT="/Users/michaelporter/Documents/at"

# *** Linux ***
# export AT="/home/porter/at"

# *** Windows ***
# make sure you change -O3 to -O2 below in the gfortran options also
# export AT="C:\Documents and Settings\Michael Porter\My Documents\at"

# Here's the sequence to install the Acoustics Toolbox:
# make clean
# make install

# -openmp seems to slow things down a lot ... probably memory contention

# ______________________________________________________________________________

# *** ifort
# need -heap-arrays to avoid stack overflows for big runs ...
export FC=ifort
export FFLAGS= -O3 -parallel -axSSE4.2 -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I$(AT)/misc
# export FFLAGS= -O3 -parallel -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I$(AT)/misc
export FFLAGS= -O3 -openmp -openmp_report -parallel -axSSE4.2 -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I$(AT)/misc

# -parallel off for export:
# export FFLAGS= -O3 -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I$(AT)/misc
# diagnostics on:
# export FFLAGS= -O2 -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I$(AT)/misc -check -traceback
# profiling:
# export FFLAGS= -g -nologo -inline-level=2 -assume byterecl -threads -heap-arrays -I$(AT)/misc -check -traceback

# ______________________________________________________________________________

# *** gfortran
# -O2 was the highest level of optimization that worked under Windows
# -fopenmp slows things down
export FC=gfortran
# export FFLAGS= -march=native -Wall -std=gnu -O2 -I$(AT)/misc
# export FFLAGS= -fopenmp -march=native -Wall -std=gnu -O3 -I$(AT)/misc
# export FFLAGS= -fopenmp -march=native -Wall -std=gnu -O3 -I$(AT)/misc -ffast-math -funroll-all-loops -msse3 -fomit-frame-pointer -mtune=native -Q --help=optimizers
export FFLAGS= -march=native -Waliasing -Wampersand -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -std=gnu -O3 -ffast-math -funroll-all-loops -fomit-frame-pointer -mtune=native -I$(AT)/misc

# diagnostics on:
# export FFLAGS= -march=native -Wall -std=gnu -O -fcheck=all -I$(AT)/misc
# profiling:
# export FFLAGS= -gp -Wall -std=f95 -fbounds-check -I$(AT)/misc

# ______________________________________________________________________________

# *** g95
# export FC=g95
# export FFLAGS= -Wall -std=f2003 -O3 -I$(AT)/misc
# diagnostics on:
# export FFLAGS= -Wall -std=f2003 -ftrace=full -fbounds-check -I$(AT)/misc

# ______________________________________________________________________________

export RM=rm
export CC=gcc
export CFLAGS=-g

# KRAKEL is commented out below because it requires the LAPACK library
# If you have LAPACK, you can uncomment it
install:
	(cd misc;	make -k all)
	(cd tslib;	make -k install)
	(cd Kraken;	make -k install)
	#(cd Krakel;	make -k install)
	(cd Scooter;	make -k install)
	(cd Bellhop;	make -k install)
	@echo " "
	@echo "*************************************"
	@echo "****** Acoustics Toolbox installed **"
	@echo "*************************************"

clean:
	(rm -f */*.o lib/*.a module/*.mod)
	(cd misc;	make -k -i clean)
	(cd tslib;	make -k -i clean)
	(cd Kraken;	make -k -i clean)
	(cd Krakel;	make -k -i clean)
	(cd Scooter;	make -k -i clean)
	(cd Bellhop;	make -k -i clean)
	(cd tests;	make -k -i clean)
