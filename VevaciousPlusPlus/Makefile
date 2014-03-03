#
# Makefile for VevaciousPlusPlus.
#     Author: Ben O'Leary (benjamin.oleary@gmail.com)
#      Copyright 2014 Ben O'Leary
#
#      This file is part of VevaciousPlusPlus, released under the GNU General
#      Public License. Please see the accompanying README.VevaciousPlusPlus.txt
#      file for a full list of files, brief documentation on how to use these
#      classes, and further details on the license.
#

SHELL=/bin/sh

CC=g++
AR=ar
ARFLAGS=rcs
LHPCDIR=./LHPC/
LHPCLIBNAME=LHPC
LHPCINCLUDEDIR=$(LHPCDIR)/include/
LHPCLIBDIR=$(LHPCDIR)/lib/
EIGENDIR=./Eigen/
CCFLAGS=-Wall -O3 -I./include/ -I$(LHPCINCLUDEDIR) -I$(EIGENDIR)
CCOBJECTFLAGS=-c
CCEXEFLAGS=
SOURCES=source/PotentialEvaluation/HomotopyContinuationReadyPotential.cpp \
 source/PotentialEvaluation/MassSquaredMatrix.cpp \
 source/PotentialEvaluation/ParameterFunctionoid.cpp \
 source/PotentialEvaluation/PolynomialCoefficient.cpp \
 source/PotentialEvaluation/PolynomialSum.cpp \
 source/PotentialEvaluation/PolynomialTerm.cpp \
 source/PotentialEvaluation/PotentialFromPolynomialAndMasses.cpp \
 source/PotentialEvaluation/PotentialFunction.cpp \
 source/PotentialEvaluation/RunningParameterManager.cpp \
 source/PotentialEvaluation/SlhaFunctionoid.cpp \
 source/PotentialMinimization/HomotopyContinuationAndGradient.cpp \
 source/PotentialMinimization/PotentialMinimizer.cpp \
 source/PotentialMinimization/PotentialMinimum.cpp \
 source/TunnelingCalculation/BounceWithSplines.cpp \
 source/TunnelingCalculation/TunnelingCalculator.cpp \
 source/VevaciousPlusPlus.cpp
OBJECTS=$(SOURCES:.cpp=.o)
BINDIR=./bin/
LIBNAME=VevaciousPlusPlus
LIBDIR=./lib/
LDFLAGS= -L$(LIBDIR) -l$(LIBNAME)
LIBRARYFILE=$(LIBDIR)/lib$(LIBNAME).a
LDFLAGS= -L$(LIBDIR) -l$(LIBNAME) -L$(LHPCLIBDIR) -l$(LHPCLIBNAME)
MAINCPP=source/VevaciousPlusPlusMain.cpp
EXEFILE=$(BINDIR)/VevaciousPlusPlus.exe

all: $(OBJECTS) $(LIBRARYFILE) $(EXEFILE)

.PHONY: all clean

$(OBJECTS): %.o: %.cpp
	$(CC) $(CCFLAGS) $(CCOBJECTFLAGS) $< -o $@
	
$(LIBDIR):
	mkdir $@

$(LIBRARYFILE): $(LIBDIR) $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)

$(BINDIR):
	mkdir $@
	
$(EXEFILE): $(BINDIR) $(LIBRARYFILE)
	$(CC) $(CCFLAGS) $(CCEXEFLAGS) -o $@ $(MAINCPP) $(LDFLAGS)

clean:
	\rm -f $(OBJECTS) $(LIBRARYFILE) $(EXEFILE)

deleteobjects:
	\rm -f $(OBJECTS)