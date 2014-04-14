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
#LHPCDIR=./LHPC/
LHPCDIR=../../../machine-specific/LHPC/
LHPCLIBNAME=LHPC
LHPCINCLUDEDIR=$(LHPCDIR)/include/
LHPCLIBDIR=$(LHPCDIR)/lib/
#EIGENDIR=./Eigen/
EIGENDIR=../../../machine-specific/Eigen/
#MINUITLIBNAME=lcg_Minuit
MINUITLIBNAME=Minuit2
#MINUITINCLUDEDIR=./MINUIT/include/
MINUITINCLUDEDIR=../../../machine-specific/Minuit2/include/
#MINUITLIBDIR=./MINUIT/include/
MINUITLIBDIR=../../../machine-specific/Minuit2/lib/
DEBUGFLAGS=
INCLUDES=-I./include/ -I$(LHPCINCLUDEDIR) -I$(EIGENDIR) -I$(MINUITINCLUDEDIR)
CCFLAGS=-Wall -O3 -fopenmp $(DEBUGFLAGS) $(INCLUDES)
CCOBJECTFLAGS=-c
CCEXEFLAGS=
PARAMETERFUNCTIONOIDS=source/PotentialEvaluation/ParameterFunctionoids/
SOURCES=$(PARAMETERFUNCTIONOIDS)/BinaryOperationFunctionoid.cpp \
 $(PARAMETERFUNCTIONOIDS)/ConstantFunctionoid.cpp \
 $(PARAMETERFUNCTIONOIDS)/SlhaFunctionoid.cpp \
 $(PARAMETERFUNCTIONOIDS)/UnaryOperationFunctionoid.cpp \
 $(PARAMETERFUNCTIONOIDS)/ParameterFunctionoid.cpp \
 source/PotentialEvaluation/ComplexMassSquaredMatrix.cpp \
 source/PotentialEvaluation/FixedScaleOneLoopPotential.cpp \
 source/PotentialEvaluation/HomotopyContinuationReadyPolynomial.cpp \
 source/PotentialEvaluation/HomotopyContinuationReadyPotential.cpp \
 source/PotentialEvaluation/MassesSquaredCalculator.cpp \
 source/PotentialEvaluation/PolynomialSum.cpp \
 source/PotentialEvaluation/PolynomialTerm.cpp \
 source/PotentialEvaluation/PotentialFromPolynomialAndMasses.cpp \
 source/PotentialEvaluation/PotentialFunction.cpp \
 source/PotentialEvaluation/ProductOfPolynomialSums.cpp \
 source/PotentialEvaluation/RealMassesSquaredMatrix.cpp \
 source/PotentialEvaluation/RgeImprovedOneLoopPotential.cpp \
 source/PotentialEvaluation/RunningParameterManager.cpp \
 source/PotentialEvaluation/SlhaManager.cpp \
 source/PotentialEvaluation/SymmetricComplexMassMatrix.cpp \
 source/PotentialEvaluation/ThermalFunctions.cpp \
 source/PotentialMinimization/Hom4ps2AndMinuit.cpp \
 source/PotentialMinimization/HomotopyContinuationAndGradient.cpp \
 source/PotentialMinimization/MinuitManager.cpp \
 source/PotentialMinimization/MinuitMinimum.cpp \
 source/PotentialMinimization/PotentialForMinuit.cpp \
 source/PotentialMinimization/PotentialMinimizer.cpp \
 source/PotentialMinimization/PotentialMinimum.cpp \
 source/TunnelingCalculation/BounceWithSplines.cpp \
 source/TunnelingCalculation/TunnelingCalculator.cpp \
 source/VevaciousPlusPlus.cpp
OBJECTS=$(SOURCES:.cpp=.o)
BINDIR=./bin/
LIBNAME=VevaciousPlusPlus
LIBDIR=./lib/
LIBRARYFILE=$(LIBDIR)/lib$(LIBNAME).a
LDFLAGS= -L$(LIBDIR) -l$(LIBNAME) \
 -L$(LHPCLIBDIR) -l$(LHPCLIBNAME) \
  $(MINUITLIBDIR)/lib$(MINUITLIBNAME).a
# -l$(LHPCLIBNAME) -L$(MINUITLIBDIR) -l:lib$(MINUITLIBNAME).a
# There is some silly screwing around to force linking to the static Minuit2 
# library.
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
