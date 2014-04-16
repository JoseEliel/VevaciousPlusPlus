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
MASSESSQUAREDCALCULATORS=source/PotentialEvaluation/MassesSquaredCalculators/
PARAMETERFUNCTIONOIDS=source/PotentialEvaluation/ParameterFunctionoids/
POTENTIALFUNCTIONS=source/PotentialEvaluation/PotentialFunctions/
HOMOTOPYCONTINUATION=source/PotentialMinimization/HomotopyContinuation/
MINUITMINIMIZATION=source/PotentialMinimization/MinuitMinimization/
SOURCES=$(MASSESSQUAREDCALCULATORS)/ComplexMassSquaredMatrix.cpp \
 $(MASSESSQUAREDCALCULATORS)/MassesSquaredCalculator.cpp \
 $(MASSESSQUAREDCALCULATORS)/RealMassesSquaredMatrix.cpp \
 $(MASSESSQUAREDCALCULATORS)/SymmetricComplexMassMatrix.cpp \
 $(PARAMETERFUNCTIONOIDS)/BinaryOperationFunctionoid.cpp \
 $(PARAMETERFUNCTIONOIDS)/ConstantFunctionoid.cpp \
 $(PARAMETERFUNCTIONOIDS)/SlhaFunctionoid.cpp \
 $(PARAMETERFUNCTIONOIDS)/UnaryOperationFunctionoid.cpp \
 $(PARAMETERFUNCTIONOIDS)/ParameterFunctionoid.cpp \
 $(POTENTIALFUNCTIONS)/FixedScaleOneLoopPotential.cpp \
 $(POTENTIALFUNCTIONS)/PotentialFromPolynomialAndMasses.cpp \
 $(POTENTIALFUNCTIONS)/PotentialFunction.cpp \
 $(POTENTIALFUNCTIONS)/RgeImprovedOneLoopPotential.cpp \
 source/PotentialEvaluation/PolynomialSum.cpp \
 source/PotentialEvaluation/PolynomialTerm.cpp \
 source/PotentialEvaluation/ProductOfPolynomialSums.cpp \
 source/PotentialEvaluation/RunningParameterManager.cpp \
 source/PotentialEvaluation/SumOfProductOfPolynomialSums.cpp \
 source/PotentialEvaluation/SlhaManager.cpp \
 source/PotentialEvaluation/ThermalFunctions.cpp \
 $(HOMOTOPYCONTINUATION)/BasicPolynomialHomotopyContinuation.cpp \
 $(HOMOTOPYCONTINUATION)/Hom4ps2Runner.cpp \
 $(HOMOTOPYCONTINUATION)/HomotopyContinuationSolver.cpp \
 $(HOMOTOPYCONTINUATION)/HomotopyContinuationTargetSystem.cpp \
 $(HOMOTOPYCONTINUATION)/PolynomialGradientTargetSystem.cpp \
 $(MINUITMINIMIZATION)/MinuitManager.cpp \
 $(MINUITMINIMIZATION)/MinuitMinimum.cpp \
 $(MINUITMINIMIZATION)/PotentialForMinuit.cpp \
 source/PotentialMinimization/HomotopyContinuationAndGradient.cpp \
 source/PotentialMinimization/HomotopyContinuationAndMinuit.cpp \
 source/PotentialMinimization/PotentialMinimizer.cpp \
 source/PotentialMinimization/PotentialMinimum.cpp \
 source/TunnelingCalculation/BounceWithSplines.cpp \
 source/TunnelingCalculation/CosmoTransitionsRunner.cpp \
 source/TunnelingCalculation/MinuitBounceActionMinimizer.cpp \
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
