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
#BOOSTDIR=./Boost/
BOOSTDIR=../../../machine-specific/Boost/
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
INCLUDES=-I./include/ \
 -I$(LHPCINCLUDEDIR) -I$(EIGENDIR) -I$(BOOSTDIR) -I$(MINUITINCLUDEDIR)
CCFLAGS=-Wall -O3 -fopenmp $(DEBUGFLAGS) $(INCLUDES)
CCOBJECTFLAGS=-c
CCEXEFLAGS=
BOUNCEACTIONPATHFINDING= source/BounceActionEvaluation/BounceActionPathFinding/
PATHPARAMETERIZATION= source/BounceActionEvaluation/PathParameterization/
MASSESSQUAREDCALCULATORS= source/PotentialEvaluation/MassesSquaredCalculators/
PARAMETERFUNCTIONOIDS= source/PotentialEvaluation/ParameterFunctionoids/
POTENTIALFUNCTIONS= source/PotentialEvaluation/PotentialFunctions/
POTENTIALEVALUATION= source/PotentialMinimization/
HOMOTOPYCONTINUATION= $(POTENTIALEVALUATION)/HomotopyContinuation/
GRADIENTBASEDMINIMIZATION= $(POTENTIALEVALUATION)/GradientBasedMinimization/
BOUNCEACTIONTUNNELING= source/TunnelingCalculation/BounceActionTunneling/
SOURCES= source/BasicFunctions/PolynomialSum.cpp \
 source/BasicFunctions/PolynomialTerm.cpp \
 source/BasicFunctions/ProductOfPolynomialSums.cpp \
 source/BasicFunctions/SimplePolynomial.cpp \
 source/BasicFunctions/SumOfProductOfPolynomialSums.cpp \
 source/BounceActionEvaluation/BounceActionCalculator.cpp \
 $(BOUNCEACTIONPATHFINDING)/FullPathVaryingMinuit.cpp \
 $(BOUNCEACTIONPATHFINDING)/MinuitNodePotentialMinimizer.cpp \
 $(BOUNCEACTIONPATHFINDING)/MinuitPathBounceMinimizer.cpp \
 $(BOUNCEACTIONPATHFINDING)/MinuitPathFinder.cpp \
 $(BOUNCEACTIONPATHFINDING)/MinuitPathPotentialMinimizer.cpp \
 $(BOUNCEACTIONPATHFINDING)/SingleNodeVaryingMinuit.cpp \
 source/BounceActionEvaluation/BubbleProfile.cpp \
 source/BounceActionEvaluation/BubbleRadialValueDescription.cpp \
 source/BounceActionEvaluation/BubbleShootingOnSpline.cpp \
 source/BounceActionEvaluation/OdeintBubbleDerivatives.cpp \
 source/BounceActionEvaluation/OdeintBubbleObserver.cpp \
 $(PATHPARAMETERIZATION)/LinearSplinePathSegment.cpp \
 $(PATHPARAMETERIZATION)/LinearSplineThroughNodes.cpp \
 $(PATHPARAMETERIZATION)/LinearSplineThroughNodesFactory.cpp \
 $(PATHPARAMETERIZATION)/NodesFromParameterization.cpp \
 $(PATHPARAMETERIZATION)/NodesOnBisectingPlanes.cpp \
 $(PATHPARAMETERIZATION)/NodesOnParallelPlanes.cpp \
 $(PATHPARAMETERIZATION)/NodesOnPlanes.cpp \
 $(PATHPARAMETERIZATION)/PathFromNodesFactory.cpp \
 $(PATHPARAMETERIZATION)/PolynomialThroughNodes.cpp \
 $(PATHPARAMETERIZATION)/PolynomialThroughNodesFactory.cpp \
 $(PATHPARAMETERIZATION)/QuadraticSplinePathSegment.cpp \
 $(PATHPARAMETERIZATION)/QuadraticSplineThroughNodes.cpp \
 $(PATHPARAMETERIZATION)/QuadraticSplineThroughNodesFactory.cpp \
 $(PATHPARAMETERIZATION)/TunnelPath.cpp \
 $(PATHPARAMETERIZATION)/TunnelPathFactory.cpp \
 source/BounceActionEvaluation/SplinePotential.cpp \
 source/PotentialEvaluation/MassesSquaredCalculator.cpp \
 $(MASSESSQUAREDCALCULATORS)/ComplexMassSquaredMatrix.cpp \
 $(MASSESSQUAREDCALCULATORS)/RealMassesSquaredMatrix.cpp \
 $(MASSESSQUAREDCALCULATORS)/SymmetricComplexMassMatrix.cpp \
 source/PotentialEvaluation/ParameterFunctionoid.cpp \
 $(PARAMETERFUNCTIONOIDS)/BinaryOperationFunctionoid.cpp \
 $(PARAMETERFUNCTIONOIDS)/ConstantFunctionoid.cpp \
 $(PARAMETERFUNCTIONOIDS)/SlhaFunctionoid.cpp \
 $(PARAMETERFUNCTIONOIDS)/UnaryOperationFunctionoid.cpp \
 source/PotentialEvaluation/PotentialFunction.cpp \
 $(POTENTIALFUNCTIONS)/FixedScaleOneLoopPotential.cpp \
 $(POTENTIALFUNCTIONS)/PotentialFromPolynomialAndMasses.cpp \
 $(POTENTIALFUNCTIONS)/RgeImprovedOneLoopPotential.cpp \
 source/PotentialEvaluation/SlhaManagement/RunningParameterManager.cpp \
 source/PotentialEvaluation/SlhaManagement/SlhaManager.cpp \
 source/PotentialEvaluation/SlhaManagement/SlhaUpdatePropagator.cpp \
 source/PotentialEvaluation/ThermalFunctions.cpp \
 $(HOMOTOPYCONTINUATION)/BasicPolynomialHomotopyContinuation.cpp \
 $(HOMOTOPYCONTINUATION)/FieldPolynomialsWithScale.cpp \
 $(HOMOTOPYCONTINUATION)/Hom4ps2Runner.cpp \
 $(HOMOTOPYCONTINUATION)/HomotopyContinuationSolver.cpp \
 $(HOMOTOPYCONTINUATION)/HomotopyContinuationTargetSystem.cpp \
 $(HOMOTOPYCONTINUATION)/PolynomialGradientTargetSystem.cpp \
 $(GRADIENTBASEDMINIMIZATION)/MinuitMinimum.cpp \
 $(GRADIENTBASEDMINIMIZATION)/MinuitPotentialMinimizer.cpp \
 $(GRADIENTBASEDMINIMIZATION)/PotentialForMinuit.cpp \
 source/PotentialMinimization/GradientFromStartingPoints.cpp \
 source/PotentialMinimization/GradientMinimizer.cpp \
 source/PotentialMinimization/PotentialMinimizer.cpp \
 source/PotentialMinimization/PotentialMinimum.cpp \
 source/PotentialMinimization/StartingPointFinder.cpp \
 source/TunnelingCalculation/BounceActionTunneler.cpp \
 $(BOUNCEACTIONTUNNELING)/BounceActionForMinuit.cpp \
 $(BOUNCEACTIONTUNNELING)/BounceAlongPathWithThreshold.cpp \
 $(BOUNCEACTIONTUNNELING)/BouncePathFinder.cpp \
 $(BOUNCEACTIONTUNNELING)/CosmoTransitionsRunner.cpp \
 $(BOUNCEACTIONTUNNELING)/MinuitBounceActionMinimizer.cpp \
 $(BOUNCEACTIONTUNNELING)/ModifiedBounceForMinuit.cpp \
 $(BOUNCEACTIONTUNNELING)/PathFieldsAndPotential.cpp \
 $(BOUNCEACTIONTUNNELING)/PathFromNodes.cpp \
 $(BOUNCEACTIONTUNNELING)/ThermalActionFitter.cpp \
 source/TunnelingCalculation/TunnelingCalculator.cpp \
 source/VersionInformation.cpp \
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
