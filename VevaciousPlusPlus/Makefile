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
BOOSTDIR=/home/bol/BOL/ProjectDependencies/Boost/
#LHPCDIR=./LHPC/
LHPCDIR=/home/bol/BOL/ProjectDependencies/LHPC/
LHPCLIBNAME=LHPC
LHPCINCLUDEDIR=$(LHPCDIR)/include/
LHPCLIBDIR=$(LHPCDIR)/lib/
#EIGENDIR=./Eigen/
EIGENDIR=/home/bol/BOL/ProjectDependencies/Eigen/
#MINUITLIBNAME=lcg_Minuit
MINUITLIBNAME=Minuit2
#MINUITINCLUDEDIR=./MINUIT/include/
MINUITINCLUDEDIR=/home/bol/BOL/ProjectDependencies/Minuit2/include/
#MINUITLIBDIR=./MINUIT/include/
MINUITLIBDIR=/home/bol/BOL/ProjectDependencies/Minuit2/lib/
DEBUGFLAGS=
INCLUDES=-I./include/ \
 -I$(LHPCINCLUDEDIR) -I$(EIGENDIR) -I$(BOOSTDIR) -I$(MINUITINCLUDEDIR)
CCFLAGS=-Wall -Wno-unused-local-typedefs -O3 -fopenmp $(DEBUGFLAGS) $(INCLUDES)
CCOBJECTFLAGS=-c
CCEXEFLAGS=

LAGRANGIANPARAMETERMANAGEMENT= source/LagrangianParameterManagement/
BUILDINGBLOCKS= source/PotentialEvaluation/BuildingBlocks/
BOUNCEACTIONPATHFINDING= source/BounceActionEvaluation/BounceActionPathFinding/
PATHPARAMETERIZATION= source/BounceActionEvaluation/PathParameterization/
MASSESSQUAREDCALCULATORS= source/PotentialEvaluation/MassesSquaredCalculators/
PARAMETERFUNCTIONOIDS= source/PotentialEvaluation/ParameterFunctionoids/
POTENTIALFUNCTIONS= source/PotentialEvaluation/PotentialFunctions/
SLHADERIVEDFUNCTIONOIDS= $(LAGRANGIANPARAMETERMANAGEMENT)/SlhaDerivedFunctionoids/
POTENTIALMINIMIZATION= source/PotentialMinimization/
HOMOTOPYCONTINUATION= $(POTENTIALMINIMIZATION)/HomotopyContinuation/
STARTINGPOINTGENERATION= $(POTENTIALMINIMIZATION)/StartingPointGeneration/
GRADIENTBASEDMINIMIZATION= $(POTENTIALMINIMIZATION)/GradientBasedMinimization/
BOUNCEACTIONTUNNELING= source/TunnelingCalculation/BounceActionTunneling/
NEWSOURCES=  $(LAGRANGIANPARAMETERMANAGEMENT)/LagrangianParameterManager.cpp \
 $(LAGRANGIANPARAMETERMANAGEMENT)/LhaLinearlyInterpolatedBlockEntry.cpp \
 $(LAGRANGIANPARAMETERMANAGEMENT)/LhaPolynomialFitBlockEntry.cpp \
 $(LAGRANGIANPARAMETERMANAGEMENT)/LesHouchesAccordBlockEntryManager.cpp \
 $(LAGRANGIANPARAMETERMANAGEMENT)/SlhaBlocksWithSpecialCasesManager.cpp \
 $(LAGRANGIANPARAMETERMANAGEMENT)/SlhaCompatibleWithSarahManager.cpp \
 $(BUILDINGBLOCKS)/ParametersAndFieldsProductTerm.cpp \
 $(BUILDINGBLOCKS)/ParametersAndFieldsProductSum.cpp

SOURCES= $(NEWSOURCES) \
 source/BasicFunctions/PolynomialSum.cpp \
 source/BasicFunctions/PolynomialTerm.cpp \
 source/BasicFunctions/ProductOfPolynomialSums.cpp \
 source/BasicFunctions/SumOfProductOfPolynomialSums.cpp \
 $(BOUNCEACTIONPATHFINDING)/MinuitOnHypersurfaces.cpp \
 $(BOUNCEACTIONPATHFINDING)/MinuitOnPotentialOnParallelPlanes.cpp \
 $(BOUNCEACTIONPATHFINDING)/MinuitOnPotentialPerpendicularToPath.cpp \
 $(BOUNCEACTIONPATHFINDING)/MinuitPathFinder.cpp \
 $(PATHPARAMETERIZATION)/LinearSplinePathSegment.cpp \
 $(PATHPARAMETERIZATION)/LinearSplineThroughNodes.cpp \
 $(PATHPARAMETERIZATION)/TunnelPath.cpp \
 source/BounceActionEvaluation/BounceActionCalculator.cpp \
 source/BounceActionEvaluation/BouncePathFinder.cpp \
 source/BounceActionEvaluation/BubbleProfile.cpp \
 source/BounceActionEvaluation/BubbleRadialValueDescription.cpp \
 source/BounceActionEvaluation/BubbleShootingOnPathInFieldSpace.cpp \
 source/BounceActionEvaluation/OdeintBubbleDerivatives.cpp \
 source/BounceActionEvaluation/OdeintBubbleObserver.cpp \
 source/BounceActionEvaluation/OneDimensionalPotentialAlongPath.cpp \
 source/BounceActionEvaluation/SplinePotential.cpp \
 source/BounceActionEvaluation/UndershootOvershootBubble.cpp \
 source/MinuitWrappersAndHelpers/MinuitMinimum.cpp \
 $(MASSESSQUAREDCALCULATORS)/BaseComplexMassMatrix.cpp \
 $(MASSESSQUAREDCALCULATORS)/ComplexMassSquaredMatrix.cpp \
 $(MASSESSQUAREDCALCULATORS)/RealMassesSquaredMatrix.cpp \
 $(MASSESSQUAREDCALCULATORS)/SymmetricComplexMassMatrix.cpp \
 $(MASSESSQUAREDCALCULATORS)/OldComplexMassSquaredMatrix.cpp \
 $(MASSESSQUAREDCALCULATORS)/OldRealMassesSquaredMatrix.cpp \
 $(MASSESSQUAREDCALCULATORS)/OldSymmetricComplexMassMatrix.cpp \
 $(PARAMETERFUNCTIONOIDS)/BinaryOperationFunctionoid.cpp \
 $(PARAMETERFUNCTIONOIDS)/ConstantFunctionoid.cpp \
 $(PARAMETERFUNCTIONOIDS)/SlhaFunctionoid.cpp \
 $(PARAMETERFUNCTIONOIDS)/UnaryOperationFunctionoid.cpp \
 $(POTENTIALFUNCTIONS)/FixedScaleOneLoopPotential.cpp \
 $(POTENTIALFUNCTIONS)/OldFixedScaleOneLoopPotential.cpp \
 $(POTENTIALFUNCTIONS)/PotentialFromPolynomialWithMasses.cpp \
 $(POTENTIALFUNCTIONS)/OldPotentialFromPolynomialAndMasses.cpp \
 $(POTENTIALFUNCTIONS)/RgeImprovedOneLoopPotential.cpp \
 $(POTENTIALFUNCTIONS)/OldRgeImprovedOneLoopPotential.cpp \
 source/PotentialEvaluation/MassesSquaredCalculator.cpp \
 source/PotentialEvaluation/OldMassesSquaredCalculator.cpp \
 source/PotentialEvaluation/ParameterFunctionoid.cpp \
 source/PotentialEvaluation/PotentialFunction.cpp \
 source/PotentialEvaluation/ThermalFunctions.cpp \
 $(LAGRANGIANPARAMETERMANAGEMENT)/ParameterUpdatePropagator.cpp \
 $(LAGRANGIANPARAMETERMANAGEMENT)/RunningParameterManager.cpp \
 $(LAGRANGIANPARAMETERMANAGEMENT)/SlhaManager.cpp \
 $(GRADIENTBASEDMINIMIZATION)/MinuitPotentialMinimizer.cpp \
 $(GRADIENTBASEDMINIMIZATION)/PotentialForMinuit.cpp \
 $(STARTINGPOINTGENERATION)/PolynomialAtFixedScalesSolver.cpp \
 $(STARTINGPOINTGENERATION)/PolynomialSystemSolver.cpp \
 $(HOMOTOPYCONTINUATION)/BasicPolynomialHomotopyContinuation.cpp \
 $(HOMOTOPYCONTINUATION)/FieldPolynomialsWithScale.cpp \
 $(HOMOTOPYCONTINUATION)/Hom4ps2Runner.cpp \
 $(HOMOTOPYCONTINUATION)/OldHom4ps2Runner.cpp \
 $(HOMOTOPYCONTINUATION)/HomotopyContinuationSolver.cpp \
 $(HOMOTOPYCONTINUATION)/HomotopyContinuationTargetSystem.cpp \
 $(HOMOTOPYCONTINUATION)/PolynomialGradientTargetSystem.cpp \
 source/PotentialMinimization/GradientFromStartingPoints.cpp \
 source/PotentialMinimization/GradientMinimizer.cpp \
 source/PotentialMinimization/PotentialMinimizer.cpp \
 source/PotentialMinimization/PotentialMinimum.cpp \
 $(BOUNCEACTIONTUNNELING)/BounceAlongPathWithThreshold.cpp \
 $(BOUNCEACTIONTUNNELING)/CosmoTransitionsRunner.cpp \
 $(BOUNCEACTIONTUNNELING)/ThermalActionFitter.cpp \
 source/TunnelingCalculation/BounceActionTunneler.cpp \
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
