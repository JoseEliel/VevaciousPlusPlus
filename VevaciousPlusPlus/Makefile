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
#EIGENDIR=./Eigen/
EIGENDIR=/home/bol/BOL/ProjectDependencies/Eigen/
#MINUITLIBNAME=lcg_Minuit
MINUITLIBNAME=Minuit2
#MINUITINCLUDEDIR=./MINUIT/include/
MINUITINCLUDEDIR=/home/bol/BOL/ProjectDependencies/Minuit2/include/
#MINUITLIBDIR=./MINUIT/include/
MINUITLIBDIR=/home/bol/BOL/ProjectDependencies/Minuit2/lib/
DEBUGFLAGS=
INCLUDES=-I./include/ -I./include/LHPC/ \
 -I$(EIGENDIR) -I$(BOOSTDIR) -I$(MINUITINCLUDEDIR)
CCFLAGS=-Wall -Wno-unused-local-typedefs -O3 -fopenmp -std=c++11 $(DEBUGFLAGS) $(INCLUDES)
CCOBJECTFLAGS=-c
CCEXEFLAGS=

LAGRANGIANPARAMETERMANAGEMENT= source/LagrangianParameterManagement/
BUILDINGBLOCKS= source/PotentialEvaluation/BuildingBlocks/
BOUNCEACTIONPATHFINDING= source/BounceActionEvaluation/BounceActionPathFinding/
PATHPARAMETERIZATION= source/BounceActionEvaluation/PathParameterization/
MASSESSQUAREDCALCULATORS= source/PotentialEvaluation/MassesSquaredCalculators/
POTENTIALFUNCTIONS= source/PotentialEvaluation/PotentialFunctions/
POTENTIALMINIMIZATION= source/PotentialMinimization/
HOMOTOPYCONTINUATION= $(POTENTIALMINIMIZATION)/HomotopyContinuation/
STARTINGPOINTGENERATION= $(POTENTIALMINIMIZATION)/StartingPointGeneration/
GRADIENTBASEDMINIMIZATION= $(POTENTIALMINIMIZATION)/GradientBasedMinimization/
BOUNCEACTIONTUNNELING= source/TunnelingCalculation/BounceActionTunneling/

SOURCES= source/Utilities/WarningLogger.cpp \
 $(LAGRANGIANPARAMETERMANAGEMENT)/LhaLinearlyInterpolatedBlockEntry.cpp \
 $(LAGRANGIANPARAMETERMANAGEMENT)/LhaPolynomialFitBlockEntry.cpp \
 $(LAGRANGIANPARAMETERMANAGEMENT)/LesHouchesAccordBlockEntryManager.cpp \
 $(LAGRANGIANPARAMETERMANAGEMENT)/SlhaBlocksWithSpecialCasesManager.cpp \
 $(LAGRANGIANPARAMETERMANAGEMENT)/SlhaCompatibleWithSarahManager.cpp \
 $(BUILDINGBLOCKS)/ParametersAndFieldsProductTerm.cpp \
 source/PotentialEvaluation/ThermalFunctions.cpp \
 $(MASSESSQUAREDCALCULATORS)/RealMassesSquaredMatrix.cpp \
 $(MASSESSQUAREDCALCULATORS)/ComplexMassSquaredMatrix.cpp \
 $(MASSESSQUAREDCALCULATORS)/SymmetricComplexMassMatrix.cpp \
 $(POTENTIALFUNCTIONS)/PotentialFromPolynomialWithMasses.cpp \
 $(POTENTIALFUNCTIONS)/FixedScaleOneLoopPotential.cpp \
 $(POTENTIALFUNCTIONS)/RgeImprovedOneLoopPotential.cpp \
 $(STARTINGPOINTGENERATION)/PolynomialSystemSolver.cpp \
 $(STARTINGPOINTGENERATION)/PolynomialAtFixedScalesSolver.cpp \
 $(HOMOTOPYCONTINUATION)/Hom4ps2Runner.cpp \
 $(HOMOTOPYCONTINUATION)/PHCRunner.cpp \
 source/PotentialMinimization/GradientFromStartingPoints.cpp \
 $(PATHPARAMETERIZATION)/LinearSplineThroughNodes.cpp \
 source/BounceActionEvaluation/SplinePotential.cpp \
 source/BounceActionEvaluation/UndershootOvershootBubble.cpp \
 source/BounceActionEvaluation/BubbleShootingOnPathInFieldSpace.cpp \
 $(BOUNCEACTIONPATHFINDING)/MinuitOnHypersurfaces.cpp \
 $(BOUNCEACTIONPATHFINDING)/MinuitOnPotentialOnParallelPlanes.cpp \
 $(BOUNCEACTIONPATHFINDING)/MinuitOnPotentialPerpendicularToPath.cpp \
 source/TunnelingCalculation/BounceActionTunneler.cpp \
 $(BOUNCEACTIONTUNNELING)/ThermalActionFitter.cpp \
 $(BOUNCEACTIONTUNNELING)/CosmoTransitionsRunner.cpp \
 $(BOUNCEACTIONTUNNELING)/BounceAlongPathWithThreshold.cpp \
 source/VevaciousPlusPlus.cpp
OBJECTS=$(SOURCES:.cpp=.o)
BINDIR=./bin/
LIBNAME=VevaciousPlusPlus
LIBDIR=./lib/
LIBRARYFILE=$(LIBDIR)/lib$(LIBNAME).a
LDFLAGS= -L$(LIBDIR) -l$(LIBNAME) \
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
