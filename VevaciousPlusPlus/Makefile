#
# Makefile for VevaciousPlusPlus.
#     Add your name if you modify
#     Author: Ben O'Leary (benjamin.oleary@gmail.com)
#     Date: 2014
#
#     Author: José Eliel Camargo-Molina (Elielcamargomolina@gmail.com)
#     Date: Sept. 2018
#
#      This file is part of VevaciousPlusPlus, released under the GNU General
#      Public License. Please see the accompanying README.VevaciousPlusPlus.txt
#      file for a full list of files, brief documentation on how to use these
#      classes, and further details on the license.
#

SHELL=/bin/sh

CC=g++-7
AR=ar
ARFLAGS=rcs
#BOOSTDIR=./Boost/
BOOSTDIR=/usr/local/include/boost/
#EIGENDIR=./Eigen/
EIGENDIR=/usr/local/include/eigen3/
#MINUITLIBNAME=lcg_Minuit
MINUITLIBNAME=Minuit2
#MINUITINCLUDEDIR=./MINUIT/include/
MINUITINCLUDEDIR=/usr/local/include/
#MINUITLIBDIR=./MINUIT/include/
MINUITLIBDIR=/usr/local/lib
DEBUGFLAGS=
INCLUDES=-I./include/ -I./include/LHPC/ \
 -I$(EIGENDIR) -I$(BOOSTDIR) -I$(MINUITINCLUDEDIR)
CCFLAGS=-Wall -Wno-unused-local-typedefs -O3 -fopenmp -fPIC -std=c++11 $(DEBUGFLAGS) $(INCLUDES)
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
 $(LAGRANGIANPARAMETERMANAGEMENT)/SARAHManager.cpp \
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
 source/BOSS_factory_VevaciousPlusPlus.cpp \
 source/BOSS_VevaciousPlusPlus.cpp \
 source/BOSS_wrapperutils.cpp \
 source/VevaciousPlusPlus.cpp
OBJECTS=$(SOURCES:.cpp=.o)
BINDIR=./bin/
LIBNAME=VevaciousPlusPlus
LIBDIR=./lib/
LIBRARYFILE=$(LIBDIR)/lib$(LIBNAME).a
SHAREDLIBRARYFILE=$(LIBDIR)/lib$(LIBNAME).so
LDFLAGS= -L$(LIBDIR) -l$(LIBNAME) $(MINUITLIBDIR)/lib$(MINUITLIBNAME).a
#-L$(MINUITLIBDIR) -l:lib$(MINUITLIBNAME).a

# -l$(LHPCLIBNAME) -L$(MINUITLIBDIR) -l:lib$(MINUITLIBNAME).a
# There is some silly screwing around to force linking to the static Minuit2
# library.
MAINCPP=source/VevaciousPlusPlusMain.cpp
EXEFILE=$(BINDIR)/VevaciousPlusPlus.exe

all: $(OBJECTS) $(LIBRARYFILE) $(SHAREDLIBRARYFILE) $(EXEFILE)

shared: $(SHAREDLIBRARYFILE)

.PHONY: all clean

$(OBJECTS): %.o: %.cpp
	$(CC) $(CCFLAGS) $(CCOBJECTFLAGS) $< -o $@

$(LIBDIR):
	mkdir $@

$(LIBRARYFILE): $(LIBDIR) $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)

$(SHAREDLIBRARYFILE): $(OBJECTS)
	$(CC) $(CCFLAGS) -shared -o $@ $^ $(MINUITLIBDIR)/lib$(MINUITLIBNAME).a

$(BINDIR):
	mkdir $@

$(EXEFILE): $(BINDIR) $(LIBRARYFILE) $(SHAREDLIBRARYFILE)
	$(CC) $(CCFLAGS) $(CCEXEFLAGS) -o $@ $(MAINCPP) $(LDFLAGS)

clean:
	\rm -f $(OBJECTS) $(LIBRARYFILE) $(SHAREDLIBRARYFILE) $(EXEFILE)

deleteobjects:
	\rm -f $(OBJECTS)