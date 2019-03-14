# VevaciousPlusPlus
#************************************************
#
#  Master CMake configuration script for VevaciousPlusPlus.
#
#
#  Authors (add name and date if you modify):
#
#  \author José Eliel Camargo-Molina
#          (elielcamargomolina@gmail.com)
#
#
#************************************************


cmake_minimum_required(VERSION 3.12)
project(VevaciousPlusPlus)

#############################################################################
# CHECKS
#############################################################################

# Check for C++11 support
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")
else()
    message("${CMAKE_CXX_COMPILER} does not support the C++11 standard. Please use a different C++ compiler.")
endif()


# Check for Eigen
include(cmake/FindEigen3.cmake)
find_package(Eigen3 3.1.0)
if(EIGEN3_FOUND)
    include_directories("${EIGEN3_INCLUDE_DIR}")
    message("-- Eigen version: ${EIGEN3_VERSION}")
else()
    message("${BoldRed}   Eigen v3.1.0 or greater not found. This is required for Vevacious to work. ${ColourReset}")
endif()

# Check for Boost
find_package(Boost 1.41)
if(Boost_FOUND)
    include_directories("${Boost_INCLUDE_DIR}")
else()
    message("${BoldRed} Boost libraries 1.41 or greater not found. This is required for Vevacious to work. ${ColourReset}")
endif()

#############################################################################
# EXTERNAL PROJECTS
#############################################################################
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

# MINUIT2
set(name "minuit2")
set(dl "https://www.cern.ch/mathlibs/sw/5_34_14/Minuit2/Minuit2-5.34.14.tar.gz")
set(ver "5.34.14")
set(dir "${PROJECT_SOURCE_DIR}/../${name}/${ver}")
ExternalProject_Add(${name}_${ver}
            DOWNLOAD_COMMAND wget -qO- ${dl} | tar xvz -C ${dir} --strip-components 1
            SOURCE_DIR ${dir}
            BUILD_IN_SOURCE 1
            CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} ./configure --prefix=${dir}
            BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
            INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
        )


# HOM4PS
set(name "hom4ps")
set(ver "2.0")
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(dl "http://www.math.nsysu.edu.tw/~leetsung/works/HOM4PS_soft_files/HOM4PS2_MacOSX.tar.gz")
else()
    set(dl "http://www.math.nsysu.edu.tw/~leetsung/works/HOM4PS_soft_files/HOM4PS2_64-bit.tar.gz")
endif()

set(dir "${PROJECT_SOURCE_DIR}/../${name}/${ver}")

ExternalProject_Add(${name}_${ver}
        DOWNLOAD_COMMAND wget -qO- ${dl} | tar xvz -C ${dir} --strip-components 1
        SOURCE_DIR ${dir}
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        )

# PHC
set(name "phc")
set(ver "2.4.58")

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(dl "http://www.math.uic.edu/~jan/mactel64y_phcv24p.tar.gz")
else()
    set(dl "http://www.math.uic.edu/~jan/x86_64phcv24p.tar.gz")
endif()

set(dir "${PROJECT_SOURCE_DIR}/../${name}/${ver}")


ExternalProject_Add(${name}_${ver}
        DOWNLOAD_COMMAND wget -qO- ${dl} | tar xvz -C ${dir}
        SOURCE_DIR ${dir}
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        )

set(sources source/BounceActionEvaluation/BounceActionPathFinding/MinuitOnHypersurfaces.cpp
        source/BounceActionEvaluation/BounceActionPathFinding/MinuitOnPotentialOnParallelPlanes.cpp
        source/BounceActionEvaluation/BounceActionPathFinding/MinuitOnPotentialPerpendicularToPath.cpp
        source/BounceActionEvaluation/PathParameterization/LinearSplineThroughNodes.cpp
        source/BounceActionEvaluation/BubbleShootingOnPathInFieldSpace.cpp
        source/BounceActionEvaluation/SplinePotential.cpp
        source/BounceActionEvaluation/UndershootOvershootBubble.cpp
        source/LagrangianParameterManagement/LesHouchesAccordBlockEntryManager.cpp
        source/LagrangianParameterManagement/LhaLinearlyInterpolatedBlockEntry.cpp
        source/LagrangianParameterManagement/LhaPolynomialFitBlockEntry.cpp
        source/LagrangianParameterManagement/SARAHManager.cpp
        source/LagrangianParameterManagement/SlhaBlocksWithSpecialCasesManager.cpp
        source/LagrangianParameterManagement/SlhaCompatibleWithSarahManager.cpp
        source/PotentialEvaluation/BuildingBlocks/ParametersAndFieldsProductTerm.cpp
        source/PotentialEvaluation/MassesSquaredCalculators/ComplexMassSquaredMatrix.cpp
        source/PotentialEvaluation/MassesSquaredCalculators/RealMassesSquaredMatrix.cpp
        source/PotentialEvaluation/MassesSquaredCalculators/SymmetricComplexMassMatrix.cpp
        source/PotentialEvaluation/PotentialFunctions/FixedScaleOneLoopPotential.cpp
        source/PotentialEvaluation/PotentialFunctions/PotentialFromPolynomialWithMasses.cpp
        source/PotentialEvaluation/PotentialFunctions/RgeImprovedOneLoopPotential.cpp
        source/PotentialEvaluation/ThermalFunctions.cpp
        source/PotentialMinimization/HomotopyContinuation/Hom4ps2Runner.cpp
        source/PotentialMinimization/HomotopyContinuation/PHCRunner.cpp
        source/PotentialMinimization/StartingPointGeneration/PolynomialAtFixedScalesSolver.cpp
        source/PotentialMinimization/StartingPointGeneration/PolynomialSystemSolver.cpp
        source/PotentialMinimization/GradientFromStartingPoints.cpp
        source/TunnelingCalculation/BounceActionTunneling/BounceAlongPathWithThreshold.cpp
        source/TunnelingCalculation/BounceActionTunneling/CosmoTransitionsRunner.cpp
        source/TunnelingCalculation/BounceActionTunneling/ThermalActionFitter.cpp
        source/TunnelingCalculation/BounceActionTunneler.cpp
        source/Utilities/WarningLogger.cpp
        source/VevaciousPlusPlus.cpp
        source/VevaciousPlusPlusMain.cpp)


#############################################################################
# OBJECT LIBRARY (For making the objects used in executable and dynamic library)
#############################################################################

add_library( objlib OBJECT ${sources})

#############################################################################
# VEVACIOUS EXECUTABLE
#############################################################################

# Setting compiler flags

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-local-typedefs -O3 -fPIC ")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-local-typedefs -O3 -fPIC -fopenmp")
endif()



SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

set(Minuit_name "minuit2")
set(Minuit_lib_name "Minuit2")
set(Minuit_ver "5.34.14")
set(Minuit_include "${PROJECT_SOURCE_DIR}/../${Minuit_name}/${Minuit_ver}/include/")

include_directories("${Minuit_include}")
include_directories(include)
include_directories(include/LHPC)


add_executable(VevaciousPlusPlus
        $<TARGET_OBJECTS:objlib>)

add_dependencies(VevaciousPlusPlus ${Minuit_name}_${Minuit_ver})

add_dependencies(VevaciousPlusPlus objlib)

add_dependencies(VevaciousPlusPlus VevaciousPlusPlus-lib)

# Linking to minuit 

set(Minuit_lib "${PROJECT_SOURCE_DIR}/../${Minuit_name}/${Minuit_ver}/lib/")

target_link_libraries(VevaciousPlusPlus ${Minuit_lib}/libminuit2.a)



#############################################################################
# VEVACIOUS LIBRARY
#############################################################################

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/)

add_library(VevaciousPlusPlus-lib SHARED $<TARGET_OBJECTS:objlib>)

# Linking to minuit 

target_link_libraries(VevaciousPlusPlus-lib ${Minuit_lib}/libminuit2.a)

add_dependencies(VevaciousPlusPlus objlib)

# Setting correct name for library 

set_target_properties(VevaciousPlusPlus-lib
        PROPERTIES OUTPUT_NAME VevaciousPlusPlus)
