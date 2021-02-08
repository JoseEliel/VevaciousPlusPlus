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


cmake_minimum_required(VERSION 2.8.12)
project(VevaciousPlusPlus)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

cmake_policy(SET CMP0042 NEW)

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

# Set the general compiler flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-local-typedefs -O3 -fPIC -fopenmp")


#############################################################################
# EXTERNAL PROJECTS
#############################################################################
include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

message("${CMAKE_C_COMPILER}")

if(WITHIN_GAMBIT)

# Nothing

else()

# MINUIT2
set(Minuit_name "minuit2")
set(Minuit_lib_name "Minuit2")
set(Minuit_ver "5.34.14")
set(Minuit_include "${PROJECT_SOURCE_DIR}/${Minuit_name}/${Minuit_ver}/include/")
#set(Minuit_dl "https://www.cern.ch/mathlibs/sw/5_34_14/Minuit2/Minuit2-5.34.14.tar.gz")
# Official Minuit2 website seems to be down. Temporary fix. 
set(Minuid_dl "https://web.archive.org/web/20160322211910/http://project-mathlibs.web.cern.ch/project-mathlibs/sw/Minuit2/Minuit2.tar.gz")
set(Minuit_dir "${PROJECT_SOURCE_DIR}/${Minuit_name}/${Minuit_ver}")
ExternalProject_Add(${Minuit_name}_${Minuit_ver}
            #DOWNLOAD_COMMAND wget -qO- ${Minuit_dl} | tar xvz -C ${Minuit_dir} --strip-components 1
            URL ${Minuit_dl} 
            SOURCE_DIR ${Minuit_dir}
            BUILD_IN_SOURCE 1
            CONFIGURE_COMMAND ./configure --prefix=${Minuit_dir} CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} CFLAGS=${CMAKE_C_FLAGS} CXXFLAGS=${CMAKE_CXX_FLAGS}
            BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
            INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
        )
endif()

# HOM4PS
set(name "hom4ps")
set(hom4ps_ver "2.0")
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(dl "http://www.math.nsysu.edu.tw/~leetsung/works/HOM4PS_soft_files/HOM4PS2_MacOSX.tar.gz")
else()
    set(dl "http://www.math.nsysu.edu.tw/~leetsung/works/HOM4PS_soft_files/HOM4PS2_64-bit.tar.gz")
endif()

set(dir "${PROJECT_SOURCE_DIR}/${name}/${hom4ps_ver}")

ExternalProject_Add(${name}_${hom4ps_ver}
        #DOWNLOAD_COMMAND wget -qO- ${dl} | tar xvz -C ${dir} --strip-components 1
        URL ${dl} 
        SOURCE_DIR ${dir}
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        )

# PHC
set(name "phc")
set(phc_ver "2.4.66")

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(dl "http://www.math.uic.edu/~jan/mactel64y_phcv24p.tar.gz")
else()
    set(dl "http://www.math.uic.edu/~jan/x86_64phcv24p.tar.gz")
endif()

set(dir "${PROJECT_SOURCE_DIR}/${name}/${phc_ver}")


ExternalProject_Add(${name}_${phc_ver}
        #DOWNLOAD_COMMAND wget -qO- ${dl} | tar xvz -C ${dir}
        URL ${dl} 
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
        source/PotentialEvaluation/PotentialFunctions/TreeLevelPotential.cpp
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
# GENERAL INCLUDE DIRS
#############################################################################

include_directories("${Minuit_include}")
include_directories(include)
include_directories(include/LHPC)


#############################################################################
# OBJECT LIBRARY (For making the objects used in executable and dynamic library)
#############################################################################

add_library(objlib OBJECT ${sources})

if(WITHIN_GAMBIT)
  
  #Nothing    

else()

  add_dependencies(objlib ${Minuit_name}_${Minuit_ver})

endif()


#############################################################################
# VEVACIOUS EXECUTABLE
#############################################################################

add_executable(VevaciousPlusPlus
        $<TARGET_OBJECTS:objlib>)

# Adding dependencies, including Homotopy continuation codes
# This is skipped if within the GAMBIT ecosystem


if(WITHIN_GAMBIT)
  
  #Nothing    

else()

  add_dependencies(VevaciousPlusPlus hom4ps_${hom4ps_ver})

  #add_dependencies(VevaciousPlusPlus phc_${phc_ver})

  add_dependencies(VevaciousPlusPlus objlib)

  add_dependencies(VevaciousPlusPlus VevaciousPlusPlus-lib)

endif()


# Linking to minuit

# Linking to minuit
if(WITHIN_GAMBIT)
    #Nothing 
else()
    set(Minuit_lib "${PROJECT_SOURCE_DIR}/${Minuit_name}/${Minuit_ver}/lib/")
endif()

target_link_libraries(VevaciousPlusPlus ${Minuit_lib}/libMinuit2.a)



#############################################################################
# VEVACIOUS LIBRARY
#############################################################################

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/)

add_library(VevaciousPlusPlus-lib SHARED $<TARGET_OBJECTS:objlib>)

# Linking to minuit

target_link_libraries(VevaciousPlusPlus-lib ${Minuit_lib}/libMinuit2.a)

add_dependencies(VevaciousPlusPlus-lib objlib)

# Setting correct name for library

set_target_properties(VevaciousPlusPlus-lib
        PROPERTIES OUTPUT_NAME VevaciousPlusPlus)




#############################################################################
# Writing Paths to Initialization Files
#############################################################################

if(EXISTS "${PROJECT_SOURCE_DIR}/InitializationFiles/MSSMInitialization.xml")

message("Initialization files have been modified or already copied from the templates, skipping writing new files")

else()

message("Writing Input and Initialization files from templates")

set(vevacious_path ${PROJECT_SOURCE_DIR})

# Here I copy all Initialization File templates into the InitializationFiles folder
file(COPY ${vevacious_path}/Template_InitializationFiles/ DESTINATION ${vevacious_path}/InitializationFiles/)

# Here I make a list of all xml recursively within InitializationFiles
file(GLOB_RECURSE allinitxmls "${vevacious_path}/InitializationFiles/*.xml" )

# This replaces the ${vevacious_path} variable in all xmls in InitializationFiles.
foreach(xml ${allinitxmls})
    message("${xml}")
    configure_file( ${xml}  ${xml} )
endforeach()

# Here I copy all the input files template into bin folder
file(COPY ${vevacious_path}/Template_inputFiles/ DESTINATION ${vevacious_path}/bin/)

# Here I make a list of all xml recursively within the bin folder
file(GLOB_RECURSE allinputxmls "${vevacious_path}/bin/*.xml" )

# This replaces the ${vevacious_path} variable in all xmls in InitializationFiles.
foreach(xml ${allinputxmls})
    message("${xml}")
    configure_file( ${xml}  ${xml} )
endforeach()

endif()
