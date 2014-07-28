/*
 * VevaciousPlusPlus.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef VEVACIOUSPLUSPLUS_HPP_
#define VEVACIOUSPLUSPLUS_HPP_

#include "CommonIncludes.hpp"
#include "VersionInformation.hpp"
#include "SlhaManagement/RunningParameterManager.hpp"
#include "PotentialMinimization/PotentialMinimizer.hpp"
#include "TunnelingCalculation/TunnelingCalculator.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialEvaluation/PotentialFunctions/PotentialFromPolynomialAndMasses.hpp"
#include "PotentialEvaluation/PotentialFunctions/FixedScaleOneLoopPotential.hpp"
#include "PotentialEvaluation/PotentialFunctions/RgeImprovedOneLoopPotential.hpp"
#include "PotentialMinimization/StartingPointFinder.hpp"
#include "PotentialMinimization/HomotopyContinuation/Hom4ps2Runner.hpp"
#include "PotentialMinimization/GradientMinimizer.hpp"
#include "PotentialMinimization/GradientBasedMinimization/MinuitPotentialMinimizer.hpp"
#include "TunnelingCalculation/BounceActionTunneling/CosmoTransitionsRunner.hpp"
#include "TunnelingCalculation/BounceActionTunneling/BounceAlongPathWithThreshold.hpp"
#include "BounceActionEvaluation/BounceActionCalculator.hpp"
#include "BounceActionEvaluation/BubbleShootingOnSpline.hpp"
#include "BounceActionEvaluation/BouncePathFinder.hpp"
#include "BounceActionEvaluation/BounceActionPathFinding/MinuitNodePotentialMinimizer.hpp"
#include "BounceActionEvaluation/BounceActionPathFinding/MinuitPathBounceMinimizer.hpp"
#include "BounceActionEvaluation/BounceActionPathFinding/MinuitPathPotentialMinimizer.hpp"
#include "BounceActionEvaluation/PathParameterization/LinearSplineThroughNodesFactory.hpp"
#include "BounceActionEvaluation/PathParameterization/PolynomialThroughNodesFactory.hpp"
#include "BounceActionEvaluation/PathParameterization/QuadraticSplineThroughNodesFactory.hpp"
#include "BounceActionEvaluation/PathParameterization/PolynomialsFromCoefficientsFactory.hpp"

namespace VevaciousPlusPlus
{
  class VevaciousPlusPlus
  {
  public:
    // This is the constructor for those who know what they are doing, to allow
    // the main function of the program to decide the components and pass them
    // in to the constructor, allowing for custom components without having to
    // edit the VevaciousPlusPlus files. Those wishing to just use the default
    // possibilities can use the other constructor, which takes the name of an
    // initialization file, and then creates the components based on the data
    // in that file.
    VevaciousPlusPlus( SlhaManager& slhaManager,
                       PotentialMinimizer& potentialMinimizer,
                       TunnelingCalculator& tunnelingCalculator );

    // This is the constructor that we expect to be used in normal use: it
    // reads in an initialization file in XML with name given by
    // initializationFileName and then assembles the appropriate components for
    // the objects. The body of the constructor is just a lot of statements
    // reading in XML elements and creating new instances of components. It'd
    // be great to be able to use std::unique_ptrs but we're sticking to
    // allowing non-C++11-compliant compilers.
    VevaciousPlusPlus( std::string const& initializationFileName );

    virtual ~VevaciousPlusPlus();

    void RunPoint( std::string const& parameterFilename );

    void WriteXmlResults( std::string const& xmlFilename );

    void WriteSlhaResults( std::string const& slhaFilename,
                           bool const writeWarnings = true );


  protected:
    SlhaManager* slhaManager;
    SlhaManager* deleterForSlhaManager;
    PotentialFromPolynomialAndMasses* potentialFunction;
    PotentialFromPolynomialAndMasses* deleterForPotentialFunction;
    PotentialMinimizer* potentialMinimizer;
    PotentialMinimizer* deleterForPotentialMinimizer;
    TunnelingCalculator* tunnelingCalculator;
    TunnelingCalculator* deleterForTunnelingCalculator;
    time_t currentTime;
  };

} /* namespace VevaciousPlusPlus */
#endif /* VEVACIOUSPLUSPLUS_HPP_ */
