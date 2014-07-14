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
#include "PotentialMinimization/PotentialMinimizer.hpp"
#include "TunnelingCalculation/TunnelingCalculator.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialEvaluation/PotentialFunctions/PotentialFromPolynomialAndMasses.hpp"
#include "PotentialMinimization/StartingPointFinder.hpp"
#include "PotentialMinimization/HomotopyContinuation/Hom4ps2Runner.hpp"
#include "PotentialMinimization/GradientMinimizer.hpp"
#include "PotentialMinimization/GradientBasedMinimization/MinuitPotentialMinimizer.hpp"
#include "TunnelingCalculation/BounceActionTunneling/CosmoTransitionsRunner.hpp"
#include "TunnelingCalculation/BounceActionTunneling/BounceAlongPathWithThreshold.hpp"
#include "TunnelingCalculation/BounceActionTunneling/BounceActionCalculator.hpp"
#include "TunnelingCalculation/BounceActionTunneling/BubbleShootingOnSpline.hpp"
#include "TunnelingCalculation/BounceActionTunneling/BouncePathFinder.hpp"
#include "TunnelingCalculation/BounceActionTunneling/MinuitNodePotentialMinimizer.hpp"
#include "TunnelingCalculation/BounceActionTunneling/MinuitPathBounceMinimizer.hpp"
#include "TunnelingCalculation/BounceActionTunneling/MinuitPathPotentialMinimizer.hpp"

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

    // This is the constructor which creates the components based on the data
    // in the initialization file called initializationFileName.
    VevaciousPlusPlus( std::string const& initializationFileName );

    virtual ~VevaciousPlusPlus();

    void RunPoint( std::string const& parameterFilename );

    void WriteXmlResults( std::string const& xmlFilename );

    void WriteSlhaResults( std::string const& slhaFilename,
                           bool const writeWarnings = true );


  protected:
    PotentialFromPolynomialAndMasses* potentialFunction;
    PotentialFromPolynomialAndMasses* deleterForPotentialFunction;
    RunningParameterManager runningParameterManager;
    PotentialMinimizer* potentialMinimizer;
    PotentialMinimizer* deleterForPotentialMinimizer;
    TunnelingCalculator* tunnelingCalculator;
    TunnelingCalculator* deleterForTunnelingCalculator;
    time_t currentTime;
  };

} /* namespace VevaciousPlusPlus */
#endif /* VEVACIOUSPLUSPLUS_HPP_ */
