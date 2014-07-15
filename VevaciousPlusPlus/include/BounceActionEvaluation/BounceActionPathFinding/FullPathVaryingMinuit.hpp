/*
 * FullPathVaryingMinuit.hpp
 *
 *  Created on: Jul 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef FULLPATHVARYINGMINUIT_HPP_
#define FULLPATHVARYINGMINUIT_HPP_

#include "CommonIncludes.hpp"
#include "MinuitPathFinder.hpp"
#include "Minuit2/MnMigrad.h"
#include "PotentialMinimization/GradientBasedMinimization/MinuitMinimum.hpp"
#include "TunnelPathFactory.hpp"

namespace VevaciousPlusPlus
{

  class FullPathVaryingMinuit : public MinuitPathFinder
  {
  public:
    FullPathVaryingMinuit( TunnelPathFactory* pathFactory,
                           std::string const& xmlArguments );
    virtual ~FullPathVaryingMinuit();


    // This resets the BouncePathFinder so that it sets up currentPath as its
    // initial path between the given vacua. It also resets pathCanBeImproved
    // and sets pathTemperature appropriately.
    virtual void SetInitialPath( PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 TunnelPath const* startingPath = NULL,
                                 double const pathTemperature = 0.0 );

    // This allows Minuit2 to adjust the full path a set number of times to try
    // to minimize the sum of potentials at a set of nodes or bounce action
    // along the adjusted path, and then sets the path.
    virtual void ImprovePath();


  protected:
    TunnelPathFactory* pathFactory;
    MinuitMinimum currentMinuitResult;
  };




  // This allows Minuit2 to adjust the full path a set number of times to try
  // to minimize the sum of potentials at a set of nodes or bounce action
  // along the adjusted path, and then sets the path.
  inline void FullPathVaryingMinuit::ImprovePath()
  {
    ROOT::Minuit2::MnMigrad mnMigrad( (*this),
                                      currentMinuitResult.VariableValues(),
                                      currentMinuitResult.VariableErrors(),
                                      minuitStrategy );
    MinuitMinimum minuitMinimum( mnMigrad( movesPerImprovement,
                                           currentMinuitTolerance ) );
    currentMinuitResult = minuitMinimum;
    SetCurrentPathPointer( (*pathFactory)( minuitMinimum.VariableValues(),
                                           pathTemperature ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* FULLPATHVARYINGMINUIT_HPP_ */
