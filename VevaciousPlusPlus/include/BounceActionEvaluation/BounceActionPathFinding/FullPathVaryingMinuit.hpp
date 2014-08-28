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
#include "../PathParameterization/TunnelPathFactory.hpp"

namespace VevaciousPlusPlus
{

  class FullPathVaryingMinuit : public MinuitPathFinder
  {
  public:
    FullPathVaryingMinuit( TunnelPathFactory* const pathFactory,
                           size_t const movesPerImprovement = 100,
                           unsigned int const minuitStrategy = 1,
                           double const minuitToleranceFraction = 0.5 );
    virtual ~FullPathVaryingMinuit();


    // This resets the BouncePathFinder so that it sets up currentPath as its
    // initial path between the given vacua. It also resets pathCanBeImproved
    // and sets pathTemperature appropriately.
    virtual TunnelPath const*
    SetInitialPath( PotentialMinimum const& falseVacuum,
                    PotentialMinimum const& trueVacuum,
                    TunnelPath const* startingPath = NULL,
                    double const pathTemperature = 0.0 );

    // This allows Minuit2 to adjust the full path a set number of times to try
    // to minimize the sum of potentials at a set of nodes or bounce action
    // along the adjusted path, and then sets the path.
    virtual TunnelPath const* ImprovePath();


  protected:
    TunnelPathFactory* const pathFactory;
    MinuitMinimum currentMinuitResult;
  };




  // This resets the BouncePathFinder so that it sets up currentPath as its
  // initial path between the given vacua. It also resets pathCanBeImproved
  // and sets pathTemperature appropriately.
  inline TunnelPath const*
  FullPathVaryingMinuit::SetInitialPath( PotentialMinimum const& falseVacuum,
                                         PotentialMinimum const& trueVacuum,
                                         TunnelPath const* startingPath,
                                         double const pathTemperature )
  {
    pathFactory->SetVacua( falseVacuum,
                           trueVacuum );
    this->pathTemperature = pathTemperature;
    // First we need the initial parameterization of the path: either provided,
    // or the default zero parameterization; then we can set the current path.
    std::vector< double > const*
    pathParameterization( &(pathFactory->ZeroParameterization()) );
    if( startingPath != NULL )
    {
      pathParameterization = &(startingPath->PathParameterization());
    }
    // However, we still need to set up the initial state for Minuit.
    currentMinuitTolerance = ( minuitToleranceFraction
                               * ( falseVacuum.PotentialValue()
                                   - trueVacuum.PotentialValue() ) );
    currentMinuitResult = MinuitMinimum( *pathParameterization,
                                         pathFactory->InitialStepSizes() );
    return (*pathFactory)( *pathParameterization,
                           pathTemperature );
  }

  // This allows Minuit2 to adjust the full path a set number of times to try
  // to minimize the sum of potentials at a set of nodes or bounce action
  // along the adjusted path, and then sets the path.
  inline TunnelPath const* FullPathVaryingMinuit::ImprovePath()
  {
    ROOT::Minuit2::MnMigrad mnMigrad( (*this),
                                      currentMinuitResult.VariableValues(),
                                      currentMinuitResult.VariableErrors(),
                                      minuitStrategy );
    MinuitMinimum minuitMinimum( currentMinuitResult.VariableValues().size(),
                                 mnMigrad( movesPerImprovement,
                                           currentMinuitTolerance ) );
    currentMinuitResult = minuitMinimum;
    return (*pathFactory)( minuitMinimum.VariableValues(),
                           pathTemperature );
  }

} /* namespace VevaciousPlusPlus */
#endif /* FULLPATHVARYINGMINUIT_HPP_ */
