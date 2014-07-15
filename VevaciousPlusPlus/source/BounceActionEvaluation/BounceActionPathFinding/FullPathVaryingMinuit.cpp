/*
 * FullPathVaryingMinuit.cpp
 *
 *  Created on: Jul 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/FullPathVaryingMinuit.hpp"

namespace VevaciousPlusPlus
{

  FullPathVaryingMinuit::FullPathVaryingMinuit( TunnelPathFactory* pathFactory,
                                           std::string const& xmlArguments  ) :
    MinuitPathFinder( xmlArguments ),
    pathFactory( pathFactory ),
    currentMinuitResult()
  {
    // This constructor is just an initialization list.
  }

  FullPathVaryingMinuit::~FullPathVaryingMinuit()
  {
    delete pathFactory;
  }


  // This resets the BouncePathFinder so that it sets up currentPath as its
  // initial path between the given vacua. It also resets pathCanBeImproved
  // and sets pathTemperature appropriately.
  void
  FullPathVaryingMinuit::SetInitialPath( PotentialMinimum const& falseVacuum,
                                         PotentialMinimum const& trueVacuum,
                                         TunnelPath const* startingPath,
                                         double const pathTemperature )
  {
    this->pathTemperature = pathTemperature;
    // First we need the initial parameterization of the path: either provided,
    // or the default "zero parameterization".
    std::vector< double > const*
    pathParameterization( &(pathFactory->ZeroParameterization()) );
    if( startingPath != NULL )
    {
      pathParameterization = &(startingPath->PathParameterization());
    }

    // Now we can set the current path.
    SetCurrentPathPointer( (*pathFactory)( *pathParameterization,
                                           pathTemperature ) );
    // However, we still need to set up the initial state for Minuit.
    currentMinuitTolerance = ( minuitToleranceFraction
                               * ( falseVacuum.PotentialValue()
                                   - trueVacuum.PotentialValue() ) );
    currentMinuitResult = MinuitMinimum( *pathParameterization,
                                         pathFactory->InitialStepSizes() );
  }


} /* namespace VevaciousPlusPlus */
