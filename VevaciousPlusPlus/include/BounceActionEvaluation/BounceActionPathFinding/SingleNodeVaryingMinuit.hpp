/*
 * SingleNodeVaryingMinuit.hpp
 *
 *  Created on: Jul 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SINGLENODEVARYINGMINUIT_HPP_
#define SINGLENODEVARYINGMINUIT_HPP_

#include "CommonIncludes.hpp"
#include "MinuitPathFinder.hpp"
#include "Minuit2/MnMigrad.h"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/GradientBasedMinimization/MinuitMinimum.hpp"
#include "../PathParameterization/PathFromNodesFactory.hpp"
#include "../PathParameterization/NodesFromParameterization.hpp"

namespace VevaciousPlusPlus
{

  class SingleNodeVaryingMinuit : public MinuitPathFinder
  {
  public:
    SingleNodeVaryingMinuit( PotentialFunction const& potentialFunction,
                             PathFromNodesFactory* const pathFactory,
                             size_t const movesPerImprovement = 100,
                             unsigned int const minuitStrategy = 1,
                             double const minuitToleranceFraction = 0.5 );
    virtual ~SingleNodeVaryingMinuit();


    // This resets the BouncePathFinder and return a pointer to an initial path
    // between the given vacua, and the management of this memory is entirely
    // up to the calling code (if we allowed ourselves to require a
    // C++11-compliant compiler, this function would return a
    // unique_ptr< TunnelPath >). It should also reset
    // pathCanBeImproved and set pathTemperature appropriately.
    virtual TunnelPath const*
    SetInitialPath( PotentialMinimum const& falseVacuum,
                    PotentialMinimum const& trueVacuum,
                    TunnelPath const* startingPath = NULL,
                    double const pathTemperature = 0.0 );

    // This allows Minuit2 to move each node a set number of times to try to
    // minimize the potential at that node, and then sets the path from the set
    // of nodes.
    virtual TunnelPath const* ImprovePath();


  protected:
    PotentialFunction const& potentialFunction;
    PathFromNodesFactory* const pathFactory;
    NodesFromParameterization& pathNodes;
    std::vector< MinuitMinimum > currentMinuitResults;
    size_t currentNodeIndex;
    bool currentlyTuning;
  };

} /* namespace VevaciousPlusPlus */
#endif /* SINGLENODEVARYINGMINUIT_HPP_ */
