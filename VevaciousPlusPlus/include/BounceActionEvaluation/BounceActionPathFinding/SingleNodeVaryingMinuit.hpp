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
#include "PathFromNodesFactory.hpp"
#include "NodesFromParameterization.hpp"

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


    // This resets the BouncePathFinder so that it sets up currentPath as its
    // initial path between the given vacua. It also resets pathCanBeImproved
    // and sets pathTemperature appropriately.
    virtual void SetInitialPath( PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 TunnelPath const* startingPath = NULL,
                                 double const pathTemperature = 0.0 );

    // This allows Minuit2 to move each node a set number of times to try to
    // minimize the potential at that node, and then sets the path from the set
    // of nodes.
    virtual void ImprovePath();


  protected:
    PotentialFunction const& potentialFunction;
    PathFromNodesFactory* const pathFactory;
    NodesFromParameterization& pathNodes;
    std::vector< MinuitMinimum > currentMinuitResults;
    size_t currentNodeIndex;
  };




  // This allows Minuit2 to move each node a set number of times to try to
  // minimize the potential at that node or bounce action along the adjusted
  // path, and then sets the path from the set of nodes.
  inline void SingleNodeVaryingMinuit::ImprovePath()
  {
    for( size_t nodeIndex( 1 );
         nodeIndex <= pathNodes.NumberOfVaryingNodes();
         ++nodeIndex )
    {
      // We have to set which node Minuit will vary in each iteration.
      currentNodeIndex = nodeIndex;
      MinuitMinimum const& nodeState( currentMinuitResults[ nodeIndex - 1 ] );
      ROOT::Minuit2::MnMigrad mnMigrad( (*this),
                                        nodeState.VariableValues(),
                                        nodeState.VariableErrors(),
                                        minuitStrategy );
      MinuitMinimum minuitMinimum( mnMigrad( movesPerImprovement,
                                             currentMinuitTolerance ) );
      currentMinuitResults[ nodeIndex - 1 ] = minuitMinimum;
      pathNodes.SetNodeInAdjustmentOrder( nodeIndex,
                                          minuitMinimum.VariableValues() );
    }
    SetCurrentPathPointer( (*pathFactory)( pathNodes.PathNodes(),
                                           pathTemperature ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* SINGLENODEVARYINGMINUIT_HPP_ */
