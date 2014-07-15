/*
 * SingleNodeVaryingMinuit.cpp
 *
 *  Created on: Jul 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/SingleNodeVaryingMinuit.hpp"

namespace VevaciousPlusPlus
{

  SingleNodeVaryingMinuit::SingleNodeVaryingMinuit(
                                    PotentialFunction const& potentialFunction,
                                             PathFromNodesFactory* pathFactory,
                                           std::string const& xmlArguments  ) :
    MinuitPathFinder( xmlArguments ),
    potentialFunction( potentialFunction ),
    pathFactory( pathFactory ),
    pathNodes( pathFactory->NodesFromParameterization() ),
    currentMinuitResults(),
    currentNodeIndex( 0 )
  {
    // This constructor is just an initialization list.
  }

  SingleNodeVaryingMinuit::~SingleNodeVaryingMinuit()
  {
    delete pathFactory;
  }


  // This resets the BouncePathFinder so that it sets up currentPath as its
  // initial path between the given vacua. It also resets pathCanBeImproved
  // and sets pathTemperature appropriately.
  void
  SingleNodeVaryingMinuit::SetInitialPath( PotentialMinimum const& falseVacuum,
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

    // Now we need the initial path node set.
    std::vector< std::vector< double > >
    nodeSet( pathNodes.NumberOfPathNodes(),
             std::vector< double >( pathNodes.NumberOfFields() ) );
    pathFactory->SetVacua( falseVacuum,
                           trueVacuum );
    pathNodes.PathNodeSet( nodeSet,
                           *pathParameterization );
    // Now we can set the current path.
    SetCurrentPathPointer( (*pathFactory)( nodeSet,
                                           pathTemperature ) );
    // However, we still need to set up the initial states for Minuit as well
    // as ensure that pathNodes records the nodes too, so that we can vary the
    // nodes individually.
    currentMinuitTolerance = ( minuitToleranceFraction
                               * ( falseVacuum.PotentialValue()
                                   - trueVacuum.PotentialValue() ) );
    std::vector< double > nodeParameterization;
    std::vector< double > initialStepSizes;
    pathNodes.SetSingleNodeStepSizes( initialStepSizes );

    // Now we note all the new positions of the nodes, in field space and also
    // by parameterization.
    currentMinuitResults.clear();
    for( size_t nodeIndex( 1 );
         nodeIndex <= pathNodes.NumberOfVaryingNodes();
         ++nodeIndex )
    {
      pathNodes.SetNodeInAdjustmentOrder( nodeIndex,
                                          nodeSet[ nodeIndex ] );
      pathNodes.ExtractSingleNodeParameterization( nodeParameterization,
                                                   nodeIndex,
                                                   *pathParameterization );
      currentMinuitResults.push_back( MinuitMinimum( nodeParameterization,
                                                     initialStepSizes ) );
    }
  }

} /* namespace VevaciousPlusPlus */
