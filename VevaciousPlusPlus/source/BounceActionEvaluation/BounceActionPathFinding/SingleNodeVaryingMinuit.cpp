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
                                              size_t const movesPerImprovement,
                                             unsigned int const minuitStrategy,
                                          double const minuitToleranceFraction,
                                             double const nodeMoveThreshold ) :
    MinuitPathFinder( movesPerImprovement,
                      minuitStrategy,
                      minuitToleranceFraction ),
    potentialFunction( potentialFunction ),
    pathFactory( pathFactory ),
    pathNodes( pathFactory->GetNodesFromParameterization() ),
    currentMinuitResults(),
    currentNodeIndex( 0 ),
    nodeMoveThresholdSquared( nodeMoveThreshold * nodeMoveThreshold )
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
    pathFactory->SetVacua( falseVacuum,
                           trueVacuum );
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
    pathNodes.PathNodeSet( nodeSet,
                           *pathParameterization );
    // Now we can set the current path.
    SetCurrentPathPointer( (*pathFactory)( nodeSet,
                                 pathNodes.ParameterizationForNodes( nodeSet ),
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
    currentMinuitResults.resize( nodeSet.size() );
    for( size_t adjustmentIndex( pathNodes.AdjustmentOrderStartIndex() );
         adjustmentIndex <= pathNodes.AdjustmentOrderEndIndex();
         ++adjustmentIndex )
    {
      size_t const
      pathIndex( pathNodes.PathIndexFromAdjustmentIndex( adjustmentIndex ) );
      pathNodes.SetNodeFromNodeVector( pathIndex,
                                       nodeSet[ pathIndex ] );
      pathNodes.ExtractSingleNodeParameterization( nodeParameterization,
                                                   pathIndex,
                                                   *pathParameterization );
      currentMinuitResults[ pathIndex ] = MinuitMinimum( nodeParameterization,
                                                         initialStepSizes );
    }
  }

  // This allows Minuit2 to move each node a set number of times to try to
  // minimize the potential at that node or bounce action along the adjusted
  // path, and then sets the path from the set of nodes.
  void SingleNodeVaryingMinuit::ImprovePath()
  {
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "SingleNodeVaryingMinuit::ImprovePath() called.";
    std::cout << std::endl;*/

    bool pathConverged( true );
    std::vector< double > minuitNode( pathNodes.NumberOfFields() );

    for( size_t adjustmentIndex( pathNodes.AdjustmentOrderStartIndex() );
         adjustmentIndex <= pathNodes.AdjustmentOrderEndIndex();
         ++adjustmentIndex )
    {
      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "adjustmentIndex = " << adjustmentIndex;
      std::cout << std::endl;*/

      currentNodeIndex
      = pathNodes.PathIndexFromAdjustmentIndex( adjustmentIndex );

      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "currentNodeIndex = " << currentNodeIndex;
      std::cout << std::endl;*/

      MinuitMinimum const&
      nodeState( currentMinuitResults[ currentNodeIndex ] );
      ROOT::Minuit2::MnMigrad mnMigrad( (*this),
                                        nodeState.VariableValues(),
                                        nodeState.VariableErrors(),
                                        minuitStrategy );
      MinuitMinimum
      minuitMinimum( pathNodes.ZeroNodeParameterization().size(),
                     mnMigrad( movesPerImprovement,
                               currentMinuitTolerance ) );

      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "minuitMinimum = " << minuitMinimum.AsDebuggingString();
      std::cout << std::endl;*/

      currentMinuitResults[ currentNodeIndex ] = minuitMinimum;
      pathNodes.NodeProposal( minuitNode,
                              currentNodeIndex,
                              minuitMinimum.VariableValues() );
      // Now we check to see if the node from Minuit2 is inconsistent with the
      // hope that this call of ImprovePath():
      if( pathConverged )
      {
        if( pathNodes.PreviousNodesDependOnSubsequentNodesInAdjustmentOrder() )
        {
          if( NodeMovedMoreThanThreshold( minuitNode ) )
          {
            pathConverged = false;
          }
        }
        else
        {
          if( !(minuitMinimum.IsValidMinimum()) )
          {
            pathConverged = false;
          }
        }
      }

      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "Before setting node [ " << currentNodeIndex << " ], pathConverged = "
      << pathConverged;
      std::cout << std::endl;*/

      // After possibly comparing the node from Minuit2 to the last values it
      // had, it can be updated (as long as Minuit2 didn't go crazy and just
      // try only NAN parameterizations):
      if( !(NanParameterFromMinuit( minuitMinimum.VariableValues()) ) )
      {
        pathNodes.SetNodeFromNodeVector( currentNodeIndex,
                                         minuitNode );
      }
    }

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "After checking every node, pathConverged = " << pathConverged;
    std::cout << std::endl;*/

    if( pathConverged )
    {
      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "pathNodes.HasFurtherRefinementMode() = "
      << pathNodes.HasFurtherRefinementMode();
      std::cout << std::endl;*/

      if( pathNodes.HasFurtherRefinementMode() )
      {
        pathNodes.ConvertToRefinementOrder();
        pathCanBeImproved = true;
      }
      else
      {
        pathCanBeImproved = false;
      }
    }

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "After checking for further refinement mode, pathCanBeImproved = "
    << pathCanBeImproved;
    std::cout << std::endl;*/

    SetCurrentPathPointer( (*pathFactory)( pathNodes.PathNodes(),
                   pathNodes.ParameterizationForNodes( pathNodes.PathNodes() ),
                                           pathTemperature ) );
  }

} /* namespace VevaciousPlusPlus */
