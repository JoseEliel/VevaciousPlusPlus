/*
 * MinuitNodePotentialMinimizer.cpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/BounceActionTunneling/MinuitNodePotentialMinimizer.hpp"

namespace VevaciousPlusPlus
{

  MinuitNodePotentialMinimizer::MinuitNodePotentialMinimizer(
                                    PotentialFunction const& potentialFunction,
                                             PathFromNodesFactory* pathFactory,
                                          NodesFromParameterization* pathNodes,
                                           std::string const& xmlArguments  ) :
    BouncePathFinder(),
    ROOT::Minuit2::FCNBase(),
    potentialFunction( potentialFunction ),
    pathFactory( pathFactory ),
    pathNodes( pathNodes ),
    currentMinuitResults(),
    movesPerNodePerImprovement( 100 ),
    minuitStrategy( 1 ),
    minuitToleranceFraction( 1.0 ),
    currentMinuitTolerance( NAN ),
    currentNodeIndex( 0 )
  {
    BOL::AsciiXmlParser argumentParser;
    while( argumentParser.readNextElement() )
    {
      if( argumentParser.currentElementNameMatches( "MinuitStrategy" ) )
      {
        minuitStrategy = BOL::StringParser::stringToInt(
                            argumentParser.getTrimmedCurrentElementContent() );
      }
      else if( argumentParser.currentElementNameMatches( "MinuitTolerance" ) )
      {
        minuitToleranceFraction = BOL::StringParser::stringToDouble(
                            argumentParser.getTrimmedCurrentElementContent() );
      }
    }
  }

  MinuitNodePotentialMinimizer::~MinuitNodePotentialMinimizer()
  {
    delete pathFactory;
    delete pathNodes;
  }


  // This should reset the BouncePathFinder so that it sets up currentPath as
  // its initial path between the given vacua. It should also reset
  // pathCanBeImproved and set pathTemperature appropriately.
  void MinuitNodePotentialMinimizer::SetInitialPath(
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                                TunnelPath const* startingPath,
                                                 double const pathTemperature )
  {
    this->pathTemperature = pathTemperature;
    // First we need the initial parameterization of the path: either provided,
    // or the default "zero parameterization".
    std::vector< double > const*
    pathParameterization( &(pathNodes->ZeroParameterization()) );
    if( startingPath != NULL )
    {
      pathParameterization = &(startingPath->PathParameterization());
    }

    // Now we need the initial path node set.
    std::vector< std::vector< double > >
    nodeSet( pathNodes->NumberOfPathNodes(),
             std::vector< double >( pathNodes->NumberOfFields() ) );
    pathFactory->SetVacua( falseVacuum,
                           trueVacuum );
    pathNodes->PathNodeSet( nodeSet,
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
    pathNodes->SetSingleNodeStepSizes( initialStepSizes );

    // Now we note all the new positions of the nodes, in field space and also
    // by parameterization.
    currentMinuitResults.clear();
    for( size_t nodeIndex( 1 );
         nodeIndex <= pathNodes->NumberOfVaryingNodes();
         ++nodeIndex )
    {
      pathNodes->SetNodeInAdjustmentOrder( nodeIndex,
                                           nodeSet[ nodeIndex ] );
      pathNodes->ExtractSingleNodeParameterization( nodeParameterization,
                                                    nodeIndex,
                                                    *pathParameterization );
      currentMinuitResults.push_back( MinuitMinimum( nodeParameterization,
                                                     initialStepSizes ) );
    }
  }

} /* namespace VevaciousPlusPlus */
