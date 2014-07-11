/*
 * MinuitNodePotentialMinimizer.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITNODEPOTENTIALMINIMIZER_HPP_
#define MINUITNODEPOTENTIALMINIMIZER_HPP_

#include "CommonIncludes.hpp"
#include "BouncePathFinder.hpp"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/GradientBasedMinimization/MinuitMinimum.hpp"
#include "PathFromNodesFactory.hpp"
#include "NodesFromParameterization.hpp"
#include "LinearSplineThroughNodesFactory.hpp"
#include "QuadraticSplineThroughNodesFactory.hpp"
#include "PolynomialThroughNodesFactory.hpp"

namespace VevaciousPlusPlus
{

  class MinuitNodePotentialMinimizer : public BouncePathFinder,
                                       public ROOT::Minuit2::FCNBase
  {
  public:
    MinuitNodePotentialMinimizer( PotentialFunction const& potentialFunction,
                                  std::string const& xmlArguments );
    virtual ~MinuitNodePotentialMinimizer();


    // This class implements operator() and Up() inherited from
    // ROOT::Minuit2::FCNBase. There is no deadly diamond of doom because
    // BouncePathFinder does not have an operator() or Up() at all.

    // It may seem unwise to have this object call Minuit on itself, but really
    // it's just a handy way of keeping the minimization function within the
    // class that ends up finding its minimum. In this case,
    // nodeParameterization is just the parameterization of the node.
    virtual double
    operator()( std::vector< double > const& nodeParameterization ) const;

    // This implements Up() for ROOT::Minuit2::FCNBase just to stick to a basic
    // value.
    virtual double Up() const { return 1.0; }

    // This should reset the BouncePathFinder so that it sets up currentPath as
    // its initial path between the given vacua. It should also reset
    // pathCanBeImproved and set pathTemperature appropriately.
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
    PathFromNodesFactory* pathFactory;
    NodesFromParameterization* pathNodes;
    std::vector< MinuitMinimum > currentMinuitResults;
    size_t movesPerNodePerImprovement;
    unsigned int minuitStrategy;
    double minuitToleranceFraction;
    double currentMinuitTolerance;
    size_t currentNodeIndex;
  };




  // It may seem unwise to have this object call Minuit on itself, but really
  // it's just a handy way of keeping the minimization function within the
  // class that ends up finding its minimum. The node index is set by
  // ImprovePath before the Minuit minimization and the temperature is set by
  // SetInitialPath.
  inline double MinuitNodePotentialMinimizer::operator()(
                      std::vector< double > const& nodeParameterization ) const
  {
    std::vector< double > nodeVector( pathNodes->NumberOfFields() );
    pathNodes->NodeProposal( nodeVector,
                             currentNodeIndex,
                             nodeParameterization );
    return potentialFunction( nodeVector,
                              pathTemperature );
  }

  // This allows Minuit2 to move each node a set number of times to try to
  // minimize the potential at that node, and then sets the path from the set
  // of nodes.
  inline void MinuitNodePotentialMinimizer::ImprovePath()
  {
    for( size_t nodeIndex( 1 );
         nodeIndex <= pathNodes->NumberOfVaryingNodes();
         ++nodeIndex )
    {
      // We have to set which node Minuit will vary in each iteration.
      currentNodeIndex = nodeIndex;
      MinuitMinimum const& nodeState( currentMinuitResults[ nodeIndex - 1 ] );
      ROOT::Minuit2::MnMigrad mnMigrad( (*this),
                                        nodeState.VariableValues(),
                                        nodeState.VariableErrors(),
                                        minuitStrategy );
      MinuitMinimum minuitMinimum( mnMigrad( movesPerNodePerImprovement,
                                             currentMinuitTolerance ) );
      currentMinuitResults[ nodeIndex - 1 ] = minuitMinimum;
      pathNodes->SetNodeInAdjustmentOrder( nodeIndex,
                                           minuitMinimum.VariableValues() );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINUITNODEPOTENTIALMINIMIZER_HPP_ */
