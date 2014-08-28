/*
 * MinuitNodePotentialMinimizer.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITNODEPOTENTIALMINIMIZER_HPP_
#define MINUITNODEPOTENTIALMINIMIZER_HPP_

#include "CommonIncludes.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "../PathParameterization/PathFromNodesFactory.hpp"
#include "SingleNodeVaryingMinuit.hpp"

namespace VevaciousPlusPlus
{
  class MinuitNodePotentialMinimizer : public SingleNodeVaryingMinuit
  {
  public:
    MinuitNodePotentialMinimizer( PathFromNodesFactory* const pathFactory,
                                  PotentialFunction const& potentialFunction,
                                  size_t const movesPerImprovement = 100,
                                  unsigned int const minuitStrategy = 1,
                                  double const minuitToleranceFraction = 0.5,
                                  double const nodeMoveThreshold = 0.01 );
    virtual ~MinuitNodePotentialMinimizer();


    // It may seem unwise to have this object call Minuit on itself, but really
    // it's just a handy way of keeping the minimization function within the
    // class that ends up finding its minimum. In this case,
    // nodeParameterization is just the parameterization of the node.
    virtual double
    operator()( std::vector< double > const& nodeParameterization ) const;
  };




  // It may seem unwise to have this object call Minuit on itself, but really
  // it's just a handy way of keeping the minimization function within the
  // class that ends up finding its minimum. The node index is set by
  // ImprovePath before the Minuit minimization and the temperature is set by
  // SetInitialPath.
  inline double MinuitNodePotentialMinimizer::operator()(
                      std::vector< double > const& nodeParameterization ) const
  {
    if( NanParameterFromMinuit( nodeParameterization ) )
    {
      return functionValueForNanInput;
    }
    std::vector< double > nodeVector( pathNodes.NumberOfFields() );
    pathNodes.NodeProposal( nodeVector,
                            currentNodeIndex,
                            nodeParameterization );
    return potentialFunction( nodeVector,
                              pathTemperature );
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINUITNODEPOTENTIALMINIMIZER_HPP_ */
