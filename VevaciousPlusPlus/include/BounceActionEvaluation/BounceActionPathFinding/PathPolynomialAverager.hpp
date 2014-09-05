/*
 * PathPolynomialAverager.hpp
 *
 *  Created on: Sep 5, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef PATHPOLYNOMIALAVERAGER_HPP_
#define PATHPOLYNOMIALAVERAGER_HPP_

#include "CommonIncludes.hpp"
#include "MinuitBetweenPaths.hpp"
#include "../BounceActionCalculator.hpp"
#include "../PathParameterization/TunnelPath.hpp"
#include "../PathParameterization/LinearSplineThroughNodes.hpp"

namespace VevaciousPlusPlus
{

  class PathPolynomialAverager : public MinuitBetweenPaths
  {
  public:
    PathPolynomialAverager( size_t const degreeOfPolynomial,
                    BounceActionCalculator const* const bounceActionCalculator,
                           unsigned int const minuitStrategy,
                           double const minuitToleranceFraction,
                           size_t const movesPerImprovement );
    virtual ~PathPolynomialAverager();


  protected:
    size_t const degreeOfPolynomial;


    // This takes minuitParameters as a set of coefficients for weightings and
    // creates a node path of weighted averages of each node in curvedPath with
    // the corresponding node in straightPath, where the weighting for the
    // straightPath node is
    // minuitParameters[ 0 ] + ( minuitParameters[ 1 ] * the fraction along the
    // path of the node) + ( minuitParameters[ 1 ] * that fraction squared )
    // + ..., and the weighting for the curvedPath node is 1.0 - the other
    // weighting, then it returns the path of straight lines between the
    // composed nodes.
    TunnelPath const* PathForParameterization(
                         std::vector< double > const& minuitParameters ) const;

    // This prepares currentMinuitResult to be a flat 90% of curvedPath plus
    // 10% of straightPath, with the coefficients given step sizes of 0.1 each.
    void PrepareMinuitStartingPoint();
  };




  // This prepares currentMinuitResult to be a flat 90% of curvedPath plus
  // 10% of straightPath, with the coefficients given step sizes of 0.1 each.
  inline void PathPolynomialAverager::PrepareMinuitStartingPoint()
  {
    currentMinuitResult
    = MinuitMinimum( std::vector< double >( ( degreeOfPolynomial + 1 ),
                                            0.0 ),
                     std::vector< double >( ( degreeOfPolynomial + 1 ),
                                            0.1 ) );
    currentMinuitResult.VariableValues()[ 0 ] = 0.9;
  }

} /* namespace VevaciousPlusPlus */
#endif /* PATHPOLYNOMIALAVERAGER_HPP_ */
