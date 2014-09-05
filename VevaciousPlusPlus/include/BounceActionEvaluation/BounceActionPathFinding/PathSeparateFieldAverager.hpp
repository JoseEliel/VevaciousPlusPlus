/*
 * PathSeparateFieldAverager.hpp
 *
 *  Created on: Sep 5, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef PATHSEPARATEFIELDAVERAGER_HPP_
#define PATHSEPARATEFIELDAVERAGER_HPP_

#include "CommonIncludes.hpp"
#include "MinuitBetweenPaths.hpp"
#include "../BounceActionCalculator.hpp"
#include "../PathParameterization/TunnelPath.hpp"
#include "../PathParameterization/LinearSplineThroughNodes.hpp"

namespace VevaciousPlusPlus
{

  class PathSeparateFieldAverager : public MinuitBetweenPaths
  {
  public:
    PathSeparateFieldAverager(
                          BounceActionCalculator* const bounceActionCalculator,
                               unsigned int const minuitStrategy,
                               double const minuitToleranceFraction,
                               size_t const movesPerImprovement );
    virtual ~PathSeparateFieldAverager();


  protected:
    // This takes minuitParameters as the set of fractions for each field which
    // should be taken as weightings for the straightPath node's value for that
    // field while the weighting for the field from the node from curvedPath is
    // 1.0 - that weighting, then it returns the path of straight lines between
    // the composed nodes.
    TunnelPath const* PathForParameterization(
                         std::vector< double > const& minuitParameters ) const;

    // This prepares currentMinuitResult to be 90% of curvedPath in each field,
    // with the coefficients given step sizes of 0.1 each.
    void PrepareMinuitStartingPoint();
  };




  // This prepares currentMinuitResult to be a flat 90% of curvedPath plus
  // 10% of straightPath, with the coefficients given step sizes of 0.1 each.
  inline void PathSeparateFieldAverager::PrepareMinuitStartingPoint()
  {
    currentMinuitResult
    = MinuitMinimum( std::vector< double >( numberOfFields,
                                            0.1 ),
                     std::vector< double >( numberOfFields,
                                            0.1 ) );
  }
} /* namespace VevaciousPlusPlus */
#endif /* PATHSEPARATEFIELDAVERAGER_HPP_ */
