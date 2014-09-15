/*
 * MinuitBetweenPaths.cpp
 *
 *  Created on: Sep 2, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinuitBetweenPaths.hpp"

namespace VevaciousPlusPlus
{

  MinuitBetweenPaths::MinuitBetweenPaths(
                          BounceActionCalculator* const bounceActionCalculator,
                                          unsigned int const minuitStrategy,
                                          double const minuitToleranceFraction,
                                          size_t const movesPerImprovement ) :
    ROOT::Minuit2::FCNBase(),
    curvedPath( NULL ),
    pathTemperature( NAN ),
    bounceActionCalculator( bounceActionCalculator ),
    numberOfFields( 0 ),
    minuitStrategy( minuitStrategy ),
    movesPerImprovement( movesPerImprovement ),
    minuitToleranceFraction( minuitToleranceFraction ),
    currentMinuitTolerance( NAN ),
    currentMinuitResult()
  {
    // This constructor is just an initialization list.
  }

  MinuitBetweenPaths::~MinuitBetweenPaths()
  {
    delete bounceActionCalculator;
  }


  // This sets up straightPath to go from curvedPath.front() to
  // curvedPath.back() in a straight line with as many nodes as curvedPath has,
  // as well as updating pathTemperature and currentMinuitTolerance.
  void MinuitBetweenPaths::UpdateNodes(
                        std::vector< std::vector< double > > const& curvedPath,
                                        double const pathTemperature )
  {
    bounceActionCalculator->ResetVacua( curvedPath.front(),
                                        curvedPath.back(),
                                        pathTemperature );
    this->curvedPath = &curvedPath;
    numberOfFields = curvedPath.front().size();
    this->pathTemperature = pathTemperature;
    LinearSplineThroughNodes comparisonPath( curvedPath,
                                             std::vector< double >( 0 ),
                                             pathTemperature );
    currentMinuitTolerance = ( minuitToleranceFraction
                               * (*bounceActionCalculator)( comparisonPath ) );
    PrepareMinuitStartingPoint();
  }

} /* namespace VevaciousPlusPlus */
