/*
 * MinuitPathFinder.cpp
 *
 *  Created on: Jul 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinuitPathFinder.hpp"

namespace VevaciousPlusPlus
{

  double const MinuitPathFinder::functionValueForNanInput(
                                        std::numeric_limits< double >::max() );

  MinuitPathFinder::MinuitPathFinder( size_t const movesPerImprovement,
                                      unsigned int const minuitStrategy,
                                      double const minuitToleranceFraction ) :
    BouncePathFinder(),
    ROOT::Minuit2::FCNBase(),
    movesPerImprovement( movesPerImprovement ),
    minuitStrategy( minuitStrategy ),
    minuitToleranceFraction( minuitToleranceFraction ),
    currentMinuitTolerance( NAN )
  {
    // This constructor is just an initialization list.
  }

  MinuitPathFinder::~MinuitPathFinder()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
