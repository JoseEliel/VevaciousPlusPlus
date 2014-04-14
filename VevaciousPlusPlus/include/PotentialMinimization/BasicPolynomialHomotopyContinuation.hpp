/*
 * BasicPolynomialHomotopyContinuation.hpp
 *
 *  Created on: Apr 14, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BASICPOLYNOMIALHOMOTOPYCONTINUATION_HPP_
#define BASICPOLYNOMIALHOMOTOPYCONTINUATION_HPP_

#include "../StandardIncludes.hpp"
#include "HomotopyContinuationSolver.hpp"
#include "../PotentialEvaluation/PotentialFunctions/PotentialFunctions.hpp"


namespace VevaciousPlusPlus
{

  class BasicPolynomialHomotopyContinuation : public HomotopyContinuationSolver
  {
  public:
    BasicPolynomialHomotopyContinuation(
                     HomotopyContinuationReadyPotential& polynomialPotential );
    virtual
    ~BasicPolynomialHomotopyContinuation();


    // This should find all the extrema of homotopyContinuationPotential and
    // put the purely real solutions into purelyRealSolutionSets.
    virtual void FindTreeLevelExtrema(
                std::vector< std::vector< double > >& purelyRealSolutionSets );
  };

} /* namespace VevaciousPlusPlus */
#endif /* BASICPOLYNOMIALHOMOTOPYCONTINUATION_HPP_ */
