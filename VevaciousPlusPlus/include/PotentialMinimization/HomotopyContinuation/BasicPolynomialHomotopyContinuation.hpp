/*
 * BasicPolynomialHomotopyContinuation.hpp
 *
 *  Created on: Apr 14, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BASICPOLYNOMIALHOMOTOPYCONTINUATION_HPP_
#define BASICPOLYNOMIALHOMOTOPYCONTINUATION_HPP_

#include "CommonIncludes.hpp"
#include "HomotopyContinuationSolver.hpp"
#include "HomotopyContinuationTargetSystem.hpp"


namespace VevaciousPlusPlus
{

  class BasicPolynomialHomotopyContinuation : public HomotopyContinuationSolver
  {
  public:
    BasicPolynomialHomotopyContinuation(
                     HomotopyContinuationTargetSystem& polynomialPotential );
    virtual ~BasicPolynomialHomotopyContinuation();


    // This should find all the extrema of homotopyContinuationPotential and
    // put the purely real solutions into purelyRealSolutionSets.
    virtual void FindTreeLevelExtrema(
                std::vector< std::vector< double > >& purelyRealSolutionSets );
  };

} /* namespace VevaciousPlusPlus */
#endif /* BASICPOLYNOMIALHOMOTOPYCONTINUATION_HPP_ */
