/*
 * HomotopyContinuationAndGradient.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef HOMOTOPYCONTINUATIONANDGRADIENT_HPP_
#define HOMOTOPYCONTINUATIONANDGRADIENT_HPP_

#include "CommonIncludes.hpp"
#include "PotentialMinimizer.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"

namespace VevaciousPlusPlus
{
  class HomotopyContinuationAndGradient : public PotentialMinimizer
  {
  public:
    HomotopyContinuationAndGradient(
                                   PotentialFunction const& potentialFunction,
                      HomotopyContinuationSolver& homotopyContinuationSolver );
    virtual
    ~HomotopyContinuationAndGradient();


  protected:
    HomotopyContinuationSolver& homotopyContinuationSolver;
    std::vector< std::vector< double > > purelyRealSolutionSets;
  };

} /* namespace VevaciousPlusPlus */
#endif /* HOMOTOPYCONTINUATIONANDGRADIENT_HPP_ */
