/*
 * HomotopyContinuationSolver.hpp
 *
 *  Created on: Apr 14, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef HOMOTOPYCONTINUATIONSOLVER_HPP_
#define HOMOTOPYCONTINUATIONSOLVER_HPP_

#include "../../CommonIncludes.hpp"
#include "HomotopyContinuationTargetSystem.hpp"


namespace VevaciousPlusPlus
{

  class HomotopyContinuationSolver
  {
  public:
    HomotopyContinuationSolver(
     HomotopyContinuationTargetSystem const& homotopyContinuationPotential );
    virtual
    ~HomotopyContinuationSolver();


    // This should find all the extrema of homotopyContinuationPotential and
    // put the purely real solutions into purelyRealSolutionSets.
    virtual void FindTreeLevelExtrema(
            std::vector< std::vector< double > >& purelyRealSolutionSets ) = 0;


  protected:
    HomotopyContinuationTargetSystem const& homotopyContinuationPotential;
  };

} /* namespace VevaciousPlusPlus */
#endif /* HOMOTOPYCONTINUATIONSOLVER_HPP_ */
