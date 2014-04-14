/*
 * HomotopyContinuationSolver.hpp
 *
 *  Created on: Apr 14, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef HOMOTOPYCONTINUATIONSOLVER_HPP_
#define HOMOTOPYCONTINUATIONSOLVER_HPP_

#include "../StandardIncludes.hpp"
#include "../PotentialEvaluation/PotentialFunctions/PotentialFunctions.hpp"


namespace VevaciousPlusPlus
{

  class HomotopyContinuationSolver
  {
  public:
    HomotopyContinuationSolver(
     HomotopyContinuationReadyPotential const& homotopyContinuationPotential );
    virtual
    ~HomotopyContinuationSolver();


    // This should find all the extrema of homotopyContinuationPotential and
    // put the purely real solutions into purelyRealSolutionSets.
    virtual void FindTreeLevelExtrema(
            std::vector< std::vector< double > >& purelyRealSolutionSets ) = 0;


  protected:
    HomotopyContinuationReadyPotential const& homotopyContinuationPotential;
  };

} /* namespace VevaciousPlusPlus */
#endif /* HOMOTOPYCONTINUATIONSOLVER_HPP_ */
