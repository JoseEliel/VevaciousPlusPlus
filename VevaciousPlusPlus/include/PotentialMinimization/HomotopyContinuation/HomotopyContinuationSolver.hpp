/*
 * HomotopyContinuationSolver.hpp
 *
 *  Created on: Apr 14, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef HOMOTOPYCONTINUATIONSOLVER_HPP_
#define HOMOTOPYCONTINUATIONSOLVER_HPP_

#include "CommonIncludes.hpp"
#include "../StartingPointFinder.hpp"
#include "HomotopyContinuationTargetSystem.hpp"


namespace VevaciousPlusPlus
{

  class HomotopyContinuationSolver : public StartingPointFinder
  {
  public:
    HomotopyContinuationSolver(
       HomotopyContinuationTargetSystem const& homotopyContinuationPotential );
    virtual
    ~HomotopyContinuationSolver();


    // This should find all the extrema of homotopyContinuationPotential and
    // put the purely real solutions into startingPoints (inherited from
    // StartingPointFinder, not actually over-written here).
    // virtual void
    // operator()( std::vector< std::vector< double > >& startingPoints ) = 0;


  protected:
    HomotopyContinuationTargetSystem const& homotopyContinuationPotential;
  };

} /* namespace VevaciousPlusPlus */
#endif /* HOMOTOPYCONTINUATIONSOLVER_HPP_ */
