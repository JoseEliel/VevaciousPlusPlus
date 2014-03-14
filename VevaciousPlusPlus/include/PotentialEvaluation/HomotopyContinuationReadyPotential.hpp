/*
 * HomotopyContinuationReadyPotential.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef HOMOTOPYCONTINUATIONREADYPOTENTIAL_HPP_
#define HOMOTOPYCONTINUATIONREADYPOTENTIAL_HPP_

#include "../StandardIncludes.hpp"
#include "PotentialFunction.hpp"

namespace VevaciousPlusPlus
{
  class HomotopyContinuationReadyPotential : public PotentialFunction
  {
  public:
    HomotopyContinuationReadyPotential();
    virtual
    ~HomotopyContinuationReadyPotential();


    // This should evaluate the target system and place the values in
    // destinationVector.
    virtual void
    HomotopyContinuationSystemValues(
                                   std::vector< double > solutionConfiguration,
                                std::vector< double >& destinationVector ) = 0;

    // This should evaluate the derivatives of the target system and place the
    // values in destinationMatrix.
    virtual void
    HomotopyContinuationSystemGradients(
                                   std::vector< double > solutionConfiguration,
                 std::vector< std::vector< double > >& destinationMatrix ) = 0;
  };

} /* namespace VevaciousPlusPlus */
#endif /* HOMOTOPYCONTINUATIONREADYPOTENTIAL_HPP_ */
