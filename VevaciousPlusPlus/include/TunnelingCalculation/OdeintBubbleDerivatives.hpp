/*
 * OdeintBubbleDerivatives.hpp
 *
 *  Created on: Jun 5, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef ODEINTBUBBLEDERIVATIVES_HPP_
#define ODEINTBUBBLEDERIVATIVES_HPP_

#include "../CommonIncludes.hpp"
#include "../PotentialEvaluation/SimplePolynomial.hpp"
#include "BubbleRadialValueDescription.hpp"

namespace VevaciousPlusPlus
{

  class OdeintBubbleDerivatives
  {
  public:
    OdeintBubbleDerivatives(
                        PathFieldsAndPotential const& pathFieldsAndPotential );
    virtual
    ~OdeintBubbleDerivatives();

    // This puts the first and second derivatives based on
    // auxiliaryAndFirstDerivative into firstAndSecondDerivatives, in the form
    // required for the Boost odeint package.
    void operator()( std::vector< double > const& auxiliaryAndFirstDerivative,
                     std::vector< double >& firstAndSecondDerivatives,
                     double const radialValue );


  protected:
    SplinePotential const& potentialSpline;
    std::vector< SimplePolynomial > const& firstDerivatives;
    size_t const numberOfFields;
    std::vector< SimplePolynomial > secondDerivatives;
    double const dampingFactor;
  };

} /* namespace VevaciousPlusPlus */
#endif /* ODEINTBUBBLEDERIVATIVES_HPP_ */
