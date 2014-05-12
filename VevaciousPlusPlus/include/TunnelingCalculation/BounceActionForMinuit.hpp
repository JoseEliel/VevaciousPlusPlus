/*
 * BounceActionForMinuit.hpp
 *
 *  Created on: May 12, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BOUNCEACTIONFORMINUIT_HPP_
#define BOUNCEACTIONFORMINUIT_HPP_

#include "../../CommonIncludes.hpp"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "../../PotentialEvaluation.hpp"

namespace VevaciousPlusPlus
{

  class BounceActionForMinuit : public ROOT::Minuit2::FCNBase
  {
  public:
    BounceActionForMinuit();
    virtual
    ~BounceActionForMinuit();


    // This implements operator() for FCNBase, the function that MINUIT will
    // minimize. The values of splineCoefficients should be sets of n
    // coefficients for polynomials for each of the n fields, plus the final
    // element of splineCoefficients should be the temperature.
    virtual double
    operator()( std::vector< double > const& splineCoefficients ) const;

    // This implements Up() for FCNBase just to stick to a basic value.
    virtual double Up() const { return 1.0; }


  protected:
    BubbleProfiler const& bubbleProfiler;
  };

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEACTIONFORMINUIT_HPP_ */
