/*
 * PotentialForMinuit.hpp
 *
 *  Created on: Apr 9, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POTENTIALFORMINUIT_HPP_
#define POTENTIALFORMINUIT_HPP_

#include "../StandardIncludes.hpp"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "../PotentialEvaluation/PotentialFunction.hpp"

namespace VevaciousPlusPlus
{

  class PotentialForMinuit : public ROOT::Minuit2::FCNBase
  {
  public:
    PotentialForMinuit( PotentialFunction& minimizationFunction );
    virtual
    ~PotentialForMinuit();


    // This implements operator() for FCNBase, the function that MINUIT will
    // minimize.
    virtual double
    operator()( std::vector< double > const& fieldConfiguration ) const
    { return minimizationFunction( fieldConfiguration,
                                   currentTemperature ); }

    // This implements Up() for FCNBase just to stick to a basic value.
    virtual double Up() const { return 1.0; }

    // The potential is minimized at a fixed temperature, so this sets the
    // temperature.
    void SetTemperature( double const currentTemperature )
    { this->currentTemperature = currentTemperature; }


  protected:
    PotentialFunction& minimizationFunction;
    double currentTemperature;
  };

} /* namespace VevaciousPlusPlus */
#endif /* POTENTIALFORMINUIT_HPP_ */
