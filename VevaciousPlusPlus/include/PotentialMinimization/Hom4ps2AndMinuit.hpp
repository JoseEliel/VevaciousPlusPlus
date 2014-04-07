/*
 * Hom4ps2AndMinuit.hpp
 *
 *  Created on: Apr 7, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef HOM4PS2ANDMINUIT_HPP_
#define HOM4PS2ANDMINUIT_HPP_

#include "../StandardIncludes.hpp"
#include "HomotopyContinuationAndGradient.hpp"
#include "../PotentialEvaluation/HomotopyContinuationReadyPolynomial.hpp"

namespace VevaciousPlusPlus
{

  class Hom4ps2AndMinuit : public HomotopyContinuationAndGradient
  {
  public:
    Hom4ps2AndMinuit( HomotopyContinuationReadyPolynomial& polynomialPotential,
                      std::string const pathToHom4ps2 );
    virtual
    ~Hom4ps2AndMinuit();


    // This should find all the minima of the potential evaluated at a
    // temperature given by temperatureInGev, and record them in foundMinima.
    // It should also set dsbVacuum, and record the minima lower than dsbVacuum
    // in panicVacua, and of those, it should set panicVacuum to be the prime
    // candidate for tunneling out of dsbVacuum (by default, taken to be the
    // minimum in panicVacua closest to dsbVacuum).
    virtual void FindMinima( double const temperatureInGev = 0.0 );

    // This should find the minimum at temperature temperatureInGev nearest to
    // minimumToAdjust (which is assumed to be a minimum of the potential at a
    // different temperature).
    virtual PotentialMinimum
    AdjustMinimumForTemperature( PotentialMinimum const& minimumToAdjust,
                                 double const temperatureInGev );


  protected:
    HomotopyContinuationReadyPolynomial& polynomialPotential;
    std::string const pathToHom4ps2;
    // MINUIT thing!
  };

} /* namespace VevaciousPlusPlus */
#endif /* HOM4PS2ANDMINUIT_HPP_ */
