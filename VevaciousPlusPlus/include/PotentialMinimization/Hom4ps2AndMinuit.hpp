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
#include "BOLlib/include/BOLlib.hpp"
#include "PotentialForMinuit.hpp"
#include "MinuitManager.hpp"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnPrint.h"

namespace VevaciousPlusPlus
{

  class Hom4ps2AndMinuit : public HomotopyContinuationAndGradient
  {
  public:
    Hom4ps2AndMinuit( HomotopyContinuationReadyPolynomial& polynomialPotential,
                      std::string const& pathToHom4ps2,
                      std::string const homotopyType = "1" );
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
    PotentialForMinuit potentialForMinuit;
    std::string const pathToHom4ps2;
    std::string const homotopyType;
    BOL::StringParser variableNamer;
    std::vector< std::complex< long double > > complexSolutions;
    std::vector< std::string > variableNames;
    std::map< std::string, unsigned int > nameToIndexMap;
    std::vector< std::vector< double > > purelyRealSolutionSets;
    // MINUIT thing!

    void WriteHom4p2Input( std::string const& hom4ps2InputFilename );

    void ParseHom4ps2Output( std::string const& hom4ps2OutputFilename );
  };

} /* namespace VevaciousPlusPlus */
#endif /* HOM4PS2ANDMINUIT_HPP_ */
