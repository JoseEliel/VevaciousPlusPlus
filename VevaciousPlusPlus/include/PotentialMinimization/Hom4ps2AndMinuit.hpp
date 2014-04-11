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


    // This finds all the extrema of polynomialPotential with HOM4PS2, then
    // uses them as starting points for Minuit2 to minimize potentialForMinuit,
    // evaluated at a temperature given by temperatureInGev, and records the
    // found minima in foundMinima. It also sets dsbVacuum (using
    // polynomialPotential.DsbFieldValues() as a starting point for Minuit2),
    // and records the minima lower than dsbVacuum in panicVacua, and of those,
    // it sets panicVacuum to be the minimum in panicVacua closest to
    // dsbVacuum.
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
    MinuitManager minuitManager;
    // MINUIT thing!

    // This uses HOM4PS2 to fill purelyRealSolutionSets with all the extrema
    // of polynomialPotential.TargetPolynomialGradient().
    void FindTreeLevelExtrema();

    void WriteHom4p2Input( std::string const& hom4ps2InputFilename );

    void ParseHom4ps2Output( std::string const& hom4ps2OutputFilename );

    // This uses Minuit2 to minimize potentialForMinuit starting from the
    // values in purelyRealSolutionSets.
    void RollAndSortExtrema();
  };




  // This finds all the extrema of polynomialPotential with HOM4PS2, then uses
  // them as starting points for Minuit2 to minimize potentialForMinuit,
  // evaluated at a temperature given by temperatureInGev, and records the
  // found minima in foundMinima. It also sets dsbVacuum (using
  // polynomialPotential.DsbFieldValues() as a starting point for Minuit2), and
  // records the minima lower than dsbVacuum in panicVacua, and of those, it
  // sets panicVacuum to be the minimum in panicVacua closest to dsbVacuum.
  inline void Hom4ps2AndMinuit::FindMinima( double const temperatureInGev )
  {
    FindTreeLevelExtrema();
    potentialForMinuit.SetTemperature( temperatureInGev );
    dsbVacuum = minuitManager( polynomialPotential.DsbFieldValues() );
    RollAndSortExtrema();
  }
} /* namespace VevaciousPlusPlus */
#endif /* HOM4PS2ANDMINUIT_HPP_ */
