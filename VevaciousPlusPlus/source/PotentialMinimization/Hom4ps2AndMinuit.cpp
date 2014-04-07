/*
 * Hom4ps2AndMinuit.cpp
 *
 *  Created on: Apr 7, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  Hom4ps2AndMinuit::Hom4ps2AndMinuit(
                      HomotopyContinuationReadyPolynomial& polynomialPotential,
                                      std::string const pathToHom4ps2 ) :
    HomotopyContinuationAndGradient( polynomialPotential ),
    polynomialPotential( polynomialPotential ),
    pathToHom4ps2( pathToHom4ps2 )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "Hom4ps2AndMinuit::Hom4ps2AndMinuit( ..., \"" << pathToHom4ps2
    << "\" )";
    std::cout << std::endl;/**/
  }

  Hom4ps2AndMinuit::~Hom4ps2AndMinuit()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "Hom4ps2AndMinuit::~Hom4ps2AndMinuit()";
    std::cout << std::endl;/**/
  }



  // This should find all the minima of the potential evaluated at a
  // temperature given by temperatureInGev, and record them in foundMinima.
  // It should also set dsbVacuum, and record the minima lower than dsbVacuum
  // in panicVacua, and of those, it should set panicVacuum to be the prime
  // candidate for tunneling out of dsbVacuum (by default, taken to be the
  // minimum in panicVacua closest to dsbVacuum).
  void Hom4ps2AndMinuit::FindMinima( double const temperatureInGev )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "Hom4ps2AndMinuit::FindMinima( " << temperatureInGev << " )";
    std::cout << std::endl;/**/
  }

  // This should find the minimum at temperature temperatureInGev nearest to
  // minimumToAdjust (which is assumed to be a minimum of the potential at a
  // different temperature).
  PotentialMinimum Hom4ps2AndMinuit::AdjustMinimumForTemperature(
                                       PotentialMinimum const& minimumToAdjust,
                                                double const temperatureInGev )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "Hom4ps2AndMinuit::AdjustMinimumForTemperature( ..., "
    << temperatureInGev << " )";
    std::cout << std::endl;
    return PotentialMinimum( minimumToAdjust );/**/
  }

} /* namespace VevaciousPlusPlus */
