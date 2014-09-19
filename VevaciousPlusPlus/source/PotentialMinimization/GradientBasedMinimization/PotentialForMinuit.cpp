/*
 * PotentialForMinuit.cpp
 *
 *  Created on: Apr 9, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialMinimization/GradientBasedMinimization/PotentialForMinuit.hpp"

namespace VevaciousPlusPlus
{

  PotentialForMinuit::PotentialForMinuit(
                              PotentialFunction const& minimizationFunction ) :
    ROOT::Minuit2::FCNBase(),
    minimizationFunction( minimizationFunction ),
    fieldOrigin( minimizationFunction.NumberOfFieldVariables(),
                 0.0 ),
    functionAtOrigin( minimizationFunction( fieldOrigin ) ),
    currentTemperature( 0.0 )
  {
    // This constructor is just an initialization list.
  }

  PotentialForMinuit::~PotentialForMinuit()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
