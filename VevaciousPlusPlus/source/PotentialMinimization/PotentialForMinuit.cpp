/*
 * PotentialForMinuit.cpp
 *
 *  Created on: Apr 9, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"
#include <omp.h>

namespace VevaciousPlusPlus
{

  PotentialForMinuit::PotentialForMinuit(
                                    PotentialFunction& minimizationFunction ) :
    ROOT::Minuit2::FCNBase(),
    minimizationFunction( minimizationFunction ),
    currentTemperature( 0.0 )
  {
    // This constructor is just an initialization list.
  }

  PotentialForMinuit::~PotentialForMinuit()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
