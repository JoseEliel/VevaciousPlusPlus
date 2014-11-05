/*
 * BubbleRadialValueDescription.cpp
 *
 *  Created on: May 22, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BubbleRadialValueDescription.hpp"

namespace VevaciousPlusPlus
{

  BubbleRadialValueDescription::BubbleRadialValueDescription() :
    radialValue( -1.0 ),
    auxiliaryValue( -2.0 ),
    auxiliarySlope( 0.0 )
  {
    // This constructor is just an initialization list.
  }

  BubbleRadialValueDescription::BubbleRadialValueDescription(
                                                      double const radialValue,
                                                   double const auxiliaryValue,
                                                double const auxiliarySlope ) :
    radialValue( radialValue ),
    auxiliaryValue( auxiliaryValue ),
    auxiliarySlope( auxiliarySlope )
  {
    // This constructor is just an initialization list.
  }

  BubbleRadialValueDescription::BubbleRadialValueDescription(
                             BubbleRadialValueDescription const& copySource ) :
    radialValue( copySource.radialValue ),
    auxiliaryValue( copySource.auxiliaryValue ),
    auxiliarySlope( copySource.auxiliarySlope )
  {
    // This constructor is just an initialization list.
  }

  BubbleRadialValueDescription::~BubbleRadialValueDescription()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
