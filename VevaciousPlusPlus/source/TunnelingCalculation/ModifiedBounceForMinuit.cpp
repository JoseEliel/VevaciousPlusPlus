/*
 * ModifiedBounceForMinuit.cpp
 *
 *  Created on: May 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{
  ModifiedBounceForMinuit::ModifiedBounceForMinuit(
                                    PotentialFunction const& potentialFunction,
                                         size_t const numberOfVaryingPathNodes,
                                       size_t const numberOfSplinesInPotential,
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                double const falseVacuumEvaporationTemperature,
                                      size_t const undershootOvershootAttempts,
                                   size_t const maximumMultipleOfLongestLength,
                                  double const initialFractionOfShortestLength,
                                              double const minimumScaleSquared,
                                   double const shootingCloseEnoughThreshold) :
    BounceActionForMinuit( potentialFunction,
                           numberOfVaryingPathNodes,
                           numberOfSplinesInPotential,
                           falseVacuum,
                           trueVacuum,
                           falseVacuumEvaporationTemperature,
                           undershootOvershootAttempts,
                           maximumMultipleOfLongestLength,
                           initialFractionOfShortestLength,
                           minimumScaleSquared,
                           shootingCloseEnoughThreshold )
  {
    // This constructor is just an initialization list.
  }

  ModifiedBounceForMinuit::~ModifiedBounceForMinuit()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
