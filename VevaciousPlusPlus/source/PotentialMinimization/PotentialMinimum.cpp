/*
 * PotentialMinimum.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialMinimization/PotentialMinimum.hpp"

namespace VevaciousPlusPlus
{

  PotentialMinimum::PotentialMinimum(
                               std::vector< double > const& fieldConfiguration,
                                      double const potentialDepth ) :
    MinuitMinimum()
  {
    variableValues = fieldConfiguration;
    variableErrors = std::vector< double >( fieldConfiguration.size(),
                                            0.0 );
    functionValue = potentialDepth;
    functionError = 0.0;
  }

  PotentialMinimum::PotentialMinimum( MinuitMinimum const& minuitMinimum ) :
    MinuitMinimum( minuitMinimum )
  {
    // This constructor is just an initialization list.
  }

  PotentialMinimum::PotentialMinimum() :
    MinuitMinimum()
  {
    // This constructor is just an initialization list.
  }

  PotentialMinimum::PotentialMinimum( PotentialMinimum const& copySource ) :
    MinuitMinimum( copySource )
  {
    // This constructor is just an initialization list.
  }

  PotentialMinimum::~PotentialMinimum()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
