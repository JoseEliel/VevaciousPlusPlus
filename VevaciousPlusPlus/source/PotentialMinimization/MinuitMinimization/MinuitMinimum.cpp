/*
 * MinuitMinimum.cpp
 *
 *  Created on: Apr 11, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialMinimization/MinuitMinimization/MinuitMinimum.hpp"

namespace VevaciousPlusPlus
{

  MinuitMinimum::MinuitMinimum( size_t numberOfVariables,
                        ROOT::Minuit2::FunctionMinimum const& minuitMinimum ) :
    variableValues( numberOfVariables ),
    variableErrors( numberOfVariables ),
    functionValue( minuitMinimum.Fval() ),
    functionError( minuitMinimum.Edm() )
  {
    for( size_t variableIndex( 0 );
         variableIndex < numberOfVariables;
         ++variableIndex )
    {
      variableValues[ variableIndex ]
      = minuitMinimum.UserParameters().Value( variableIndex );
      variableErrors[ variableIndex ]
      = minuitMinimum.UserParameters().Error( variableIndex );
    }
  }

  MinuitMinimum::~MinuitMinimum()
  {
    // This does nothing.
  }

  MinuitMinimum::MinuitMinimum() :
    variableValues(),
    variableErrors(),
    functionValue( NAN ),
    functionError( NAN )
  {
    // This constructor is just an initialization list.
  }

} /* namespace VevaciousPlusPlus */
