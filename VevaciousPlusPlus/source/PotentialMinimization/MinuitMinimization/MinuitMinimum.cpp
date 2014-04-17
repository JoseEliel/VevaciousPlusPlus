/*
 * MinuitMinimum.cpp
 *
 *  Created on: Apr 11, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  MinuitMinimum::MinuitMinimum( unsigned int numberOfFields,
                        ROOT::Minuit2::FunctionMinimum const& minuitMinimum ) :
    variableValues( numberOfFields ),
    variableErrors( numberOfFields ),
    functionValue( minuitMinimum.Fval() ),
    functionError( minuitMinimum.Edm() )
  {
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      variableValues[ fieldIndex ]
      = minuitMinimum.UserParameters().Value( fieldIndex );
      variableErrors[ fieldIndex ]
      = minuitMinimum.UserParameters().Error( fieldIndex );
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
