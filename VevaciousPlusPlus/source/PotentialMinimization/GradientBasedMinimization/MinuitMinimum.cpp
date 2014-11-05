/*
 * MinuitMinimum.cpp
 *
 *  Created on: Apr 11, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialMinimization/GradientBasedMinimization/MinuitMinimum.hpp"

namespace VevaciousPlusPlus
{

  MinuitMinimum::MinuitMinimum( size_t numberOfVariables,
                        ROOT::Minuit2::FunctionMinimum const& minuitMinimum ) :
    variableValues( numberOfVariables ),
    variableErrors( numberOfVariables ),
    functionValue( minuitMinimum.Fval() ),
    functionError( minuitMinimum.Edm() ),
    isValidMinimum( minuitMinimum.IsValid() )
  {
    ROOT::Minuit2::MnUserParameters const&
    userParameters( minuitMinimum.UserParameters() );
    for( size_t variableIndex( 0 );
         variableIndex < numberOfVariables;
         ++variableIndex )
    {
      variableValues[ variableIndex ] = userParameters.Value( variableIndex );
      variableErrors[ variableIndex ] = userParameters.Error( variableIndex );
    }
  }

  MinuitMinimum::MinuitMinimum( std::vector< double > const& variableValues,
                                std::vector< double > const& variableErrors ) :
    variableValues( variableValues ),
    variableErrors( variableErrors ),
    functionValue( 0.0 ),
    functionError( -1.0 ),
    isValidMinimum( false )
  {
    // This constructor is just an initialization list.
  }

  MinuitMinimum::MinuitMinimum() :
    variableValues(),
    variableErrors(),
    functionValue( 0.0 ),
    functionError( -1.0 ),
    isValidMinimum( false )
  {
    // This constructor is just an initialization list.
  }

  MinuitMinimum::MinuitMinimum( MinuitMinimum const& copySource ) :
    variableValues( copySource.variableValues ),
    variableErrors( copySource.variableErrors ),
    functionValue( copySource.functionValue ),
    functionError( copySource.functionError ),
    isValidMinimum( copySource.isValidMinimum )
  {
    // This constructor is just an initialization list.
  }

  MinuitMinimum::~MinuitMinimum()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
