/*
 * MinuitMinimum.hpp
 *
 *  Created on: Apr 11, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITMINIMUM_HPP_
#define MINUITMINIMUM_HPP_

#include <cstddef>
#include "Minuit2/FunctionMinimum.h"
#include <vector>
#include <string>
#include <sstream>

namespace VevaciousPlusPlus
{
  // This is a class to convert Minuit2's clumsy ROOT::Minuit2::FunctionMinimum
  // instances into easy collections of doubles and vectors. Currently it just
  // ignores covariance and names.
  class MinuitMinimum
  {
  public:
    MinuitMinimum( size_t numberOfVariables,
                   ROOT::Minuit2::FunctionMinimum const& minuitMinimum );

    MinuitMinimum( std::vector< double > const& variableValues,
                   std::vector< double > const& variableErrors ) :
      variableValues( variableValues ),
      variableErrors( variableErrors ),
      functionValue( 0.0 ),
      functionError( -1.0 ),
      isValidMinimum( false ) {}

    MinuitMinimum( std::vector< double > const& variableValues,
                   double const functionValue ) :
      variableValues( variableValues ),
      variableErrors( variableValues.size(),
                      0.0 ),
      functionValue( functionValue ),
      functionError( 0.0 ),
      isValidMinimum( false ) {}

    MinuitMinimum( std::vector< double > const& variableValues,
                   std::vector< double > const& variableErrors,
                   double const functionValue,
                   double const functionError ) :
      variableValues( variableValues ),
      variableErrors( variableErrors ),
      functionValue( functionValue ),
      functionError( functionError ),
      isValidMinimum( false ) {}

    MinuitMinimum() : variableValues(),
                      variableErrors(),
                      functionValue( 0.0 ),
                      functionError( -1.0 ),
                      isValidMinimum( false ) {}

    MinuitMinimum( MinuitMinimum const& copySource ) :
      variableValues( copySource.variableValues ),
      variableErrors( copySource.variableErrors ),
      functionValue( copySource.functionValue ),
      functionError( copySource.functionError ),
      isValidMinimum( copySource.isValidMinimum ) {}

    virtual ~MinuitMinimum() {}


    std::vector< double >& VariableValues() { return variableValues; }

    std::vector< double > const& VariableValues() const
    { return variableValues; }

    std::vector< double > const& VariableErrors() const
    { return variableErrors; }

    double FunctionValue() const { return functionValue; }

    double FunctionError() const { return functionError; }

    bool IsValidMinimum() const { return isValidMinimum; }


    // This is mainly for debugging.
    std::string AsDebuggingString() const;


  protected:
    std::vector< double > variableValues;
    std::vector< double > variableErrors;
    double functionValue;
    double functionError;
    bool isValidMinimum;
  };





  inline MinuitMinimum::MinuitMinimum( size_t numberOfVariables,
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

  // This is mainly for debugging.
  inline std::string MinuitMinimum::AsDebuggingString() const
  {
    std::stringstream stringBuilder;
    stringBuilder << "functionValue = " << functionValue << std::endl
    << "functionError = " << functionError << std::endl;
    for( size_t variableIndex( 0 );
         variableIndex < variableValues.size();
         ++variableIndex )
    {
      stringBuilder << "variableValues[ " << variableIndex << " ] = "
      << variableValues[ variableIndex ] << ", variableErrors[ "
      << variableIndex << " ] = " << variableErrors[ variableIndex ]
      << std::endl;
    }
    stringBuilder << "isValidMinimum = " << isValidMinimum << std::endl;
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINUITMINIMUM_HPP_ */
