/*
 * MinuitMinimum.hpp
 *
 *  Created on: Apr 11, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITMINIMUM_HPP_
#define MINUITMINIMUM_HPP_

#include "CommonIncludes.hpp"
#include "Minuit2/FunctionMinimum.h"

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
                   std::vector< double > const& variableErrors );
    MinuitMinimum();
    MinuitMinimum( MinuitMinimum const& copySource );
    virtual ~MinuitMinimum();


    std::vector< double >& VariableValues(){ return variableValues; }
    std::vector< double > const& VariableValues() const
    { return variableValues; }

    std::vector< double > const& VariableErrors() const
    { return variableErrors; }

    double FunctionValue() const{ return functionValue; }
    double FunctionError() const{ return functionError; }


    // This is mainly for debugging.
    std::string AsDebuggingString() const;


  protected:
    std::vector< double > variableValues;
    std::vector< double > variableErrors;
    double functionValue;
    double functionError;
  };




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
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINUITMINIMUM_HPP_ */
