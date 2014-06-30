/*
 * BinaryOperationFunctionoid.cpp
 *
 *  Created on: Mar 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/ParameterFunctionoids/BinaryOperationFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  BinaryOperationFunctionoid::BinaryOperationFunctionoid(
                                       double (*binaryOperation)( double const,
                                                                double const ),
                                  ParameterFunctionoid* const firstFunctionoid,
                                 ParameterFunctionoid* const secondFunctionoid,
                                             std::string const& creationString,
                                     std::string const& pythonParameterName ) :
    ParameterFunctionoid( creationString,
                          pythonParameterName ),
    binaryOperation( binaryOperation ),
    firstFunctionoid( firstFunctionoid ),
    secondFunctionoid( secondFunctionoid )
  {
    // This constructor is just an initialization list.
  }

  BinaryOperationFunctionoid::~BinaryOperationFunctionoid()
  {
    // This does nothing.
  }


  // This is for creating a Python version of the potential.
  std::string BinaryOperationFunctionoid::PythonParameterEvaluation() const
  {
    std::stringstream stringBuilder;
    stringBuilder << std::setprecision( 12 ) << pythonParameterName << " = ( "
    << firstFunctionoid->PythonParameterName();
    if( binaryOperation == &PlusFunction )
    {
      stringBuilder << " + ";
    }
    else if( binaryOperation == &MinusFunction )
    {
      stringBuilder << " - ";
    }
    else if( binaryOperation == &MultiplyFunction )
    {
      stringBuilder << " * ";
    }
    else if( binaryOperation == &DivideFunction )
    {
      stringBuilder << " / ";
    }
    else if( binaryOperation == &IfNonZeroFunction )
    {
      stringBuilder << " if ( " << firstFunctionoid->PythonParameterName()
      << " != 0.0 ) else ";
    }
    else if( binaryOperation == &pow )
    {
      stringBuilder << "**";
    }
    else
    {
      std::string errorMessage(
                 "Unknown binary operation for conversion to Python from \"" );
      errorMessage.append( creationString );
      errorMessage.append( "\"" );
      throw std::runtime_error( errorMessage );
    }
    stringBuilder << secondFunctionoid->PythonParameterName() << " )";
    return stringBuilder.str();
  }
} /* namespace VevaciousPlusPlus */
