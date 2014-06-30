/*
 * UnaryOperationFunctionoid.cpp
 *
 *  Created on: Mar 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/ParameterFunctionoids/UnaryOperationFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  UnaryOperationFunctionoid::UnaryOperationFunctionoid(
                                            double (*unaryOperation)( double ),
                                ParameterFunctionoid* const functionoidPointer,
                                             std::string const& creationString,
                                     std::string const& pythonParameterName ) :
    ParameterFunctionoid( creationString,
                          pythonParameterName ),
    unaryOperation( unaryOperation ),
    functionoidPointer( functionoidPointer )
  {
    // This constructor is just an initialization list.
  }

  UnaryOperationFunctionoid::~UnaryOperationFunctionoid()
  {
    // This does nothing.
  }


  // This is for creating a Python version of the potential.
  std::string UnaryOperationFunctionoid::PythonParameterEvaluation() const
  {
    std::stringstream stringBuilder;
    stringBuilder << std::setprecision( 12 ) << pythonParameterName << " = ";
    if( unaryOperation == &sqrt )
    {
      stringBuilder << "math.sqrt( ";
    }
    else if(  unaryOperation == &exp )
    {
      stringBuilder << "math.exp( ";
    }
    else if(  unaryOperation == &log )
    {
      stringBuilder << "math.log( ";
    }
    else if(  unaryOperation == &sin )
    {
      stringBuilder << "math.sin( ";
    }
    else if(  unaryOperation == &cos )
    {
      stringBuilder << "math.cos( ";
    }
    else if(  unaryOperation == &tan )
    {
      stringBuilder << "math.tan( ";
    }
    else if(  unaryOperation == &asin )
    {
      stringBuilder << "math.asin( ";
    }
    else if(  unaryOperation == &acos )
    {
      stringBuilder << "math.acos( ";
    }
    else if(  unaryOperation == &atan )
    {
      stringBuilder << "math.atan( ";
    }
    else if(  unaryOperation == &sinh )
    {
      stringBuilder << "math.sinh( ";
    }
    else if(  unaryOperation == &cosh )
    {
      stringBuilder << "math.cosh( ";
    }
    else if(  unaryOperation == &tanh )
    {
      stringBuilder << "math.tanh( ";
    }
    else
    {
      std::string errorMessage(
                  "Unknown unary operation for conversion to Python from \"" );
      errorMessage.append( creationString );
      errorMessage.append( "\"" );
      throw std::runtime_error( errorMessage );
    }
    stringBuilder << functionoidPointer->PythonParameterName() << " )";
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
