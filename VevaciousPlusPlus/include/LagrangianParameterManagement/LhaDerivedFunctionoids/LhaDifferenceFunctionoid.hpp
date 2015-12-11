/*
 * LhaDifferenceFunctionoid.hpp
 *
 *  Created on: Nov 27, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef LHADIFFERENCEFUNCTIONOID_HPP_
#define LHADIFFERENCEFUNCTIONOID_HPP_

#include "LagrangianParameterManagement/LhaSourcedParameterFunctionoid.hpp"
#include <vector>
#include <string>
#include <sstream>

namespace VevaciousPlusPlus
{

  class LhaDifferenceFunctionoid : public LhaSourcedParameterFunctionoid
  {
  public:
    LhaDifferenceFunctionoid( size_t const indexInValuesVector,
                              size_t const subtractorIndex,
                              size_t const subtractedIndex ) :
      LhaSourcedParameterFunctionoid( indexInValuesVector ),
      subtractorIndex( subtractorIndex ),
      subtractedIndex( subtractedIndex ) {}

    virtual ~LhaDifferenceFunctionoid() {}


    // This returns the parameter value at firstChoiceIndex if it is non-zero,
    // otherwise it returns the parameter value at secondChoiceIndex.
    virtual double operator()( double const logarithmOfScale,
                        std::vector< double > const& interpolatedValues ) const
    { return operator()( interpolatedValues[ subtractorIndex ],
                         interpolatedValues[ subtractedIndex ] ); }

    // This returns the first parameter value if it is non-zero, otherwise the
    // second parameter value.
    double operator()( double const subtractorValue,
                       double const subtractedValue ) const
    { return ( subtractorValue - subtractedValue ); }

    // This is for creating a Python version of the potential.
    virtual std::string
    PythonParameterEvaluation( int const indentationSpaces ) const;

    // This is mainly for debugging.
    virtual std::string AsDebuggingString() const;


  protected:
    size_t const subtractorIndex;
    size_t const subtractedIndex;
  };





  // This is for creating a Python version of the potential.
  inline std::string LhaDifferenceFunctionoid::PythonParameterEvaluation(
                                            int const indentationSpaces ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << std::setprecision( 12 )
    << PythonIndent( indentationSpaces )
    << "parameterValues[ " << IndexInValuesVector()
    << " ] = ( parameterValues[ " << subtractorIndex
    << " ] - parameterValues[ " << subtractedIndex << " ] )";
    return stringBuilder.str();
  }

  // This is mainly for debugging.
  inline std::string LhaDifferenceFunctionoid::AsDebuggingString() const
  {
    std::stringstream stringBuilder;
    stringBuilder
    << "IndexInValuesVector() = " << IndexInValuesVector() << std::endl;
    stringBuilder << "subtractorIndex = " << subtractorIndex << std::endl;
    stringBuilder << "subtractedIndex = " << subtractedIndex << std::endl;
    return stringBuilder.str();
  }

}

#endif /* LHADIFFERENCEFUNCTIONOID_HPP_ */
