/*
 * LhaSourcedParameterFunctionoid.hpp
 *
 *  Created on: Oct 28, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef LHASOURCEDPARAMETERFUNCTIONOID_HPP_
#define LHASOURCEDPARAMETERFUNCTIONOID_HPP_

#include <cstddef>
#include <vector>
#include <string>

namespace VevaciousPlusPlus
{

  class LhaSourcedParameterFunctionoid
  {
  public:
    LhaSourcedParameterFunctionoid( size_t const indexInValuesVector ) :
      indexInValuesVector( indexInValuesVector ) {}

    virtual ~LhaSourcedParameterFunctionoid() {}


    size_t IndexInValuesVector() const { return indexInValuesVector; }

    // This should return the value of the functionoid for the given logarithm
    // of the scale, using the values of the parameters directly interpolated
    // from the values explicitly given in the SLHA file, given by
    // interpolatedValues (which should have correct values in the elements
    // with index lower than indexInValuesVector).
    virtual double operator()( double const logarithmOfScale,
                   std::vector< double > const& interpolatedValues ) const = 0;

    // This should return a string for creating a Python version of the
    // potential, indented by indentationSpaces spaces.
    virtual std::string
    PythonParameterEvaluation( int const indentationSpaces ) const = 0;

    // This is mainly for debugging.
    virtual std::string AsDebuggingString() const = 0;


  protected:
    static std::string PythonIndent( int const indentationSpaces )
    { return std::string( indentationSpaces,
                          ' ' ); }

    static std::string PythonNewlineThenIndent( int const indentationSpaces )
    { return ( std::string( "\n" ) + std::string( indentationSpaces,
                                                  ' ' ) ); }

    size_t indexInValuesVector;
  };

} /* namespace VevaciousPlusPlus */

#endif /* LHASOURCEDPARAMETERFUNCTIONOID_HPP_ */
