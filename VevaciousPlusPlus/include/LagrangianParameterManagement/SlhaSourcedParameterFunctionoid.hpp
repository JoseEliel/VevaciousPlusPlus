/*
 * SlhaSourcedParameterFunctionoid.hpp
 *
 *  Created on: Oct 28, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHASOURCEDPARAMETERFUNCTIONOID_HPP_
#define SLHASOURCEDPARAMETERFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class SlhaSourcedParameterFunctionoid
  {
  public:
    SlhaSourcedParameterFunctionoid( size_t const indexInValuesVector );
    virtual ~SlhaSourcedParameterFunctionoid();


    size_t IndexInValuesVector() const { return indexInValuesVector; }

    // This should return the value of the parameter at the requested scale
    // (exp( logarithmOfScale )) as an alternative to using a set of parameters
    // already evaluated at the scale. It is almost certain to be much slower
    // if used to obtain parameters repeatedly for different scales at the same
    // parameter point, but is more efficient for parameters which do not need
    // to be evaluated for the potential but still depend on Lagrangian
    // parameters, for example when evaluating the VEVs of the DSB vacuum for
    // the parameter point.
    virtual double operator()( double const logarithmOfScale ) const = 0;

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

#endif /* SLHASOURCEDPARAMETERFUNCTIONOID_HPP_ */
