/*
 * SlhaDerivedFunctionoid.hpp
 *
 *  Created on: Oct 23, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHADERIVEDFUNCTIONOID_HPP_
#define SLHADERIVEDFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class SlhaDerivedFunctionoid
  {
  public:
    SlhaDerivedFunctionoid();
    virtual ~SlhaDerivedFunctionoid();


    size_t IndexInValuesVector() const { return indexInValuesVector; }

    // This should return the value of the parameter at the requested scale
    // (exp( logarithmOfScale )) as an alternative to using a set of parameters
    // already evaluated at the scale. It is almost certain to be much slower
    // if used to obtain parameters repeatedly for different scales at the same
    // parameter point, but is more efficient for parameters which do not need
    // to be evaluated for the potential but still depend on Lagrangian
    // parameters, for example when evaluating the VEVs of the DSB vacuum for
    // the parameter point.
    virtual double OnceOffValue( double const logarithmOfScale ) const = 0;

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

    // This is false until the functionoid has been assigned an index as an
    // active Lagrangian parameter.
    bool IsActive() const { return isActive; }

    // This sets indexInValuesVector to activeIndex and sets isActive to true.
    void MakeActive( size_t const activeIndex );


  protected:
    // This returns firstValue if it is non-zero, otherwise secondValue.
    static double FirstIfNonzeroOtherwiseSecond( double const firstValue,
                                                 double const secondValue )
    { return ( ( firstValue == 0.0 ) ? secondValue : firstValue ); }

    // This evaluates firstFunctionoid( logarithmOfScale ) and if it is
    // non-zero returns that value. Otherwise it returns
    // secondFunctionoid( logarithmOfScale ).
    static double FirstIfNonzeroOtherwiseSecond(
                  SlhaInterpolatedParameterFunctionoid const& firstFunctionoid,
                 SlhaInterpolatedParameterFunctionoid const& secondFunctionoid,
                                               double const logarithmOfScale );

    static std::string PythonIndent( int const indentationSpaces )
    { return std::string( indentationSpaces,
                          ' ' ); }

    static std::string PythonNewlineThenIndent( int const indentationSpaces )
    { return ( std::string( "\n" ) + std::string( indentationSpaces,
                                                  ' ' ) ); }

    size_t indexInValuesVector;
    bool isActive;
  };





  // This sets indexInValuesVector to activeIndex and sets isActive to true.
  inline void SlhaDerivedFunctionoid::MakeActive( size_t const activeIndex )
  {
    indexInValuesVector = activeIndex;
    isActive = true;
  }

  // This evaluates firstFunctionoid( logarithmOfScale ) and if it is
  // non-zero returns that value. Otherwise it returns
  // secondFunctionoid( logarithmOfScale ).
  inline double
  SlhaDerivedFunctionoid::FirstIfNonzeroOtherwiseSecond(
                  SlhaInterpolatedParameterFunctionoid const& firstFunctionoid,
                 SlhaInterpolatedParameterFunctionoid const& secondFunctionoid,
                                                 double const evaluationScale )
  {
    double firstFunctionoidValue( firstFunctionoid( evaluationScale ) );
    if( firstFunctionoidValue == 0.0 )
    {
      return secondFunctionoid( evaluationScale );
    }
    else
    {
      return firstFunctionoidValue;
    }
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHADERIVEDFUNCTIONOID_HPP_ */
