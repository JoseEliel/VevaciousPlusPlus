/*
 * LhaInterpolatedParameterFunctionoid.hpp
 *
 *  Created on: Oct 23, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef LHAINTERPOLATEDPARAMETERFUNCTIONOID_HPP_
#define LHAINTERPOLATEDPARAMETERFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"
#include "LhaSourcedParameterFunctionoid.hpp"
#include "LHPC/SimpleLhaParser.hpp"

namespace VevaciousPlusPlus
{

  class LhaInterpolatedParameterFunctionoid :
                                          public LhaSourcedParameterFunctionoid
  {
  public:
    LhaInterpolatedParameterFunctionoid( size_t const indexInValuesVector,
                                        LHPC::SimpleLhaParser const& lhaParser,
                                          std::string const& parameterName )  :
      LhaSourcedParameterFunctionoid( indexInValuesVector ),
      parameterName( parameterName ),
      lhaParser( &lhaParser ) {}

    LhaInterpolatedParameterFunctionoid(
                    LhaInterpolatedParameterFunctionoid const& copySource )  :
      LhaSourcedParameterFunctionoid( copySource.indexInValuesVector ),
      parameterName( copySource.parameterName ),
      lhaParser( copySource.lhaParser ) {}

    virtual ~LhaInterpolatedParameterFunctionoid() {}


    // This should return the value of the functionoid for the given logarithm
    // of the scale.
    virtual double operator()( double const logarithmOfScale ) const = 0;

    // This should return the value of the functionoid for the given logarithm
    // of the scale. It should ignore the values of the other parameters.
    virtual double operator()( double const logarithmOfScale,
                   std::vector< double > const& interpolatedValues ) const = 0;

    // This should re-calculate the coefficients of the polynomial of the
    // logarithm of the scale used in evaluating the functionoid.
    virtual void UpdateForNewLhaParameters() = 0;

    // This is for creating a Python version of the potential.
    virtual std::string
    PythonParameterEvaluation( int const indentationSpaces ) const = 0;


  protected:
    std::string parameterName;
    LHPC::SimpleLhaParser const* lhaParser;
  };

} /* namespace VevaciousPlusPlus */

#endif /* LHAINTERPOLATEDPARAMETERFUNCTIONOID_HPP_ */
