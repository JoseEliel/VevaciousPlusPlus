/*
 * SlhaTrilinearDiagonalFunctionoid.hpp
 *
 *  Created on: Oct 30, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHATRILINEARDIAGONALFUNCTIONOID_HPP_
#define SLHATRILINEARDIAGONALFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"
#include "LagrangianParameterManagement/LhaSourcedParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaTrilinearDiagonalFunctionoid :
                                          public LhaSourcedParameterFunctionoid
  {
  public:
    SlhaTrilinearDiagonalFunctionoid( size_t const indexInValuesVector,
                                      size_t const directTrilinearIndex,
                                      size_t const trilinearOverYukawaIndex,
                                      size_t const appropriateYukawaIndex ) :
      LhaSourcedParameterFunctionoid( indexInValuesVector ),
      directTrilinearIndex( directTrilinearIndex ),
      trilinearOverYukawaIndex( trilinearOverYukawaIndex ),
      appropriateYukawaIndex( appropriateYukawaIndex ) {}
    virtual ~SlhaTrilinearDiagonalFunctionoid() {}

    // This returns the trilinear coupling evaluated at the scale given
    // through logarithmOfScale either from its direct value printed in the TL,
    // TE, TQ, TD, or TU blocks, or by multiplying the A-value with the
    // appropriate Yukawa coupling.
    virtual double operator()( double const logarithmOfScale,
                        std::vector< double > const& interpolatedValues ) const
    { return operator()( interpolatedValues[ directTrilinearIndex ],
                         interpolatedValues[ trilinearOverYukawaIndex ],
                         interpolatedValues[ appropriateYukawaIndex ] ); }

    // This returns the trilinear coupling evaluated at the scale given
    // through logarithmOfScale either from its direct value printed in the TL,
    // TE, TQ, TD, or TU blocks, or by multiplying the A-value with the
    // appropriate Yukawa coupling.
    double operator()( double const directTrilinear,
                       double const trilinearOverYukawa,
                       double const appropriateYukawa ) const
    { return ( ( directTrilinear != 0.0 ) ?
               directTrilinear :
               ( trilinearOverYukawa * appropriateYukawa ) ); }

    // This is for creating a Python version of the potential.
    virtual std::string
    PythonParameterEvaluation( int const indentationSpaces ) const;

    // This is mainly for debugging.
    virtual std::string AsDebuggingString() const;


  protected:
    size_t const directTrilinearIndex;
    size_t const trilinearOverYukawaIndex;
    size_t const appropriateYukawaIndex;
  };





  // This is for creating a Python version of the potential.
  inline std::string
  SlhaTrilinearDiagonalFunctionoid::PythonParameterEvaluation(
                                            int const indentationSpaces ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << std::setprecision( 12 )
    << PythonIndent( indentationSpaces )
    << "parameterValues[ " << IndexInValuesVector()
    << " ] = FirstIfNonzeroOtherwiseSecond( parameterValues[ "
    << directTrilinearIndex
    << " ], ( parameterValues[ " << trilinearOverYukawaIndex
    << " ] * parameterValues[ " << appropriateYukawaIndex << " ] ) )";
    return stringBuilder.str();
  }

  // This is mainly for debugging.
  inline std::string
  SlhaTrilinearDiagonalFunctionoid::AsDebuggingString() const
  {
    std::stringstream stringBuilder;
    stringBuilder
    << "IndexInValuesVector() = " << IndexInValuesVector() << std::endl;
    stringBuilder
    << "directTrilinearIndex = " << directTrilinearIndex << std::endl;
    stringBuilder << "trilinearOverYukawaIndex = "
    << trilinearOverYukawaIndex << std::endl;
    stringBuilder
    << "appropriateYukawaIndex = " << appropriateYukawaIndex << std::endl;
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHATRILINEARDIAGONALFUNCTIONOID_HPP_ */
