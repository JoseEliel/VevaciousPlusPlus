/*
 * SlhaTrilinearDiagonalFunctionoid.hpp
 *
 *  Created on: Oct 30, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHATRILINEARDIAGONALFUNCTIONOID_HPP_
#define SLHATRILINEARDIAGONALFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"
#include "LagrangianParameterManagement/SlhaDerivedFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaTrilinearDiagonalFunctionoid : public SlhaDerivedFunctionoid
  {
  public:
    SlhaTrilinearDiagonalFunctionoid(
                        SlhaSourcedParameterFunctionoid const& directTrilinear,
                    SlhaSourcedParameterFunctionoid const& trilinearOverYukawa,
                    SlhaSourcedParameterFunctionoid const& appropriateYukawa );
    virtual ~SlhaTrilinearDiagonalFunctionoid();


    // This returns trilinear coupling evaluated at the scale given
    // through logarithmOfScale either from its direct value printed in the TL,
    // TE, TQ, TD, or TU blocks, or by multiplying the A-value with the
    // appropriate Yukawa coupling.
    virtual double operator()( double const logarithmOfScale ) const;

    // This returns trilinear coupling evaluated at the scale given
    // through logarithmOfScale either from its direct value printed in the TL,
    // TE, TQ, TD, or TU blocks, or by multiplying the A-value with the
    // appropriate Yukawa coupling.
    virtual double operator()( double const logarithmOfScale,
                        std::vector< double > const& interpolatedValues ) const
    { return ( ( interpolatedValues[ directTrilinearIndex ] == 0.0 ) ?
                   ( interpolatedValues[ trilinearOverYukawaIndex ]
                     * interpolatedValues[ appropriateYukawaIndex ] ):
                   interpolatedValues[ directTrilinearIndex ] ); }

    // This is for creating a Python version of the potential.
    virtual std::string
    PythonParameterEvaluation( int const indentationSpaces ) const;


  protected:
    SlhaInterpolatedParameterFunctionoid const& directTrilinear;
    size_t const directTrilinearIndex;
    SlhaInterpolatedParameterFunctionoid const& trilinearOverYukawa;
    size_t const trilinearOverYukawaIndex;
    SlhaInterpolatedParameterFunctionoid const& appropriateYukawa;
    size_t const appropriateYukawaIndex;
  };





  // This returns trilinear coupling evaluated at the scale given
  // through logarithmOfScale either from its direct value printed in the TL,
  // TE, TQ, TD, or TU blocks, or by multiplying the A-value with the
  // appropriate Yukawa coupling.
  inline double SlhaTrilinearDiagonalFunctionoid::operator()(
                                          double const logarithmOfScale ) const
  {
    double directValue( directTrilinear( logarithmOfScale ) );
    if( directValue == 0.0 )
    {
      return ( trilinearOverYukawa( logarithmOfScale )
               * appropriateYukawa( logarithmOfScale ) );
    }
    else
    {
      return directValue;
    }
  }

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

} /* namespace VevaciousPlusPlus */

#endif /* SLHATRILINEARDIAGONALFUNCTIONOID_HPP_ */
