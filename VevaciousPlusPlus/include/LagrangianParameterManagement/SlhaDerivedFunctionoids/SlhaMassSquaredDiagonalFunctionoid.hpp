/*
 * SlhaMassSquaredDiagonalFunctionoid.hpp
 *
 *  Created on: Oct 30, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHAMASSSQUAREDDIAGONALFUNCTIONOID_HPP_
#define SLHAMASSSQUAREDDIAGONALFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"
#include "LagrangianParameterManagement/SlhaDerivedFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaMassSquaredDiagonalFunctionoid : public SlhaDerivedFunctionoid
  {
  public:
    SlhaMassSquaredDiagonalFunctionoid(
                             SlhaSourcedParameterFunctionoid const& squareMass,
                           SlhaSourcedParameterFunctionoid const& linearMass );
    virtual ~SlhaMassSquaredDiagonalFunctionoid();


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
    { return ( ( interpolatedValues[ squareMassIndex ] == 0.0 ) ?
                   ( interpolatedValues[ linearMassIndex ]
                     * interpolatedValues[ linearMassIndex ] ):
                   interpolatedValues[ squareMassIndex ] ); }

    // This is for creating a Python version of the potential.
    virtual std::string
    PythonParameterEvaluation( int const indentationSpaces ) const;


  protected:
    SlhaInterpolatedParameterFunctionoid const& squareMass;
    size_t const squareMassIndex;
    SlhaInterpolatedParameterFunctionoid const& linearMass;
    size_t const linearMassIndex;
  };





  // This returns trilinear coupling evaluated at the scale given
  // through logarithmOfScale either from its direct value printed in the TL,
  // TE, TQ, TD, or TU blocks, or by multiplying the A-value with the
  // appropriate Yukawa coupling.
  inline double SlhaMassSquaredDiagonalFunctionoid::operator()(
                                          double const logarithmOfScale ) const
  {
    double directValue( squareMass( logarithmOfScale ) );
    if( directValue == 0.0 )
    {
      double linearValue( linearMass( logarithmOfScale ) );
      return ( linearValue * linearValue );
    }
    else
    {
      return directValue;
    }
  }

  // This is for creating a Python version of the potential.
  inline std::string
  SlhaMassSquaredDiagonalFunctionoid::PythonParameterEvaluation(
                                            int const indentationSpaces ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << std::setprecision( 12 )
    << PythonIndent( indentationSpaces )
    << "parameterValues[ " << IndexInValuesVector()
    << " ] = FirstIfNonzeroOtherwiseSecond( parameterValues[ "
    << squareMassIndex
    << " ], ( ( parameterValues[ " << linearMassIndex << " ] )**2 ) )";
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHAMASSSQUAREDDIAGONALFUNCTIONOID_HPP_ */
