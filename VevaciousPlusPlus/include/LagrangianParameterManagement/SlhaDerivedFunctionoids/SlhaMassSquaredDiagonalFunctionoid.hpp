/*
 * SlhaMassSquaredDiagonalFunctionoid.hpp
 *
 *  Created on: Oct 30, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHAMASSSQUAREDDIAGONALFUNCTIONOID_HPP_
#define SLHAMASSSQUAREDDIAGONALFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"
#include "LagrangianParameterManagement/SlhaSourcedParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaMassSquaredDiagonalFunctionoid :
                                         public SlhaSourcedParameterFunctionoid
  {
  public:
    SlhaMassSquaredDiagonalFunctionoid( size_t const indexInValuesVector,
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

    // This is mainly for debugging.
    virtual std::string AsDebuggingString() const;


  protected:
    SlhaSourcedParameterFunctionoid const& squareMass;
    size_t const squareMassIndex;
    SlhaSourcedParameterFunctionoid const& linearMass;
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

  // This is mainly for debugging.
  inline std::string
  SlhaMassSquaredDiagonalFunctionoid::AsDebuggingString() const
  {
    std::stringstream stringBuilder;
    stringBuilder
    << "IndexInValuesVector() = " << IndexInValuesVector() << std::endl;
    stringBuilder << "squareMassIndex = " << squareMassIndex << std::endl;
    stringBuilder
    << "squareMass = " << squareMass.AsDebuggingString() << std::endl;
    stringBuilder << "linearMassIndex = " << linearMassIndex << std::endl;
    stringBuilder
    << "linearMass = " << linearMass.AsDebuggingString() << std::endl;
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHAMASSSQUAREDDIAGONALFUNCTIONOID_HPP_ */
