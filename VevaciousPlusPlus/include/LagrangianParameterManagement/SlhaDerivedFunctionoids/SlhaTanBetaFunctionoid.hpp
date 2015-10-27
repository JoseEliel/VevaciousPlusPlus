/*
 * SlhaTanBetaFunctionoid.hpp
 *
 *  Created on: Oct 27, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHATANBETAFUNCTIONOID_HPP_
#define SLHATANBETAFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"
#include "LagrangianParameterManagement/SlhaDerivedFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaTanBetaFunctionoid : public SlhaDerivedFunctionoid
  {
  public:
    SlhaTanBetaFunctionoid(
                       SlhaInterpolatedParameterFunctionoid const& hmixTanBeta,
                   SlhaInterpolatedParameterFunctionoid const& extparTanBeta );
    virtual ~SlhaTanBetaFunctionoid();


    // This returns the value of tan(beta) by first looking for HMIX[2], then
    // in EXTPAR[25] if HMIX[2] was 0 (which the LHPC SLHA reader returns if
    // there was no explicit value for the entry in the block).
    virtual double OnceOffValue( double const logarithmOfScale ) const
    { return FirstIfNonzeroOtherwiseSecond( hmixTanBeta,
                                            extparTanBeta,
                                            logarithmOfScale ); }

    // This returns the value of tan(beta) by first looking for HMIX[2], then
    // in EXTPAR[25] if HMIX[2] was 0 (which the LHPC SLHA reader returns if
    // there was no explicit value for the entry in the block).
    virtual double operator()( double const logarithmOfScale,
                       std::vector< double > const& interpolatedValues ) const
    { return
         FirstIfNonzeroOtherwiseSecond( interpolatedValues[ hmixTanBetaIndex ],
                                  interpolatedValues[ extparTanBetaIndex ] ); }

    // This is for creating a Python version of the potential.
    virtual std::string
    PythonParameterEvaluation( int const indentationSpaces ) const;


  protected:
    SlhaInterpolatedParameterFunctionoid const& hmixTanBeta;
    size_t const hmixTanBetaIndex;
    SlhaInterpolatedParameterFunctionoid const& extparTanBeta;
    size_t const extparTanBetaIndex;
  };





  // This is for creating a Python version of the potential.
  inline std::string SlhaTanBetaFunctionoid::PythonParameterEvaluation(
                                            int const indentationSpaces ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << PythonIndent( indentationSpaces )
    << "parameterValues[ " << IndexInValuesVector()
    << " ] = FirstIfNonzeroOtherwiseSecond( parameterValues[ "
    << hmixTanBetaIndex
    << " ], parameterValues[ " << extparTanBetaIndex << " ] )";
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHATANBETAFUNCTIONOID_HPP_ */
