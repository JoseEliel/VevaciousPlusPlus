/*
 * SlhaDsbHiggsVevFunctionoid.hpp
 *
 *  Created on: Oct 26, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHADSBHIGGSVEVFUNCTIONOID_HPP_
#define SLHADSBHIGGSVEVFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"
#include "LagrangianParameterManagement/SlhaSourcedParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaDsbHiggsVevFunctionoid : public SlhaSourcedParameterFunctionoid
  {
  public:
    SlhaDsbHiggsVevFunctionoid( size_t const indexInValuesVector,
                     SlhaSourcedParameterFunctionoid const& vevEuclideanLength,
                                SlhaSourcedParameterFunctionoid const& tanBeta,
                                bool const sinNotCos );
    virtual ~SlhaDsbHiggsVevFunctionoid();


    // This returns the value of the DSB VEV for the neutral component of the
    // functionoid's Higgs doublet: HMIX[3] times either cos(beta) or sin(beta)
    // is returned depending on whether this functionoid is for the doublet
    // which gives mass to the down-type or up-type quarks respectively, with
    // beta calculated from the value of the functionoid evaluating tan(beta).
    virtual double operator()( double const logarithmOfScale ) const
    { return ( vevEuclideanLength( logarithmOfScale )
               * (*cosOrSin)( atan( tanBeta( logarithmOfScale ) ) ) ); }

    // This returns the value of the DSB VEV for the neutral component of the
    // functionoid's Higgs doublet: HMIX[3] times either cos(beta) or sin(beta)
    // is returned depending on whether this functionoid is for the doublet
    // which gives mass to the down-type or up-type quarks respectively, with
    // beta calculated from the value of the functionoid evaluating tan(beta).
    virtual double operator()( double const logarithmOfScale,
                       std::vector< double > const& interpolatedValues ) const
    {return ( interpolatedValues[ vevIndex ]
              * (*cosOrSin)( atan( interpolatedValues[ tanBetaIndex ] ) ) ); }

    // This is for creating a Python version of the potential.
    virtual std::string
    PythonParameterEvaluation( int const indentationSpaces ) const;


  protected:
    SlhaSourcedParameterFunctionoid const& vevEuclideanLength;
    size_t const vevIndex;
    SlhaSourcedParameterFunctionoid const& tanBeta;
    size_t const tanBetaIndex;
    double (*cosOrSin)( double );
    std::string cosOrSinPythonString;
  };





  // This is for creating a Python version of the potential.
  inline std::string SlhaDsbHiggsVevFunctionoid::PythonParameterEvaluation(
                                            int const indentationSpaces ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << std::setprecision( 12 )
    << PythonIndent( indentationSpaces )
    << "parameterValues[ " << IndexInValuesVector()
    << " ] = " << cosOrSinPythonString << "( math.atan( parameterValues[ "
    << tanBetaIndex << " ] ) )";
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHADSBHIGGSVEVFUNCTIONOID_HPP_ */
