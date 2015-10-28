/*
 * SlhaDsbHiggsVevFunctionoid.hpp
 *
 *  Created on: Oct 26, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHADSBHIGGSVEVFUNCTIONOID_HPP_
#define SLHADSBHIGGSVEVFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"
#include "LagrangianParameterManagement/SlhaDerivedFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaDsbHiggsVevFunctionoid : public SlhaDerivedFunctionoid
  {
  public:
    SlhaDsbHiggsVevFunctionoid(
            SlhaInterpolatedParameterFunctionoid const& sarahHiggsVevComponent,
       SlhaInterpolatedParameterFunctionoid const& slhaHiggsVevEuclideanLength,
                                SlhaTwoSourceFunctionoid const& tanBeta,
                                bool const sinNotCos );
    virtual ~SlhaDsbHiggsVevFunctionoid();


    // This returns the value of the DSB VEV for the neutral component of the
    // functionoid's Higgs doublet. If HMIX has an entry for the component of
    // the VEV directly, its value is returned. (Entry 102 in HMIX is where
    // SARAH-SPheno directly prints the VEV for the doublet which gives masses
    // to the down-type quarks and 103 for that for the up-type quarks, though
    // they are not official entries for HMIX and thus is not present in SLHA
    // files produced by any other spectrum generator.) If there is no non-zero
    // value in HMIX for the component directly, then HMIX[3] times cos(beta)
    // or sin(beta) is returned depending on whether this functionoid is for
    // the doublet which gives mass to the down-type or up-type quarks
    // respectively, with beta calculated from the value of the functionoid
    // evaluating tan(beta).
    virtual double operator()( double const logarithmOfScale ) const;

    // This returns the value of the DSB VEV for the neutral component of the
    // functionoid's Higgs doublet. If HMIX has an entry for the component of
    // the VEV directly, its value is returned. (Entry 102 in HMIX is where
    // SARAH-SPheno directly prints the VEV for the doublet which gives masses
    // to the down-type quarks and 103 for that for the up-type quarks, though
    // they are not official entries for HMIX and thus is not present in SLHA
    // files produced by any other spectrum generator.) If there is no non-zero
    // value in HMIX for the component directly, then HMIX[3] times cos(beta)
    // or sin(beta) is returned depending on whether this functionoid is for
    // the doublet which gives mass to the down-type or up-type quarks
    // respectively, with beta calculated from the value of the functionoid
    // evaluating tan(beta).
    virtual double operator()( double const logarithmOfScale,
                       std::vector< double > const& interpolatedValues ) const;

    // This is for creating a Python version of the potential.
    virtual std::string
    PythonParameterEvaluation( int const indentationSpaces ) const;


  protected:
    SlhaInterpolatedParameterFunctionoid const& sarahHiggsVevComponent;
    size_t const sarahHmixIndex;
    SlhaInterpolatedParameterFunctionoid const& slhaHiggsVevEuclideanLength;
    size_t const slhaHmixIndex;
    SlhaTwoSourceFunctionoid const& tanBeta;
    size_t const tanBetaIndex;
    double (*cosOrSin)( double );
    std::string cosOrSinPythonString;
  };





  // This returns the value of the DSB VEV for the neutral component of the
  // functionoid's Higgs doublet. If HMIX has an entry for the component of
  // the VEV directly, its value is returned. (Entry 102 in HMIX is where
  // SARAH-SPheno directly prints the VEV for the doublet which gives masses
  // to the down-type quarks and 103 for that for the up-type quarks, though
  // they are not official entries for HMIX and thus is not present in SLHA
  // files produced by any other spectrum generator.) If there is no non-zero
  // value in HMIX for the component directly, then HMIX[3] times cos(beta)
  // or sin(beta) is returned depending on whether this functionoid is for
  // the doublet which gives mass to the down-type or up-type quarks
  // respectively, with beta calculated from the value of the functionoid
  // evaluating tan(beta).
  inline double SlhaDsbHiggsVevFunctionoid::operator()(
                                          double const logarithmOfScale ) const
  {
    double sarahValue( sarahHiggsVevComponent( logarithmOfScale ) );
    if( sarahValue != 0.0 )
    {
      return sarahValue;
    }
    else
    {
      return ( slhaHiggsVevEuclideanLength( logarithmOfScale )
               * (*cosOrSin)( atan( tanBeta( logarithmOfScale ) ) ) );
    }
  }

  // This returns the value of the DSB VEV for the neutral component of the
  // functionoid's Higgs doublet. If HMIX has an entry for the component of
  // the VEV directly, its value is returned. (Entry 102 in HMIX is where
  // SARAH-SPheno directly prints the VEV for the doublet which gives masses
  // to the down-type quarks and 103 for that for the up-type quarks, though
  // they are not official entries for HMIX and thus is not present in SLHA
  // files produced by any other spectrum generator.) If there is no non-zero
  // value in HMIX for the component directly, then HMIX[3] times cos(beta)
  // or sin(beta) is returned depending on whether this functionoid is for
  // the doublet which gives mass to the down-type or up-type quarks
  // respectively, with beta calculated from the value of the functionoid
  // evaluating tan(beta).
  inline double
  SlhaDsbHiggsVevFunctionoid::operator()( double const logarithmOfScale,
                        std::vector< double > const& interpolatedValues ) const
  {
    double sarahValue( interpolatedValues[ sarahHmixIndex ] );
    if( sarahValue != 0.0 )
    {
      return sarahValue;
    }
    else
    {
      return ( interpolatedValues[ slhaHmixIndex ]
               * (*cosOrSin)( atan( interpolatedValues[ tanBetaIndex ] ) ) );
    }
  }

  // This is for creating a Python version of the potential.
  inline std::string SlhaDsbHiggsVevFunctionoid::PythonParameterEvaluation(
                                            int const indentationSpaces ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << std::setprecision( 12 )
    << PythonIndent( indentationSpaces )
    << "SlhaDsbHiggsVevFunctionoid_sarahValue = parameterValues[ "
    << sarahHmixIndex() << " ]" << PythonNewlineThenIndent( indentationSpaces )
    << "if ( SlhaDsbHiggsVevFunctionoid_sarahValue != 0.0 ):"
    << PythonNewlineThenIndent( indentationSpaces + 2 )
    << "parameterValues[ " << IndexInValuesVector()
    << " ] = SlhaDsbHiggsVevFunctionoid_sarahValue"
    << PythonNewlineThenIndent( indentationSpaces )
    << "else:"
    << PythonNewlineThenIndent( indentationSpaces + 2 )
    << "parameterValues[ " << IndexInValuesVector()
    << " ] = " << cosOrSinPythonString << "( math.atan( parameterValues[ "
    << tanBetaIndex << " ] ) )";
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHADSBHIGGSVEVFUNCTIONOID_HPP_ */
