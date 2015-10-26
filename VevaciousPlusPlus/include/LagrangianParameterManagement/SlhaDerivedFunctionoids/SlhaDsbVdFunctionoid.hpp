/*
 * SlhaDsbVdFunctionoid.hpp
 *
 *  Created on: Oct 26, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHADSBVDFUNCTIONOID_HPP_
#define SLHADSBVDFUNCTIONOID_HPP_

namespace VevaciousPlusPlus
{

  class SlhaDsbVdFunctionoid : public SlhaDerivedFunctionoid
  {
  public:
    SlhaDsbVdFunctionoid();
    virtual ~SlhaDsbVdFunctionoid();

    // This returns the value of the DSB VEV for the neutral component of the
    // down-type Higgs double, looking in several sources in order.
    // 1) It tries HMIX[102], which is where SARAH-SPheno prints the value
    // directly. (It is not an official entry for HMIX though, and is not
    // present in SLHA files produced by any other spectrum generator.)
    // 2) If there is no non-zero value from HMIX[102], then HMIX[3] times the
    // cos(beta) is returned, looking for tan(beta) first in HMIX[2], then in
    // EXTPAR[25] if HMIX[2] was 0 (which the LHPC SLHA reader returns if there
    // was no explicit value for the entry in the block).
    virtual double operator()( double const logarithmOfScale ) const;

    // This is for creating a Python version of the potential.
    virtual std::string PythonParameterEvaluation() const = 0;


  protected:
    SlhaInterpolatedParameterFunctionoid const& sarahVd;
    SlhaInterpolatedParameterFunctionoid const& slhaV;
    SlhaInterpolatedParameterFunctionoid const& hmixTanBeta;
    SlhaInterpolatedParameterFunctionoid const& extparTanBeta;
  };





  // This returns the value of the DSB VEV for the neutral component of the
  // down-type Higgs double, looking in several sources in order.
  // 1) It tries HMIX[102], which is where SARAH-SPheno prints the value
  // directly. (It is not an official entry for HMIX though, and is not
  // present in SLHA files produced by any other spectrum generator.)
  // 2) If there is no non-zero value from HMIX[102], then HMIX[3] times the
  // cos(beta) is returned, looking for tan(beta) first in HMIX[2], then in
  // EXTPAR[25] if HMIX[2] was 0 (which the LHPC SLHA reader returns if there
  // was no explicit value for the entry in the block).
  inline double
  SlhaDsbVdFunctionoid::operator()( double const logarithmOfScale ) const
  {
    double returnValue( sarahVd( logarithmOfScale ) );
    if( ( returnValue > 0.0 )
        ||
        ( returnValue < 0.0 ) )
    {
      return returnValue;
    }
    double tanBeta( hmixTanBeta( logarithmOfScale ) );
    if( !( tanBeta > 0.0 )
        &&
        !( tanBeta < 0.0 ) )
    {
      tanBeta = extparTanBeta( logarithmOfScale );
    }
    return ( slhaV( logarithmOfScale ) * cos( atan( tanBeta ) ) );
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHADSBVDFUNCTIONOID_HPP_ */
