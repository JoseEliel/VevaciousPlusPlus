/*
 * ModifiedBounceForMinuit.hpp
 *
 *  Created on: May 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MODIFIEDBOUNCEFORMINUIT_HPP_
#define MODIFIEDBOUNCEFORMINUIT_HPP_

#include "../CommonIncludes.hpp"
#include "BounceActionForMinuit.hpp"

namespace VevaciousPlusPlus
{

  class ModifiedBounceForMinuit : public BounceActionForMinuit
  {
  public:
    ModifiedBounceForMinuit( PotentialFunction const& potentialFunction,
                             size_t const numberOfVaryingPathNodes,
                             size_t const numberOfSplinesInPotential,
                             PotentialMinimum const& falseVacuum,
                             PotentialMinimum const& trueVacuum,
                             double const falseVacuumEvaporationTemperature,
                             size_t const undershootOvershootAttempts = 32,
                             size_t const maximumMultipleOfLongestLength = 16,
                           double const initialFractionOfShortestLength = 0.05,
                             double const minimumScaleSquared = 1.0,
                            double const shootingCloseEnoughThreshold = 0.01 );
    virtual
    ~ModifiedBounceForMinuit();


    // This modifies the return value if the temperature is non-zero to give
    // the dimensionless value ( [ S_3(T)/T ] + log[ S_3(T)/GeV ] ), which is
    // the quantity that should be minimized by MINUIT to minimize the DSB
    // vacuum survival probability under the approximations made in
    // arXiv:1405.7376, when we want to minimize the number of minimizations of
    // the bounce action.
    virtual double
    operator()( std::vector< double > const& pathParameterization ) const;
  };


  // This modifies the return value if the temperature is non-zero to give
  // the dimensionless value ( [ S_3(T)/T ] + log[ S_3(T)/GeV ] ), which is
  // the quantity that should be minimized by MINUIT to minimize the DSB
  // vacuum survival probability under the approximations made in
  // arXiv:1405.7376, when we want to minimize the number of minimizations of
  // the bounce action.
  inline double ModifiedBounceForMinuit::operator()(
                      std::vector< double > const& pathParameterization ) const
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "ModifiedBounceForMinuit::operator()( { ";
    for( std::vector< double >::const_iterator
         pathParameter( pathParameterization.begin() );
         pathParameter < pathParameterization.end();
         ++pathParameter )
    {
      if( pathParameter > pathParameterization.begin() )
      {
        std::cout << ", ";
      }
      std::cout << *pathParameter;
    }
    std::cout << " } ): ";/**/
    PathFieldsAndPotential
    pathFieldsAndPotential( DecodePathParameters( pathParameterization ) );
    PotentialAlongPath( pathFieldsAndPotential );
    if( pathFieldsAndPotential.NonZeroTemperature() )
    {
      double const bounceAction( BounceAction( pathFieldsAndPotential ) );
      return ( ( bounceAction / pathFieldsAndPotential.GivenTemperature() )
               + log( bounceAction ) );
    }
    else
    {
      // debugging:
      /**/double returnValue( BounceAction( pathFieldsAndPotential ) );
      std::cout << " bounce action = " << returnValue;
      return returnValue;/**/
      // return BounceAction( pathFieldsAndPotential );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* MODIFIEDBOUNCEFORMINUIT_HPP_ */
