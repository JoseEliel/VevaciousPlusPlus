/*
 * PolynomialPathFromCoefficients.hpp
 *
 *  Created on: Jul 21, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALPATHFROMCOEFFICIENTS_HPP_
#define POLYNOMIALPATHFROMCOEFFICIENTS_HPP_

#include "CommonIncludes.hpp"
#include "TunnelPath.hpp"
#include "BasicFunctions/SimplePolynomial.hpp"

namespace VevaciousPlusPlus
{

  class PolynomialPathFromCoefficients : public TunnelPath
  {
  public:
    PolynomialPathFromCoefficients( size_t const numberOfFields,
                                    std::vector< double > const& falseVacuum,
                                 std::vector< double > const& vacuumDifference,
                             std::vector< double > const& pathParameterization,
                                    double const pathTemperature,
                                    size_t const referenceFieldIndex,
                                    size_t const coefficientsBeyondLinear );
    virtual ~PolynomialPathFromCoefficients();


    // This fills fieldConfiguration with the values that the fields should
    // have when the path auxiliary is given by auxiliaryValue.
    void PutOnPathAt( std::vector< double >& fieldConfiguration,
                      double const auxiliaryValue ) const;

    // This returns the dot product with itself of the derivative of the field
    // vector with respect to the path auxiliary evaluated at auxiliaryValue.
    double SlopeSquared( double const auxiliaryValue ) const;

    // This returns the dot product of the first derivative of the field vector
    // with the second derivative, both with respect to the path auxiliary,
    // evaluated at auxiliaryValue.
    double SlopeDotAcceleration( double const auxiliaryValue ) const;

    // This is for debugging.
    std::string AsDebuggingString() const;


  protected:
    std::vector< SimplePolynomial > fieldPolynomials;
    std::vector< SimplePolynomial > firstDerivatives;
    std::vector< SimplePolynomial > secondDerivatives;
  };




  // This fills fieldConfiguration with the values that the fields should
  // have when the path auxiliary is given by auxiliaryValue.
  inline void PolynomialPathFromCoefficients::PutOnPathAt(
                                     std::vector< double >& fieldConfiguration,
                                            double const auxiliaryValue ) const
  {
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      fieldConfiguration[ fieldIndex ]
      = fieldPolynomials[ fieldIndex ]( auxiliaryValue );
    }
  }

  // This returns the dot product with itself of the derivative of the field
  // vector with respect to the path auxiliary evaluated at auxiliaryValue.
  inline double PolynomialPathFromCoefficients::SlopeSquared(
                                            double const auxiliaryValue ) const
  {
    double returnValue( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      double const
      fieldSlope( firstDerivatives[ fieldIndex ]( auxiliaryValue ) );
      returnValue += ( fieldSlope * fieldSlope );
    }
    return returnValue;
  }


  // This returns the dot product of the first derivative of the field vector
  // with the second derivative, both with respect to the path auxiliary,
  // evaluated at auxiliaryValue.
  inline double PolynomialPathFromCoefficients::SlopeDotAcceleration(
                                            double const auxiliaryValue ) const
  {
    double returnValue( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      returnValue += ( firstDerivatives[ fieldIndex ]( auxiliaryValue )
                       * secondDerivatives[ fieldIndex ]( auxiliaryValue ) );
    }
    return returnValue;
  }

  // This is for debugging.
  inline std::string PolynomialPathFromCoefficients::AsDebuggingString() const
  {
    std::stringstream returnStream;
    returnStream << "{ ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        returnStream << ", " << std::endl;
      }
      returnStream << fieldPolynomials[ fieldIndex ].AsDebuggingString();
    }
    returnStream << " }";
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALPATHFROMCOEFFICIENTS_HPP_ */
