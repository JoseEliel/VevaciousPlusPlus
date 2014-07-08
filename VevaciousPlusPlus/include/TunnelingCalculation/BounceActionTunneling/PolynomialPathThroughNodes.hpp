/*
 * PolynomialPathThroughNodes.hpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALPATHTHROUGHNODES_HPP_
#define POLYNOMIALPATHTHROUGHNODES_HPP_

#include "CommonIncludes.hpp"
#include "TunnelPath.hpp"
#include "Eigen/Dense"
#include "BasicFunctions/SimplePolynomial.hpp"

namespace VevaciousPlusPlus
{

  class PolynomialPathThroughNodes : public TunnelPath
  {
  public:
    PolynomialPathThroughNodes(
                       std::vector< std::vector< double > > const& pathNodes,
                                double const pathTemperature,
                                Eigen::MatrixXd const& pathStepsInverse );
    virtual ~PolynomialPathThroughNodes();


    // This fills fieldConfiguration with the values that the fields
    // should have when the path auxiliary is given by auxiliaryValue.
    void PutOnPathAt( std::vector< double >& fieldConfiguration,
                      double const auxiliaryValue ) const;

    // This returns the dot product with itself of the derivative of the
    // field vector with respect to the path auxiliary evaluated at
    // auxiliaryValue.
    double SlopeSquared( double const auxiliaryValue ) const;

    // This returns the dot product of the first derivative of the field
    // vector with the second derivative, both with respect to the path
    // auxiliary, evaluated at auxiliaryValue.
    double SlopeDotAcceleration( double const auxiliaryValue ) const;

    // This is for debugging.
    std::string AsDebuggingString() const;


  protected:
    std::vector< SimplePolynomial > fieldPolynomials;
    std::vector< SimplePolynomial > firstDerivatives;
    std::vector< SimplePolynomial > secondDerivatives;
  };




  // This fills fieldConfiguration with the values that the fields
  // should have when the path auxiliary is given by auxiliaryValue.
  inline void PolynomialPathThroughNodes::PutOnPathAt(
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

  // This returns the dot product with itself of the derivative of the
  // field vector with respect to the path auxiliary evaluated at
  // auxiliaryValue.
  inline double PolynomialPathThroughNodes::SlopeSquared(
                                            double const auxiliaryValue ) const
  {
    double returnValue( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      double const
      slopeComponent( firstDerivatives[ fieldIndex ]( auxiliaryValue ) );
      returnValue += ( slopeComponent * slopeComponent );
    }
    return returnValue;
  }

  // This returns the dot product of the first derivative of the field
  // vector with the second derivative, both with respect to the path
  // auxiliary, evaluated at auxiliaryValue.
  inline double PolynomialPathThroughNodes::SlopeDotAcceleration(
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

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALPATHTHROUGHNODES_HPP_ */
