/*
 * ModifiedBounceForMinuit.hpp
 *
 *  Created on: May 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MODIFIEDBOUNCEFORMINUIT_HPP_
#define MODIFIEDBOUNCEFORMINUIT_HPP_

#include "../CommonIncludes.hpp"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "Eigen/Dense"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "boost/numeric/odeint/integrate/integrate_adaptive.hpp"
#include "../PotentialEvaluation.hpp"
#include "BubbleProfiler.hpp"

namespace VevaciousPlusPlus
{

  class ModifiedBounceForMinuit : public ROOT::Minuit2::FCNBase
  {
  public:
    ModifiedBounceForMinuit( PotentialFunction const& potentialFunction,
                             unsigned int const potentialApproximationPower,
                             PotentialMinimum const& falseVacuum,
                             PotentialMinimum const& trueVacuum,
                             double const dsbEvaporationTemperature );
    virtual
    ~ModifiedBounceForMinuit();


    // This implements operator() for FCNBase, the function that MINUIT will
    // minimize. The values of splineCoefficients should be sets of n
    // coefficients for polynomials for each of the n fields, plus the final
    // element of splineCoefficients should be the temperature.
    virtual double
    operator()( std::vector< double > const& splineCoefficients ) const;

    // This implements Up() for FCNBase just to stick to a basic value.
    virtual double Up() const { return 1.0; }


  protected:
    PotentialFunction const& potentialFunction;
    unsigned int numberOfFields;
    unsigned int potentialApproximationPower;
    PotentialMinimum const& falseVacuum;
    PotentialMinimum const& trueVacuum;
    double const fieldOriginPotential;
    double const falseVacuumPotential;
    double const trueVacuumPotential;
    double const falseVacuumEvaporationTemperature;

    // This turns a flattened matrix of coefficients from splineCoefficients
    // and fills fieldsAsPolynomials appropriately. Unfortunately it does not
    // use "return value optimization" as we cannot be sure that the user will
    // compile with C++11 features enabled.
    void DecodeSplineVector( std::vector< double > const& splineCoefficients,
                          std::vector< SimplePolynomial >& fieldsAsPolynomials,
                             double& thermalFalseVacuumPotential,
                             double& thermalTrueVacuumPotential ) const;

    // This fills fieldConfiguration with the values from fieldsAsPolynomials
    // at the auxiliary variable = auxiliaryValue.
    void SetFieldConfiguration( std::vector< double >& fieldConfiguration,
                    std::vector< SimplePolynomial > const& fieldsAsPolynomials,
                                double const auxiliaryValue ) const;

    // This returns a polynomial approximation of the potential along the path
    // given by splineCoefficients.
    SimplePolynomial PotentialAlongPath(
                    std::vector< SimplePolynomial > const& fieldsAsPolynomials,
                                         double const falseVacuumPotential,
                                      double const trueVacuumPotential ) const;

    // This puts the bounce action in the thin-wall approximation into
    // bounceAction and returns true if the thin-wall approximation is
    // appropriate. Otherwise, bounceAction is left alone and false is
    // returned.
    bool ThinWallAppropriate( double potentialDifference,
                              double const tunnelingScale,
                              double const givenTemperature,
                       std::vector< SimplePolynomial > const& fieldDerivatives,
                              SimplePolynomial const& potentialApproximation,
                              double& bounceAction ) const;
  };




  // This fills fieldConfiguration with the values from fieldsAsPolynomials
  // at the auxiliary variable = auxiliaryValue.
  inline void ModifiedBounceForMinuit::SetFieldConfiguration(
                                     std::vector< double >& fieldConfiguration,
                    std::vector< SimplePolynomial > const& fieldsAsPolynomials,
                                            double const auxiliaryValue ) const
  {
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      fieldConfiguration[ fieldIndex ]
      = fieldsAsPolynomials[ fieldIndex ]( auxiliaryValue );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* MODIFIEDBOUNCEFORMINUIT_HPP_ */
