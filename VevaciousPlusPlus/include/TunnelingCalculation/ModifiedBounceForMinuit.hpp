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
#include "boost/numeric/odeint/integrate/integrate.hpp"
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
    size_t const numberOfFields;
    size_t referenceFieldIndex;
    size_t const numberOfSplineFields;
    size_t potentialApproximationPower;
    PotentialMinimum const& falseVacuum;
    PotentialMinimum const& trueVacuum;
    double const fieldOriginPotential;
    double const falseVacuumPotential;
    double const trueVacuumPotential;
    double const falseVacuumEvaporationTemperature;
    double tunnelingScaleSquared;
    double shortestLength;
    double longestLength;

    // This turns a flattened matrix of coefficients from splineCoefficients
    // and fills fieldsAsPolynomials appropriately. The coefficients are taken to
    // be in the order
    // [ c_{1,0}, c_{1,1}, ..., c_{1, (referenceFieldIndex-1)},
    //           c_{1, (referenceFieldIndex+1)}, ..., c_{1,(numberOfFields-1)},
    //   c_{2,0}, c_{2,1}, ..., c_{2, (referenceFieldIndex-1)},
    //           c_{2, (referenceFieldIndex+1)}, ..., c_{1,(numberOfFields-1)},
    //   ...
    //   c_{p,0}, c_{p,1}, ..., c_{p, (referenceFieldIndex-1)},
    //           c_{p, (referenceFieldIndex+1)}, ..., c_{p,(numberOfFields-1)},
    //   temperature ],
    // where given field [j] is then the sum of c_{i,j} * a^i and p is the
    // greatest power given implicitly by splineCoefficients. Note that the
    // given fields do not map completely to the fields of potentialFunction:
    // field [referenceFieldIndex] is skipped over by splineCoefficients, as it
    // is set to be linear in a going from the false vacuum to the true vacuum.
    // It also puts the value of the potential (minus the value at the field
    // origin at zero temperature) into thermalFalseVacuumPotential and
    // thermalTrueVacuumPotential for the false and true vacua at the
    // temperature given by splineCoefficients.back() if and only if it is
    // non-zero. Unfortunately it does not use "return value optimization" as
    // we cannot be sure that the user will compile with C++11 features
    // enabled.
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
