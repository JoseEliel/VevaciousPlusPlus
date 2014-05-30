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
                             double const dsbEvaporationTemperature,
                             // size_t const pathResolution = 128,
                            unsigned int const undershootOvershootAttempts = 8,
                        unsigned int const maximumMultipleOfLongestLength = 16,
                           double const initialFractionOfShortestLength = 0.05,
                     unsigned int const energyConservingUndershootAttempts = 4,
                            double const minimumScaleSquared = 1.0,
                            double const shootingCloseEnoughThreshold = 0.01 );
    virtual
    ~ModifiedBounceForMinuit();


    // The bounce action is calculated by the following process:
    // 1) The path in field space from the false vacuum to the true vacuum is
    //    decoded from pathParameterization to give the field configuration f
    //    as a function of a path auxiliary variable p, giving f(p) and df/dp.
    //    See the comment above DecodePathParameters for the format of
    //    pathParameterization.
    // 2) The potential is fitted as a polynomial in p of degree
    //    potentialApproximationPower, giving V(p).
    // 3) If the potential difference between the vacua is small enough, the
    //    thin-wall approximation is tried, but only returned if the resulting
    //    bubble radius turns out to be large enough compared to the wall
    //    thickness.
    // 4) If the thin-wall approximation is not returned, the one-dimensional
    //    bubble equation of motion along p is integrated to get p as a
    //    function of the radial variable r, which is the spatial radius for
    //    three-dimensional actions, or the length of the four-dimensional
    //    Euclidean vector. This gets the correct bubble profile only if the
    //    transverse equations of motion in field space are also satisfied, or
    //    for a modified potential which has additional terms that raise the
    //    energy barrier on the side of the path that has a lower energy
    //    barrier. Hence this bubble profile gives an upper bound on the bounce
    //    action for the unmodified potential, as the modified potential (which
    //    has the same true vacuum and false vacuum and is exactly the same
    //    along the path) cannot have a lower bounce action.
    //    [The p equation of motion has a strange term giving a force
    //    proportional to (dp/dr)^2 which is not a friction term as the sign
    //    is not proportional to dp/dr, rather to d^2f/dp^2. This is because p
    //    is not linear in the distance in field space, so actually this term
    //    behaves a bit like extra kinetic energy because it rolled further
    //    "down the hill" from the true vacuum initially because there is more
    //    field space distance covered by a small change in p if the derivative
    //    is that way, or like extra friction because the approach to the
    //    false vacuum is actually longer in field space distance than in p.
    //    The alternative is to divide up the path into pathResolution segments
    //    and take the path as a series of straight lines at a constant
    //    "velocity" in field space, meaning that the weird pseudo-friction
    //    force in the bubble wall equation of motion in p disappears. This
    //    would require some contortions with linked lists (probably) of
    //    segments to get f(p(r)) and df/dp.]
    // 5) The bounce action along p(r) is numerically integrated and returned.
    virtual double
    operator()( std::vector< double > const& pathParameterization ) const;

    // This implements Up() for FCNBase just to stick to a basic value.
    virtual double Up() const { return 1.0; }


  protected:
    PotentialFunction const& potentialFunction;
    size_t const numberOfFields;
    size_t referenceFieldIndex;
    size_t const numberOfParameterizationFields;
    size_t const pathResolution;
    size_t potentialApproximationPower;
    PotentialMinimum const& falseVacuum;
    PotentialMinimum const& trueVacuum;
    double const fieldOriginPotential;
    double const falseVacuumPotential;
    double const trueVacuumPotential;
    double const falseVacuumEvaporationTemperature;
    double const minimumScaleSquared;
    double const tunnelingScaleSquared;
    double shortestLength;
    double longestLength;
    unsigned int const undershootOvershootAttempts;
    unsigned int const energyConservingUndershootAttempts;
    unsigned int const maximumMultipleOfLongestLength;
    double const initialFractionOfShortestLength;
    double const shootingThresholdSquared;

    // This turns a flattened matrix of numbers parameterizing the path from
    // the false vacuum to the true vacuum through field space.
    void
    DecodePathParameters( std::vector< double > const& pathParameterization,
                          std::vector< SimplePolynomial >& fieldsAsPolynomials,
                          double& thermalFalseVacuumPotential,
                          double& thermalTrueVacuumPotential ) const;

    // This fills fieldConfiguration with the values from fieldsAsPolynomials
    // at the auxiliary variable = auxiliaryValue. It also puts the value of
    // the potential (minus the value at the field origin at zero temperature)
    // into thermalFalseVacuumPotential and thermalTrueVacuumPotential for the
    // false and true vacua at the temperature given by
    // pathParameterization.back() if and only if it is non-zero. Unfortunately
    // it does not use "return value optimization" as we cannot be sure that
    // the user will compile with C++11 features enabled.
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

    // This returns true if neither overshooting nor undershooting definitely
    // happened and the distance in field space to the top of the false vacuum
    // hill is still larger than shootingCloseEnoughThreshold times the
    // distance between the initial field configuration and the false vacuum
    // field configuration.
    bool
    WorthIntegratingFurther( OdeintBubbleObserver const& odeintBubbleObserver,
             std::vector< SimplePolynomial >const& fieldsAsPolynomials ) const;
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
