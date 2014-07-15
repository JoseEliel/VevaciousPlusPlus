/*
 * SplinePotential.hpp
 *
 *  Created on: Jun 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SPLINEPOTENTIAL_HPP_
#define SPLINEPOTENTIAL_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class SplinePotential
  {
  public:
    SplinePotential( double const minimumFalseVacuumConcavity = 1.0e-6 );
    SplinePotential( SplinePotential const& copySource );
    virtual
    ~SplinePotential();


    // This returns the value of the potential at auxiliaryValue, by finding
    // the correct segment and then returning its value at that point.
    virtual double operator()( double auxiliaryValue ) const;

    // This returns the value of the first derivative of the potential at
    // auxiliaryValue, by finding the correct segment and then returning its
    // slope at that point.
    virtual double FirstDerivative( double const auxiliaryValue ) const;

    // This returns the value of the second derivative of the potential at
    // auxiliaryValue, by finding the correct segment and then returning its
    // slope at that point.
    virtual double SecondDerivative( double const auxiliaryValue ) const;

    // This returns the value of the first derivative of the potential at
    // (definiteOvershootAuxiliary + differenceFromMaximumAuxiliary), assuming
    // that it is in the implicit final segment.
    virtual double FirstDerivativeNearPathPanic(
                           double const differenceFromMaximumAuxiliary ) const;

    // This returns the value of the second derivative of the potential at
    // (definiteOvershootAuxiliary + differenceFromMaximumAuxiliary), assuming
    // that it is in the implicit final segment.
    virtual double SecondDerivativeNearPathPanic(
                           double const differenceFromMaximumAuxiliary ) const;

    // This adds another point for the spline, assuming that it goes after the
    // previously-added point by auxiliaryDifference from the previous point,
    // to a potential value of potentialValue relative to the false vacuum.
    void AddPoint( double const auxiliaryValue,
                   double const potentialValue );

    // This sets up the spline based on auxiliaryValues and potentialValues,
    // ensuring that the potential reaches the correct values for the vacua,
    // and that the potential derivative vanishes at the vacua. It also notes
    // the first point where the potential drops below that of the false vacuum
    // in definiteUndershootAuxiliary and the first maximum after that in
    // definiteOvershootAuxiliary, and cuts off the potential at that maximum.
    void SetSpline( double const trueVacuumPotentialDifference );

    double DefiniteUndershootAuxiliary() const
    { return definiteUndershootAuxiliary; }

    double DefiniteOvershootAuxiliary() const
    { return definiteOvershootAuxiliary; }

    double StartOfFinalSegment() const{ return startOfFinalSegment; }

    double SizeOfFinalSegment() const{ return sizeOfFinalSegment; }

    // This is for debugging.
    std::string AsDebuggingString() const;


  protected:
    // There are auxiliaryValues.size() normal segments, each represented by 4
    // values: the size of the segment in the auxiliary value, the potential
    // (relative to the false vacuum at zero auxiliary value), its first
    // derivative, and its second derivative. Within the segment, the potential
    // is approximated by
    // V(p_j + d) + d * V'(p_j) + d^2 * [0.5*V''(p_j)],
    // where the auxiliary value is p = p_j + d. There is also one final
    // segment, where for an auxiliary value = 1 - d, the potential is
    // approximated by
    // V(1.0) + d^2 * [0.5*V''(1)] + d^4 * [V''''(1)/(4*3*2)].
    std::vector< double > auxiliaryValues;
    std::vector< double > potentialValues;
    std::vector< double > firstDerivatives;
    std::vector< double > halfSecondDerivatives;
    double finalPotential;
    double halfFinalSecondDerivative;
    double finalCubicCoefficient;
    double const minimumFalseVacuumConcavity;
    double definiteUndershootAuxiliary;
    double definiteOvershootAuxiliary;
    double auxiliaryUpToCurrentSegment;
    double startOfFinalSegment;
    double sizeOfFinalSegment;
  };




  // This returns the value of the first derivative of the potential at
  // (definiteOvershootAuxiliary + differenceFromMaximumAuxiliary), assuming
  // that it is in the implicit final segment.
  inline double SplinePotential::FirstDerivativeNearPathPanic(
                            double const differenceFromMaximumAuxiliary ) const
  {
    return ( ( ( 2.0 * halfFinalSecondDerivative )
                   + ( 3.0 * finalCubicCoefficient
                           * differenceFromMaximumAuxiliary ) )
             * differenceFromMaximumAuxiliary );
  }

  // This returns the value of the second derivative of the potential at
  // (definiteOvershootAuxiliary + differenceFromMaximumAuxiliary), assuming
  // that it is in the implicit final segment.
  inline double SplinePotential::SecondDerivativeNearPathPanic(
                            double const differenceFromMaximumAuxiliary ) const
  {
    return ( ( 2.0 * halfFinalSecondDerivative )
               + ( 6.0 * finalCubicCoefficient
                       * differenceFromMaximumAuxiliary ) );
  }

  // This adds another point for the spline, assuming that it goes after the
  // previously-added point by auxiliaryDifference from the previous point,
  // to a potential value of potentialValue relative to the false vacuum.
  inline void SplinePotential::AddPoint( double const auxiliaryValue,
                                         double const potentialValue )
  {
    auxiliaryValues.push_back( auxiliaryValue - auxiliaryUpToCurrentSegment );
    potentialValues.push_back( potentialValue );
    auxiliaryUpToCurrentSegment = auxiliaryValue;
  }

} /* namespace VevaciousPlusPlus */
#endif /* SPLINEPOTENTIAL_HPP_ */
