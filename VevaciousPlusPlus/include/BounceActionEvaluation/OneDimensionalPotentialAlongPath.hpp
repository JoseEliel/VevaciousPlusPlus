/*
 * OneDimensionalPotentialAlongPath.hpp
 *
 *  Created on: Nov 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef ONEDIMENSIONALPOTENTIALALONGPATH_HPP_
#define ONEDIMENSIONALPOTENTIALALONGPATH_HPP_

#include "CommonIncludes.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PathParameterization/TunnelPath.hpp"

namespace VevaciousPlusPlus
{

  class OneDimensionalPotentialAlongPath
  {
  public:
    OneDimensionalPotentialAlongPath();
    virtual ~OneDimensionalPotentialAlongPath();


    // This should return the value of the potential at auxiliaryValue.
    virtual double operator()( double auxiliaryValue ) const = 0;

    // This should returns the value of the first derivative with respect to
    // the path auxiliary of the potential at auxiliaryValue.
    virtual double FirstDerivative( double const auxiliaryValue ) const = 0;

    // This returns the value of the second derivative of the potential at
    // the false vacuum end of the path.
    virtual double SecondDerivativeAtFalseVacuum() const = 0;

    double DefiniteUndershootAuxiliary() const
    { return definiteUndershootAuxiliary; }

    double DefiniteOvershootAuxiliary() const
    { return definiteOvershootAuxiliary; }

    // This should return the value of the first derivative of the potential at
    // (definiteOvershootAuxiliary + differenceFromMaximumAuxiliary), assuming
    // that it is greater than thresholdForNearPathPanic. It is assumed that
    // differenceFromMaximumAuxiliary is negative.
    virtual double FirstDerivativeNearPathPanic(
                       double const differenceFromMaximumAuxiliary ) const = 0;

    // This should return the value of the second derivative of the potential
    // at (definiteOvershootAuxiliary + differenceFromMaximumAuxiliary),
    // assuming that it is greater than thresholdForNearPathPanic. It is
    // assumed that differenceFromMaximumAuxiliary is negative.
    virtual double SecondDerivativeNearPathPanic() const = 0;

    // This returns what should be the path auxiliary value which is close
    // enough to the path panic minimum that special care must be taken in
    // terms of how accurate a double can be.
    double ThresholdForNearPathPanic() const
    { return thresholdForNearPathPanic; }


    // This is for debugging.
    std::string AsDebuggingString() const;


  protected:
    double definiteUndershootAuxiliary;
    double definiteOvershootAuxiliary;
    double thresholdForNearPathPanic;
  };

} /* namespace VevaciousPlusPlus */
#endif /* ONEDIMENSIONALPOTENTIALALONGPATH_HPP_ */
