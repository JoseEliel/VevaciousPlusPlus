/*
 * SplinePotential.hpp
 *
 *  Created on: Jun 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SPLINEPOTENTIAL_HPP_
#define SPLINEPOTENTIAL_HPP_

#include "CommonIncludes.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PathParameterization/TunnelPath.hpp"
#include "OneDimensionalPotentialAlongPath.hpp"

namespace VevaciousPlusPlus
{

  class SplinePotential : public OneDimensionalPotentialAlongPath
  {
  public:
    SplinePotential( PotentialFunction const& potentialFunction,
                     TunnelPath const& tunnelPath,
                     size_t const numberOfPotentialSegments,
                     double const minimumSquareDistanceBetweenPathVacua );
    virtual ~SplinePotential();


    virtual bool EnergyBarrierWasResolved() const
    {return energyBarrierWasResolved; }

    // This returns the value of the potential at auxiliaryValue, by finding
    // the correct segment and then returning its value at that point.
    virtual double operator()( double auxiliaryValue ) const;

    // This returns the value of the first derivative of the potential at
    // auxiliaryValue, by finding the correct segment and then returning its
    // slope at that point.
    virtual double FirstDerivative( double const auxiliaryValue ) const;

    // This returns the value of the second derivative of the potential at
    // the false vacuum end of the path.
    virtual double SecondDerivativeAtFalseVacuum() const
    { return ( firstSegmentQuadratic + firstSegmentQuadratic ); }

    // This returns the value of the first derivative of the potential at
    // (definiteOvershootAuxiliary + differenceFromMaximumAuxiliary), assuming
    // that it is in the implicit final segment.
    virtual double FirstDerivativeNearPathPanic(
                            double const differenceFromMaximumAuxiliary ) const
    { return ( 2.0 * differenceFromMaximumAuxiliary * lastSegmentQuadratic ); }

    // This returns the value of the second derivative of the potential at
    // (definiteOvershootAuxiliary + differenceFromMaximumAuxiliary), assuming
    // that it is in the implicit final segment.
    virtual double SecondDerivativeNearPathPanic() const
    { return lastSegmentQuadratic; }

    // This is for debugging.
    virtual std::string AsDebuggingString() const;


  protected:
    double auxiliaryStep;
    double inverseOfAuxiliaryStep;
    // There are ( numberOfPotentialSegments - 2 ) normal segments, which are
    // taken as linear functions, with the false vacuum side potential of each
    // segment stored in potentialValues[ i - 2 ] and the slope in
    // firstDerivatives[ i - 2 ] for the ith segment (starting from 1). The
    // segment at the start and the segment at the end are taken as quadratics
    // without linear terms, parameterized by firstSegmentQuadratic for the
    // segment at the false vacuum and finalPotential and lastSegmentQuadratic
    // for the segment at the true vacuum, with that segment being treated as
    // being based at the end of the path (truncated at the first path panic
    // minimum).
    size_t numberOfNormalSegments;
    std::vector< double > potentialValues;
    std::vector< double > firstDerivatives;
    double firstSegmentQuadratic;
    double finalPotential;
    double lastSegmentQuadratic;
    bool energyBarrierWasResolved;
    bool tunnelingPossibleOnPath;
  };




  // This returns the value of the first derivative of the potential at
  // auxiliaryValue, by finding the correct segment and then returning its
  // slope at that point.
  inline double
  SplinePotential::FirstDerivative( double const auxiliaryValue ) const
  {
    if( ( auxiliaryValue <= 0.0 )
        ||
        ( auxiliaryValue >= definiteOvershootAuxiliary ) )
    {
      return 0.0;
    }
    if( auxiliaryValue < auxiliaryStep )
    {
      return ( 2.0 * auxiliaryValue * firstSegmentQuadratic );
    }
    if( auxiliaryValue >= thresholdForNearPathPanic )
    {
      return FirstDerivativeNearPathPanic( auxiliaryValue
                                           - definiteOvershootAuxiliary );
    }
    return
    firstDerivatives[ size_t( auxiliaryValue * inverseOfAuxiliaryStep ) - 1 ];
  }

} /* namespace VevaciousPlusPlus */
#endif /* SPLINEPOTENTIAL_HPP_ */
