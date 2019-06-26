/*
 * SplinePotential.hpp
 *
 *  Created on: Jun 10, 2014
 *      Authors: Ben O'Leary (benjamin.oleary@gmail.com)
 *      José Eliel Camargo-molina (elielcamargomolina@gmail.com)
 */

#ifndef SPLINEPOTENTIAL_HPP_
#define SPLINEPOTENTIAL_HPP_

#include "OneDimensionalPotentialAlongPath.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PathParameterization/TunnelPath.hpp"
#include <string>
#include <cstddef>
#include <vector>
#include "LHPC/Utilities/ParsingUtilities.hpp"

namespace VevaciousPlusPlus
{

  class SplinePotential : public OneDimensionalPotentialAlongPath
  {
  public:
    SplinePotential( PotentialFunction const& potentialFunction,
                     TunnelPath const& tunnelPath,
                     unsigned int const numberOfPotentialSegments,
                     double const minimumSquareDistanceBetweenPathVacua );
    virtual ~SplinePotential();


    virtual bool EnergyBarrierWasResolved() const
    { return energyBarrierWasResolved; }

    // This returns the value of the potential at auxiliaryValue, by finding
    // the correct segment and then returning its value at that point.
    virtual double operator()( double auxiliaryValue ) const;

    // This returns the value of the first derivative of the potential at
    // auxiliaryValue, by finding the correct segment and then returning its
    // slope at that point.
    virtual double FirstDerivative( double auxiliaryValue ) const;

    // This returns the value of the second derivative of the potential at
    // the false vacuum end of the path.
    virtual double SecondDerivativeAtFalseVacuum() const
    { return ( firstSegmentQuadratic + firstSegmentQuadratic ); }

    // This returns the value of the first derivative of the potential at
    // (auxiliaryOfPathPanicVacuum + differenceFromMaximumAuxiliary), assuming
    // that it is in the implicit final segment.
    virtual double FirstDerivativeNearPathPanic(
                            double const differenceFromMaximumAuxiliary ) const
    { return ( 2.0 * differenceFromMaximumAuxiliary * lastSegmentQuadratic ); }

    // This returns the value of the second derivative of the potential at
    // (definiteOvershootAuxiliary + differenceFromMaximumAuxiliary), assuming
    // that it is in the implicit final segment. The second derivative is
    // constant in the segment, so differenceFromMaximumAuxiliary is ignored.
    virtual double SecondDerivativeNearPathPanic(
                            double const differenceFromMaximumAuxiliary) const
    { return lastSegmentQuadratic; }

    // This is for debugging.
    virtual std::string AsDebuggingString() const;


  protected:
    double auxiliaryStep;
    double inverseOfAuxiliaryStep;
    // If the path is to be split into N segments, then there are ( N - 2 )
    // normal segments, which are taken as linear functions, with the false
    // vacuum side potential of each segment stored in potentialValues[ i - 2 ]
    // and the slope in firstDerivatives[ i - 2 ] for the ith segment (starting
    // from 1). The segment at the start and the segment at the end are taken
    // as quadratics without linear terms, parameterized by
    // firstSegmentQuadratic for the segment at the false vacuum and
    // finalPotential and lastSegmentQuadratic for the segment at the true
    // vacuum, with that segment being treated as being based at the end of the
    // path (truncated at the first path panic minimum).
    size_t numberOfNormalSegments;
    std::vector< double > potentialValues;
    std::vector< double > firstDerivatives;
    size_t pathFalseVacuumIndex;
    double pathFalsePotential;
    double firstSegmentQuadratic;
    double finalPotential;
    double lastSegmentQuadratic;
    bool energyBarrierWasResolved;
    bool tunnelingPossibleOnPath;
    std::vector< double > fieldConfiguration;
    PotentialFunction const& potentialFunction;
    TunnelPath const& tunnelPath;
    double const pathTemperature;
    // This is the threshold that is used to determine if a barrier between
    // false and true vacuum is actually a barrier and not some numerical
    // artifact. THe product of relativeBarrierThreshold and the biggest
    // absolute value between two points in the potential is used as
    // the threshold.
    double const relativeBarrierThreshold;

    // This returns the potential on tunnelPath at auxiliaryValue, relative to
    // the potential at the path false vacuum being zero. It puts values into
    // fieldConfiguration, so it is not const and should be used with caution.
    double CalculatePotentialDifference( double const auxiliaryValue );

    // This goes along the segment ends checking for the lowest before the
    // potential starts to rise again, returning true if it finds such a point
    // before the end of the path, false otherwise. It leaves the auxiliary
    // value of this point with lowest potential in auxiliaryOfPathFalseVacuum,
    // the lowest potential in pathFalsePotential, its field configuration in
    // fieldConfiguration, and the difference in potential at
    // ( auxiliaryOfPathFalseVacuum + auxiliaryStep ) along tunnelPath from
    // pathFalsePotential in potentialValues.front().
    bool RollForwardToPathFalseVacuum( size_t const maximumFalseIndex );

    // This goes along the segment ends adding the differences from
    // pathFalsePotential to potentialValues appropriately so that
    // potentialValues[ 0 ] is the potential at the start of the first straight
    // segment minus pathFalsePotential, potentialValues[ 1 ] at the start of
    // the second segment minus pathFalsePotential, and so on. This continues
    // until it finds the segment which starts with potential higher than that
    // of the path false vacuum and ends with potential lower than that of the
    // path false vacuum. It sets definiteUndershootAuxiliary to be the
    // auxiliary value of the start of this segment. It returns true if it
    // managed to do all the above, or false if tunneling was not possible (in
    // which case, definiteUndershootAuxiliary is a junk value).
    bool
    CheckForEndOfPositiveBarrier( size_t const maximumNumberOfNormalSegments );

    // This goes along the segment ends adding the differences from
    // pathFalsePotential to potentialValues appropriately so that
    // potentialValues[ 0 ] is the potential at the start of the first straight
    // segment minus pathFalsePotential, potentialValues[ 1 ] at the start of
    // the second segment minus pathFalsePotential, and so on. This continues
    // until it finds the path panic vacuum, which is either the start of the
    // first straight segment which has the potential at its start end lower
    // than pathFalsePotential and also has positive slope, or is the point at
    // the true vacuum end of tunnelPath. It returns false if it doesn't find
    // any point with potential lower than pathFalsePotential.
    double
    RollForwardToPathPanicVacuum( size_t const maximumNumberOfNormalSegments );
  };




  // This returns the value of the first derivative of the potential at
  // auxiliaryValue, by finding the correct segment and then returning its
  // slope at that point.
  inline double
  SplinePotential::FirstDerivative( double auxiliaryValue ) const
  {
    if( ( auxiliaryValue <= auxiliaryOfPathFalseVacuum )
        ||
        ( auxiliaryValue >= auxiliaryOfPathPanicVacuum ) )
    {
      return 0.0;
    }
    auxiliaryValue -= auxiliaryOfPathFalseVacuum;
    if( auxiliaryValue <= auxiliaryStep )
    {
      return ( 2.0 * auxiliaryValue * firstSegmentQuadratic );
    }
    size_t const auxiliarySteps( static_cast< size_t >( auxiliaryValue
                                                  * inverseOfAuxiliaryStep ) );
    if( auxiliarySteps > numberOfNormalSegments )
    {
      return FirstDerivativeNearPathPanic( auxiliaryValue
                                           + auxiliaryOfPathFalseVacuum
                                           - auxiliaryOfPathPanicVacuum );
    }
    return firstDerivatives[ auxiliarySteps - 1 ];
  }

  // This returns the potential on tunnelPath at auxiliaryValue, relative to
  // the potential at the path false vacuum being zero. It puts values into
  // fieldConfiguration, so it is not const and should be used with caution.
  inline double
  SplinePotential::CalculatePotentialDifference( double const auxiliaryValue )
  {
    tunnelPath.PutOnPathAt( fieldConfiguration,
                            auxiliaryValue );
    return ( potentialFunction( fieldConfiguration,
                                pathTemperature )
             - pathFalsePotential );
  }

  // This goes along the segment ends checking for the lowest before the
  // potential starts to rise again, returning true if it finds such a point
  // before the end of the path, false otherwise. It leaves the auxiliary
  // value of this point with lowest potential in auxiliaryOfPathFalseVacuum,
  // the lowest potential in pathFalsePotential, its field configuration in
  // fieldConfiguration, and the difference in potential at
  // ( auxiliaryOfPathFalseVacuum + auxiliaryStep ) along tunnelPath from
  // pathFalsePotential in potentialValues.front().
  inline bool SplinePotential::RollForwardToPathFalseVacuum(
                                               size_t const maximumFalseIndex )
  {
    while( pathFalseVacuumIndex <= maximumFalseIndex )
    {
      tunnelPath.PutOnPathAt( fieldConfiguration,
                              ( auxiliaryOfPathFalseVacuum + auxiliaryStep ) );
      potentialValues.front() = potentialFunction( fieldConfiguration,
                                                   pathTemperature );
      // Here we check whether the current potential value is higher than the
      // potential value at the false vacuum. If so, it means there is a barrier
      // as it has to go down eventually to reach the deeper panic vacuum.

      double maxAbsofPotentialValues =
              std::max(std::abs(potentialValues.front()), std::abs(pathFalsePotential));

      if( potentialValues.front() > pathFalsePotential &&
          potentialValues.front() - pathFalsePotential
                                       > relativeBarrierThreshold * maxAbsofPotentialValues )
      {
        potentialValues.front() -= pathFalsePotential;
        return true;
      }
      pathFalsePotential = potentialValues.front();
      auxiliaryOfPathFalseVacuum += auxiliaryStep;
      ++pathFalseVacuumIndex;
    }
    return false;
  }

} /* namespace VevaciousPlusPlus */
#endif /* SPLINEPOTENTIAL_HPP_ */
