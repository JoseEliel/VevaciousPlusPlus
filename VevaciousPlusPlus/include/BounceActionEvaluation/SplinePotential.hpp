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
    { return energyBarrierWasResolved; }

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
    size_t numberOfPotentialSegments;
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
    size_t pathFalseVacuumIndex;
    size_t pathPanicVacuumIndex;
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
    bool RollForwardToLocalMinimum();

    // This goes along the segment ends adding the differences from
    // pathFalsePotential to potentialValues appropriately so that
    // potentialValues[ 0 ] is the potential at the start of the first straight
    // segment minus pathFalsePotential, potentialValues[ 1 ] at the start of
    // the second segment minus pathFalsePotential, and so on. This continues
    // until it finds the segment which starts with potential higher than that of
    // the path false vacuum and ends with potential lower than that of the path
    // false vacuum. It sets definiteUndershootAuxiliary to be the auxiliary
    // value of the start of this segment. It returns true if it managed to do
    // all the above, or false if tunneling was not possible (in which case,
    // definiteUndershootAuxiliary is a junk value).
    bool CheckForEndOfPositiveBarrier();

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
    double FindPathPanicVacuum();
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


  // This returns the potential on tunnelPath at auxiliaryValue, relative to
  // the potential at the path false vacuum being zero. It puts values into
  // fieldConfiguration, so it is not const and should be used with caution.
  inline double
  SplinePotential::CalculatePotentialDifference( double const auxiliaryValue )
  {
    tunnelPath.PutOnPathAt( fieldConfiguration,
                            auxiliaryValue );
    return ( potentialFunction( auxiliaryValue,
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
  inline bool SplinePotential::RollForwardToLocalMinimum()
  {
    // We need at least one step between auxiliaryOfPathFalseVacuum and the
    // true vacuum which means that we cannot go any closer than 2 steps away
    // from the true vacuum end of the path.
    size_t const maximumFalseIndex( numberOfPotentialSegments - 2 );
    while( pathFalseVacuumIndex <= maximumFalseIndex )
    {
      tunnelPath.PutOnPathAt( fieldConfiguration,
                              ( auxiliaryOfPathFalseVacuum + auxiliaryStep ) );
      potentialValues.front() = potentialFunction( fieldConfiguration,
                                                   pathTemperature );
      if( potentialValues.front() > pathFalsePotential )
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

  // This goes along the segment ends adding the differences from
  // pathFalsePotential to potentialValues appropriately so that
  // potentialValues[ 0 ] is the potential at the start of the first straight
  // segment minus pathFalsePotential, potentialValues[ 1 ] at the start of
  // the second segment minus pathFalsePotential, and so on. This continues
  // until it finds the segment which starts with potential higher than that of
  // the path false vacuum and ends with potential lower than that of the path
  // false vacuum. It sets definiteUndershootAuxiliary to be the auxiliary
  // value of the start of this segment. It returns true if it managed to do
  // all the above, or false if tunneling was not possible (in which case,
  // definiteUndershootAuxiliary is a junk value).
  inline bool SplinePotential::CheckForEndOfPositiveBarrier()
  {
    // When this is called, the path false vacuum is at
    // auxiliaryOfPathFalseVacuum, and the potential at
    // ( auxiliaryOfPathFalseVacuum + auxiliaryStep ) is already in
    // potentialValues.front().
    size_t const maximumNumberOfNormalSegments( numberOfPotentialSegments
                                                - pathFalseVacuumIndex
                                                - 2 );
    definiteUndershootAuxiliary
    = ( auxiliaryOfPathFalseVacuum + auxiliaryStep );
    numberOfNormalSegments = 0;

    NEED TO CONSIDER CASE OF path false - step - barrier - step - path panic!

    while( numberOfNormalSegments < maximumNumberOfNormalSegments )
    {
      double const segmentEndPotential( CalculatePotentialDifference(
                               definiteUndershootAuxiliary + auxiliaryStep ) );
      potentialValues[ numberOfNormalSegments ] = segmentEndPotential;
      ++numberOfNormalSegments;

      if( segmentEndPotential < 0 )
      {
        // If this segment starts higher than the path false vacuum (which it
        // must, or else we would have already returned true) and ends lower,
        // then we have found the end of the barrier of greater energy than the
        // path false vacuum.
        return true;
      }

      definiteUndershootAuxiliary += auxiliaryStep;
    }

    // If the loop ended without the function already returning true, then
    // either the very end of tunnelPath is the path panic vacuum, or tunneling
    // is not possible because we did not find a point lower than the path
    // false vacuum.
    definiteUndershootAuxiliary -= auxiliaryStep;
    finalPotential
    = ( CalculatePotentialDifference( auxiliaryOfPathPanicVacuum ) );

    return ( finalPotential < 0.0 );
  }

  // This goes along the segment ends adding the differences from
  // pathFalsePotential to potentialValues appropriately so that
  // potentialValues[ n ] is the potential at the start of the nth (counting
  // from 0) straight segment minus pathFalsePotential. This continues
  // until it finds the path panic vacuum, which is either the start of the
  // first straight segment which has the potential at its start end lower
  // than pathFalsePotential and also has positive slope, or is the point at
  // the true vacuum end of tunnelPath, and returns the auxiliary value for
  // this.
  inline double SplinePotential::FindPathPanicVacuum()
  {
    GO THROUGH LOGIC AGAIN, ENSURE THAT VALUES ARE PUT INTO CORRECT PLACES IN
    potentialValues WITHOUT USING push_back().
    // When this is called, potentialValues.back() is the potential at
    // ( definiteUndershootAuxiliary + auxiliaryStep ).
    double segmentStartAuxiliary
    = ( definiteUndershootAuxiliary + auxiliaryStep + auxiliaryStep );

    // If we start too close to the end of tunnelPath, we just
    if( segmentStartAuxiliary > thresholdForNearPathPanic )
    {
      finalPotential = ( CalculatePotentialDifference( 1.0 ) );
      return 1.0;
    }

    // We look at the entire segment before adding anything to potentialValues.
    double segmentStartPotential( potentialValues.back() );
    double segmentEndPotential( CalculatePotentialDifference(
                                     segmentStartAuxiliary + auxiliaryStep ) );

    if( segmentEndPotential > segmentStartPotential )
    {
      return segmentStartAuxiliary;
    }

    while( segmentStartAuxiliary <= thresholdForNearPathPanic )
    {
      if( segmentEndPotential > segmentStartPotential )
      {
        // If this segment ends higher than it starts, it is a minimum on the
        // path which can be tunneled to from the path false vacuum.
        finalPotential = segmentStartPotential;
        return segmentStartAuxiliary;
      }

      // If we have not already returned from this function, then either we
      // are not yet through the barrier, or the end of the segment is deeper,
      // meaning that we have not yet found the path panic vacuum. Either way,
      // we record the potential at the end of this segment (as the start of
      // the next segment).
      potentialValues[ numberOfNormalSegments ] = segmentEndPotential;
      ++numberOfNormalSegments;

      // If the start of this segment was above the path false vacuum and its
      // end drops below the path false vacuum, we note that we are through the
      // barrier and that the start of the segment was the largest auxiliary
      // which we are sure is going to undershoot.
      if( !throughBarrier
          &&
          ( nextPotential < 0 ) )
      {
        throughBarrier = true;
        definiteUndershootAuxiliary = auxiliaryOfPathPanicVacuum;
      }

      auxiliaryOfPathPanicVacuum += auxiliaryStep;
    }

    // If the loop ended without the function already returning true, then
    // either the very end of tunnelPath is the path panic vacuum, or tunneling
    // is not possible because we did not find a point lower than the path
    // false vacuum.
    tunnelPath.PutOnPathAt( fieldConfiguration,
                            1.0 );
    finalPotential = ( potentialFunction( fieldConfiguration,
                                          pathTemperature )
                       - pathFalsePotential );

    // We return false if tunneling is not possible.
    if( finalPotential >= 0.0 )
    {
      return false;
    }

    // If the potential already dropped below the path false vacuum in a
    // straight segment, throughBarrier is already true.


    NEED TO PUT CODE IN HERE!
  }

} /* namespace VevaciousPlusPlus */
#endif /* SPLINEPOTENTIAL_HPP_ */
