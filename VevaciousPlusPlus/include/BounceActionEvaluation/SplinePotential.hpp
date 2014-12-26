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
    double pathFalsePotential;
    double firstSegmentQuadratic;
    double finalPotential;
    double lastSegmentQuadratic;
    bool energyBarrierWasResolved;
    bool tunnelingPossibleOnPath;


    // This goes along the segment ends checking for the lowest before the
    // potential starts to rise again, returning true if it finds such a point
    // before the end of the path, false otherwise. It leaves the auxiliary
    // value of this point with lowest potential in auxiliaryOfPathFalseVacuum,
    // the lowest potential in pathFalsePotential, its field configuration in
    // fieldConfiguration, and the potential at
    // ( auxiliaryOfPathFalseVacuum + auxiliaryStep ) along tunnelPath in
    // potentialValues.front().
    bool RollToPathLocalMinimum( PotentialFunction const& potentialFunction,
                                 TunnelPath const& tunnelPath,
                                 double const pathTemperature,
                                 std::vector< double >& fieldConfiguration );

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
    bool FindPathPanicVacuum( PotentialFunction const& potentialFunction,
                              TunnelPath const& tunnelPath,
                              double const pathTemperature,
                              std::vector< double >& fieldConfiguration );
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


  // This goes along the segment ends checking for the lowest before the
  // potential starts to rise again, returning true if it finds such a point
  // before the end of the path, false otherwise. It leaves the auxiliary
  // value of this point with lowest potential in auxiliaryOfPathFalseVacuum,
  // the lowest potential in pathFalsePotential, its field configuration in
  // fieldConfiguration, and the potential at
  // ( auxiliaryOfPathFalseVacuum + auxiliaryStep ) along tunnelPath in
  // potentialValues.front().
  inline bool SplinePotential::RollToPathLocalMinimum(
                                    PotentialFunction const& potentialFunction,
                                                  TunnelPath const& tunnelPath,
                                                  double const pathTemperature,
                                    std::vector< double >& fieldConfiguration )
  {
    // We need at least one step between auxiliaryOfPathFalseVacuum and the
    // true vacuum which means that we cannot go any closer than 2 steps away
    // from the true vacuum end of the path.
    double const maximumFalseAuxiliary( 1.0 - auxiliaryStep - auxiliaryStep );
    while( auxiliaryOfPathFalseVacuum < maximumFalseAuxiliary )
    {
      tunnelPath.PutOnPathAt( fieldConfiguration,
                              ( auxiliaryOfPathFalseVacuum + auxiliaryStep ) );
      potentialValues.front() = potentialFunction( fieldConfiguration,
                                                   pathTemperature );
      if( potentialValues.front() > pathFalsePotential )
      {
        return true;
      }
      pathFalsePotential = potentialValues.front();
      auxiliaryOfPathFalseVacuum += auxiliaryStep;
    }
    return false;
  }

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
  inline bool SplinePotential::FindPathPanicVacuum(
                                    PotentialFunction const& potentialFunction,
                                                  TunnelPath const& tunnelPath,
                                                  double const pathTemperature,
                                    std::vector< double >& fieldConfiguration )
  {
    NEED TO PUT CODE IN HERE!
  }

} /* namespace VevaciousPlusPlus */
#endif /* SPLINEPOTENTIAL_HPP_ */
