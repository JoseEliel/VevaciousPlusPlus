/*
 * SplinePotential.cpp
 *
 *  Created on: Jun 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/SplinePotential.hpp"

namespace VevaciousPlusPlus
{

  SplinePotential::SplinePotential( PotentialFunction const& potentialFunction,
                                    TunnelPath const& tunnelPath,
                                    size_t const numberOfPotentialSegments,
                         double const minimumSquareDistanceBetweenPathVacua ) :
    OneDimensionalPotentialAlongPath(),
    auxiliaryStep( 1.0 / static_cast< double >( numberOfPotentialSegments ) ),
    inverseOfAuxiliaryStep( numberOfPotentialSegments ),
    numberOfNormalSegments( numberOfPotentialSegments - 2 ),
    potentialValues( numberOfPotentialSegments,
                     0.0 ),
    firstDerivatives( numberOfNormalSegments,
                      0.0 ),
    pathFalseVacuumIndex( 0 ),
    pathFalsePotential( -1.0 ),
    firstSegmentQuadratic( -1.0 ),
    finalPotential( 0.0 ),
    lastSegmentQuadratic( -1.0 ),
    energyBarrierWasResolved( false ),
    tunnelingPossibleOnPath( false ),
    fieldConfiguration( potentialFunction.NumberOfFieldVariables() ),
    potentialFunction( potentialFunction ),
    tunnelPath( tunnelPath ),
    pathTemperature( tunnelPath.TemperatureValue() )
  {
    // First we have to find the path false minimum. The base constructor
    // already set auxiliaryOfPathPanicVacuum to zero.
    tunnelPath.PutOnPathAt( fieldConfiguration,
                            auxiliaryOfPathFalseVacuum );
    std::vector< double > const
    pathFalseEndConfiguration( fieldConfiguration );
    pathFalsePotential = potentialFunction( pathFalseEndConfiguration,
                                            pathTemperature );

    // We start at the assumed path false vacuum, then move along the resolved
    // segment ends to find any that are lower which are not separated by an
    // energy barrier which ends sufficiently far away from the given path
    // false vacuum end. The point where the energy barrier "ends" is where the
    // potential drops below what it is for the path false vacuum, and this
    // point must be separated in field space from the false vacuum end of the
    // path by at least sqrt(minimumSquareDistanceBetweenPathVacua).

    // The index runs from zero, so the segment at index
    // ( numberOfPotentialSegments - 1 ) ends in the true vacuum, and is
    // assumed to be a quadratic, so cannot start with the false vacuum, and
    // also by definition cannot end with an early panic vacuum.
    size_t const
    maximumIndexBeforeGivenTrueVauum( numberOfPotentialSegments - 2 );
    while( !energyBarrierWasResolved )
    {
      size_t const
      maximumNumberOfNormalSegments( maximumIndexBeforeGivenTrueVauum
                                     - pathFalseVacuumIndex );

      // The return value of RollForwardToLocalMinimum is false if it moves all
      // the way to the true vacuum end of the path without the potential ever
      // increasing. Here fieldConfiguration is being used as the field
      // configuration of the path false vacuum.
      tunnelingPossibleOnPath
      = RollForwardToPathFalseVacuum( maximumIndexBeforeGivenTrueVauum );
      if( !tunnelingPossibleOnPath )
      {
        // There is no point in continuing if tunneling is not possible.
        break;
      }

      // The return value of FindEndOfEnergyBarrier is false if it moves all
      // the way to the true vacuum end of the path without the potential ever
      // dropping below pathFalsePotential.
      tunnelingPossibleOnPath
      = CheckForEndOfPositiveBarrier( maximumNumberOfNormalSegments );

      if( !tunnelingPossibleOnPath )
      {
        // There is no point in continuing if tunneling is not possible.
        break;
      }

      auxiliaryOfPathPanicVacuum
      = RollForwardToPathPanicVacuum( maximumNumberOfNormalSegments );
      finalPotential = potentialValues[ numberOfNormalSegments + 1 ];

      // At this point, the path false vacuum is at auxiliaryOfPathFalseVacuum,
      // there are numberOfNormalSegments normal segments, the values in
      // potentialValues[ 0 ] to potentialValues[ numberOfNormalSegments - 1 ]
      // are the potentials at the starts of these segments (the remaining
      // elements are junk, but resizing would only save memory, which is not
      // critical, while it would take time, which is critical), and both
      // pathFalsePotential and finalPotential have been set.

      // Now we have to check that the tunneling is not over a spurious barrier
      // which might be due to numerical effects. The best we can do at the
      // moment is to check that the path panic vacuum is sufficiently far away
      // in field space from the false vacuum given by the false end of the
      // path.
      double vacuaSeparationSquared( 0.0 );
      double fieldDifference( 0.0 );
      for( size_t fieldIndex( 0 );
           fieldIndex < potentialFunction.NumberOfFieldVariables();
           ++fieldIndex )
      {
        fieldDifference = ( pathFalseEndConfiguration[ fieldIndex ]
                            - fieldConfiguration[ fieldIndex ] );
        vacuaSeparationSquared += ( fieldDifference * fieldDifference );
      }
      energyBarrierWasResolved
      = ( vacuaSeparationSquared >= minimumSquareDistanceBetweenPathVacua );

      if( !energyBarrierWasResolved )
      {
        // If the path panic vacuum wasn't sufficiently separated from the end
        // of tunnelPath, we set it to be the path false vacuum and iterate the
        // loop again. The potential in finalPotential is actually relative to
        // pathFalsePotential, hence adding it (it has a negative value) to
        // pathFalsePotential rather than setting pathFalsePotential to be
        // finalPotential.
        auxiliaryOfPathFalseVacuum = auxiliaryOfPathPanicVacuum;
        definiteUndershootAuxiliary = auxiliaryOfPathPanicVacuum;
        pathFalsePotential += finalPotential;
        continue;
      }

      // If we get here, then the energy barrier was resolved, the path false
      // vacuum is at auxiliaryOfPathFalseVacuum along tunnelPath, the
      // potential has dropped below pathFalsePotential between
      // definiteUndershootAuxiliary and
      // ( definiteUndershootAuxiliary + auxiliaryStep ), and the path panic
      // vacuum is at auxiliaryOfPathPanicVacuum along tunnelPath and is
      // sufficiently separated from the false vacuum end of tunnelPath, with
      // the potential at the path panic vacuum in finalPotential and in
      // potentialValues[ numberOfNormalSegments + 1 ]. Now we have to ensure
      // that the parabolic segments have the proper values and to set the
      // slopes of the straight segments.
      firstSegmentQuadratic = ( potentialValues.front()
                                * inverseOfAuxiliaryStep
                                * inverseOfAuxiliaryStep );
      thresholdForNearPathPanic
      = ( auxiliaryOfPathPanicVacuum - auxiliaryStep );
      for( size_t segmentIndex( 0 );
           segmentIndex < numberOfNormalSegments;
           ++segmentIndex )
      {
        firstDerivatives[ segmentIndex ]
        = ( ( potentialValues[ segmentIndex + 1 ]
              - potentialValues[ segmentIndex ] )
            * inverseOfAuxiliaryStep );
      }

      lastSegmentQuadratic
      = ( ( potentialValues[ numberOfNormalSegments ]
            - potentialValues[ numberOfNormalSegments + 1 ] )
          * inverseOfAuxiliaryStep * inverseOfAuxiliaryStep );
    }
  }

  SplinePotential::~SplinePotential()
  {
    // This does nothing.
  }


  // This returns the value of the potential at auxiliaryValue, by finding the
  // correct segment and then returning its value at that point.
  double SplinePotential::operator()( double auxiliaryValue ) const
  {
    if( auxiliaryValue >= auxiliaryOfPathPanicVacuum )
    {
      return finalPotential;
    }
    if( auxiliaryValue <= auxiliaryOfPathFalseVacuum )
    {
      return 0.0;
    }
    auxiliaryValue -= auxiliaryOfPathFalseVacuum;
    if( auxiliaryValue <= auxiliaryStep )
    {
      return ( auxiliaryValue * auxiliaryValue * firstSegmentQuadratic );
    }
    size_t const auxiliarySteps( auxiliaryValue * inverseOfAuxiliaryStep );
    if( auxiliarySteps > numberOfNormalSegments )
    {
      // We have to undo the "auxiliaryValue -= auxiliaryOfPathFalseVacuum" to
      // find the difference from the path panic point.
      double const differenceFromMaximumAuxiliary( auxiliaryOfPathPanicVacuum
                           - ( auxiliaryValue + auxiliaryOfPathFalseVacuum ) );
      return ( potentialValues[ numberOfNormalSegments + 1 ]
               + ( differenceFromMaximumAuxiliary
                   * differenceFromMaximumAuxiliary
                   * lastSegmentQuadratic ) );
    }
    double const auxiliaryDifference( auxiliaryValue
                                      - ( auxiliarySteps * auxiliaryStep ) );

    return ( potentialValues[ auxiliarySteps - 1 ]
          + ( auxiliaryDifference * firstDerivatives[ auxiliarySteps - 1 ] ) );
  }

  // This is for debugging.
  std::string SplinePotential::AsDebuggingString() const
  {
    std::stringstream returnStream;
    returnStream
    << "auxiliaryOfPathFalseVacuum = " << auxiliaryOfPathFalseVacuum
    << ", auxiliaryOfPathPanicVacuum = " << auxiliaryOfPathPanicVacuum
    << ", definiteUndershootAuxiliary = " << definiteUndershootAuxiliary
    << ", thresholdForNearPathPanic = " << thresholdForNearPathPanic
    << ", auxiliaryStep = " << auxiliaryStep
    << ", inverseOfAuxiliaryStep = " << inverseOfAuxiliaryStep
    << ", numberOfNormalSegments = " << numberOfNormalSegments
    << std::endl
    << std::endl;
    returnStream << "potential = ( UnitStep[x] * ( "
    << LHPC::ParsingUtilities::FormatNumberForMathematica(
                                                        firstSegmentQuadratic )
    << " * x^(2) ) * UnitStep["
    << LHPC::ParsingUtilities::FormatNumberForMathematica( auxiliaryStep )
    << " - x]";
    double cumulativeAuxiliary( auxiliaryStep );
    for( size_t segmentIndex( 0 );
         segmentIndex < numberOfNormalSegments;
         ++segmentIndex )
    {
      returnStream << " + "
      << "UnitStep[x - "
      << LHPC::ParsingUtilities::FormatNumberForMathematica(
                                                          cumulativeAuxiliary )
      << "] * ( "
      << LHPC::ParsingUtilities::FormatNumberForMathematica(
                                              potentialValues[ segmentIndex ] )
      << " + (x-"
      << LHPC::ParsingUtilities::FormatNumberForMathematica(
                                                          cumulativeAuxiliary )
      << ") * "
      << LHPC::ParsingUtilities::FormatNumberForMathematica(
                                             firstDerivatives[ segmentIndex ] )
      << " ) * UnitStep[";
      cumulativeAuxiliary += auxiliaryStep;
      returnStream
      << LHPC::ParsingUtilities::FormatNumberForMathematica(
                                                          cumulativeAuxiliary )
      << " - x]";
    }
    returnStream
    << " + UnitStep[x - "
    << LHPC::ParsingUtilities::FormatNumberForMathematica(
                                                          cumulativeAuxiliary )
    << "] * ( "
    << LHPC::ParsingUtilities::FormatNumberForMathematica( finalPotential )
    << " + (x-" << auxiliaryOfPathPanicVacuum << ")^2 * "
    << LHPC::ParsingUtilities::FormatNumberForMathematica(
                                                         lastSegmentQuadratic )
    << " ) * UnitStep["
    << LHPC::ParsingUtilities::FormatNumberForMathematica(
                                                   auxiliaryOfPathPanicVacuum )
    << " - x] )" << std::endl
    << std::endl;
    returnStream
    << "potentialValues = { ";
    for( size_t segmentIndex( 0 );
         segmentIndex <= numberOfNormalSegments;
         ++segmentIndex )
    {
      if( segmentIndex > 0 )
      {
        returnStream << ",";
      }
      returnStream << " "
      << LHPC::ParsingUtilities::FormatNumberForMathematica(
                                             potentialValues[ segmentIndex ] );
    }
    returnStream << std::endl << "(then {";
    for( size_t segmentIndex( numberOfNormalSegments + 1 );
         segmentIndex < potentialValues.size();
         ++segmentIndex )
    {
      if( segmentIndex > ( numberOfNormalSegments + 1 ) )
      {
        returnStream << ",";
      }
      returnStream << " "
      << LHPC::ParsingUtilities::FormatNumberForMathematica(
                                             potentialValues[ segmentIndex ] );
    }
    returnStream << " })";
    returnStream << std::endl
    << "firstDerivatives = { ";
    for( size_t segmentIndex( 0 );
         segmentIndex < numberOfNormalSegments;
         ++segmentIndex )
    {
      if( segmentIndex > 0 )
      {
        returnStream << ",";
      }
      returnStream << " "
      << LHPC::ParsingUtilities::FormatNumberForMathematica(
                                            firstDerivatives[ segmentIndex ] );
    }
    returnStream << std::endl << "(then {";
    for( size_t segmentIndex( numberOfNormalSegments );
         segmentIndex < firstDerivatives.size();
         ++segmentIndex )
    {
      if( segmentIndex > numberOfNormalSegments )
      {
        returnStream << ",";
      }
      returnStream << " "
      << LHPC::ParsingUtilities::FormatNumberForMathematica(
                                            firstDerivatives[ segmentIndex ] );
    }
    returnStream << " })";
    returnStream << std::endl;
    returnStream
    << "pathFalsePotential = " << pathFalsePotential
    << ", firstSegmentQuadratic = " << firstSegmentQuadratic
    << ", lastSegmentQuadratic = " << lastSegmentQuadratic
    << ", energyBarrierWasResolved = " << energyBarrierWasResolved
    << ", tunnelingPossibleOnPath = " << tunnelingPossibleOnPath
    << ", fieldConfiguration not really relevant"
    << ", potentialFunction not really printable" << std::endl;
    returnStream << std::endl;
    returnStream << "tunnelPath = " << tunnelPath.AsDebuggingString()
    << ", pathTemperature = " << pathTemperature
    << std::endl
    << std::endl;

    return returnStream.str();
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
  bool SplinePotential::CheckForEndOfPositiveBarrier(
                                   size_t const maximumNumberOfNormalSegments )
  {
    // When this is called, the path false vacuum is at
    // auxiliaryOfPathFalseVacuum, and the potential at
    // ( auxiliaryOfPathFalseVacuum + auxiliaryStep ) is already in
    // potentialValues.front().
    definiteUndershootAuxiliary
    = ( auxiliaryOfPathFalseVacuum + auxiliaryStep );
    numberOfNormalSegments = 0;

    while( numberOfNormalSegments < maximumNumberOfNormalSegments )
    {
      double const segmentEndPotential( CalculatePotentialDifference(
                               definiteUndershootAuxiliary + auxiliaryStep ) );
      potentialValues[ ++numberOfNormalSegments ] = segmentEndPotential;

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
    // false vacuum. The potential difference at the end of the path relative
    // to the path false vacuum is placed in
    // potentialValues[ numberOfNormalSegments + 1 ] (which is not a
    // problem, as potentialValues was initialized with at least
    // ( numberOfNormalSegments + 2 ) elements, so the index is not out
    // of range).
    potentialValues[ numberOfNormalSegments + 1 ]
    = ( CalculatePotentialDifference( auxiliaryOfPathPanicVacuum ) );

    return ( potentialValues[ numberOfNormalSegments + 1 ] < 0.0 );
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
  double SplinePotential::RollForwardToPathPanicVacuum(
                                   size_t const maximumNumberOfNormalSegments )
  {
    // The number of normal segments up to and including the first segment
    // which ended below the path false vacuum is currently held in
    // numberOfNormalSegments, from the last call of
    // CheckForEndOfPositiveBarrier().
    if( numberOfNormalSegments == 0 )
    {
      // This is the case where the path false vacuum was only found to be
      // separated from the end of the path by twice auxiliaryStep.
      return 1.0;
    }
    // Otherwise we can check at least the end of the first normal segment.
    // Hence we start by assuming that the path panic vacuum is at the end of
    // the first segment after the first which went below the path false
    // potential.
    double panicAuxiliary( definiteUndershootAuxiliary
                           + auxiliaryStep );

    while( numberOfNormalSegments < maximumNumberOfNormalSegments )
    {
      // We look at the end of latest "added" segment.
      panicAuxiliary += auxiliaryStep;
      double const
      segmentEndPotential( CalculatePotentialDifference( panicAuxiliary ) );
      if( segmentEndPotential < potentialValues[ numberOfNormalSegments ] )
      {
        // If the next segment starts lower, the path  panic vacuum has to be
        // at the beginning of that segment or further along, and we note the
        // potential at the start of that segment and mark it as a normal
        // segment.
        potentialValues[ ++numberOfNormalSegments ] = segmentEndPotential;
      }
      else
      {
        // If the next segment starts higher, then the last segment marked as
        // normal is that which ends in the path panic vacuum. The potential of
        // the path panic vacuum is kept in
        // potentialValues[ numberOfNormalSegments + 1 ], accounting for the
        // decrement of numberOfNormalSegments.
        --numberOfNormalSegments;
        return ( panicAuxiliary - auxiliaryStep );
      }
    }

    // If we get to here, the potential just kept decreasing until one
    // segment length before the end of the path. The potential difference at
    // the end of the path relative to the path false vacuum is placed in
    // potentialValues[ numberOfNormalSegments + 1 ] (which is not a
    // problem, as potentialValues was initialized with at least
    // ( numberOfNormalSegments + 2 ) elements, so the index is not out
    // of range).
    potentialValues[ numberOfNormalSegments + 1 ]
    = CalculatePotentialDifference( 1.0 );
    if( potentialValues[ numberOfNormalSegments + 1 ]
        < potentialValues[ numberOfNormalSegments ] )
    {
      // If the path end is the path panic vacuum, then
      // potentialValues[ numberOfNormalSegments + 1 ] has been correctly set.
      return 1.0;
    }
    else
    {
      // Otherwise decrementing numberOfNormalSegments leaves the potential
      // at the path panic vacuum still in
      // potentialValues[ numberOfNormalSegments + 1 ], and the value in
      // potentialValues[ numberOfNormalSegments + 2 ] (from
      // CalculatePotentialDifference( 1.0 )) is just ignored.
      --numberOfNormalSegments;
      return panicAuxiliary;
    }
  }

} /* namespace VevaciousPlusPlus */
