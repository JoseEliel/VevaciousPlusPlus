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
                                    double const falseVacuumPotential,
                                    double const trueVacuumPotential ) :
    energyBarrierResolved( false ),
    trueVacuumLowerThanPathFalseMinimum( false ),
    auxiliaryStep( 1.0 / static_cast< double >( numberOfPotentialSegments ) ),
    inverseOfAuxiliaryStep( numberOfPotentialSegments ),
    numberOfNormalSegments( numberOfPotentialSegments - 2 ),
    potentialValues( numberOfNormalSegments,
                     NAN ),
    firstDerivatives( numberOfNormalSegments,
                      NAN ),
    firstSegmentQuadratic( NAN ),
    finalPotential( trueVacuumPotential - falseVacuumPotential ),
    lastSegmentQuadratic( NAN ),
    definiteUndershootAuxiliary( NAN ),
    definiteOvershootAuxiliary( 1.0 ),
    startOfFinalSegment( NAN )
  {
    // First we have to find the path false minimum.
    double pathFalsePotential( falseVacuumPotential );
    double segmentEndAuxiliary( auxiliaryStep );
    std::vector< double >
    fieldConfiguration( potentialFunction.NumberOfFieldVariables() );
    tunnelPath.PutOnPathAt( fieldConfiguration,
                            segmentEndAuxiliary );
    double const pathTemperature( tunnelPath.TemperatureValue() );
    double segmentEndPotential( potentialFunction( fieldConfiguration,
                                                   pathTemperature ) );
    size_t segmentIndex( 0 );
    while( segmentIndex < ( numberOfPotentialSegments - 1 ) )
    {
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "segmentIndex = " << segmentIndex << ", pathFalsePotential = "
      << pathFalsePotential << ", segmentEndPotential = "
      << segmentEndPotential;
      std::cout << std::endl;/**/

      // If we find the start of the energy barrier, we note so and break from
      // the loop.
      if( segmentEndPotential > pathFalsePotential )
      {
        energyBarrierResolved = true;
        firstSegmentQuadratic = ( ( segmentEndPotential - pathFalsePotential )
                           * inverseOfAuxiliaryStep * inverseOfAuxiliaryStep );
        break;
      }
      // Otherwise we move onto the next segment to check for the start of the
      // barrier.
      pathFalsePotential = segmentEndPotential;
      ++segmentIndex;
      segmentEndAuxiliary += auxiliaryStep;
      tunnelPath.PutOnPathAt( fieldConfiguration,
                              segmentEndAuxiliary );
      segmentEndPotential = potentialFunction( fieldConfiguration,
                                               pathTemperature );
    }
    // Now we move on to resolving the energy barrier, looking out for an early
    // path panic minimum. However, this is only done if we resolved the start
    // of an energy barrier.
    if( energyBarrierResolved )
    {
      double const auxiliaryStartOfNormalSegments( segmentEndAuxiliary );
      ++segmentIndex;
      segmentEndAuxiliary += auxiliaryStep;
      // Now segmentEndAuxiliary and segmentIndex are correct for the first
      // straight segment.
      size_t normalSegmentIndex( 0 );
      potentialValues.front() = ( segmentEndPotential - pathFalsePotential );
      tunnelPath.PutOnPathAt( fieldConfiguration,
                              segmentEndAuxiliary );
      // From this point on, pathFalsePotential is always subtracted from
      // potentials along the path.
      segmentEndPotential
      = ( potentialFunction( fieldConfiguration,
                             pathTemperature ) - pathFalsePotential );
      firstDerivatives.front()
      = ( segmentEndPotential - potentialValues.front() );
      bool foundDefiniteUndershoot( false );

      while( segmentIndex < ( numberOfPotentialSegments - 2 ) )
      {
        // debugging:
        /**/std::cout << std::endl << "debugging:"
        << std::endl
        << "segmentIndex = " << segmentIndex << ", normalSegmentIndex = "
        << normalSegmentIndex << ", potentialValues[ " << normalSegmentIndex
        << " ] = " << potentialValues[ normalSegmentIndex ]
        << ", firstDerivatives[ " << normalSegmentIndex << " ] = "
        << firstDerivatives[ normalSegmentIndex ];
        std::cout << std::endl;/**/

        // If we find an early path panic minimum, we note it and break early.
        if( ( potentialValues[ normalSegmentIndex ] < 0.0 )
            &&
            ( firstDerivatives[ normalSegmentIndex ] > 0.0 ) )
        {
          trueVacuumLowerThanPathFalseMinimum = true;
          finalPotential = potentialValues[ normalSegmentIndex ];
          lastSegmentQuadratic
          = ( ( potentialValues[ normalSegmentIndex - 1 ] - finalPotential )
                           * inverseOfAuxiliaryStep * inverseOfAuxiliaryStep );
          numberOfNormalSegments = normalSegmentIndex;
          definiteOvershootAuxiliary
          = ( ( numberOfNormalSegments + 2 ) * auxiliaryStep );
          break;
        }
        // Otherwise we move onto the next segment.
        ++segmentIndex;
        ++normalSegmentIndex;
        segmentEndAuxiliary += auxiliaryStep;
        potentialValues[ normalSegmentIndex ] = segmentEndPotential;
        tunnelPath.PutOnPathAt( fieldConfiguration,
                                segmentEndAuxiliary );
        // From this point on, pathFalsePotential is always subtracted from
        // potentials along the path.
        segmentEndPotential
        = ( potentialFunction( fieldConfiguration,
                               pathTemperature ) - pathFalsePotential );
        firstDerivatives[ normalSegmentIndex ]
        = ( segmentEndPotential - potentialValues.front() );
      }
      // Now we have to check to see if the end of the path is lower than the
      // path false minimum, if we haven't already found an early path panic
      // minimum.
      if( !trueVacuumLowerThanPathFalseMinimum )
      {
        if( trueVacuumPotential < pathFalsePotential )
        {
          trueVacuumLowerThanPathFalseMinimum = true;
          finalPotential = ( trueVacuumPotential - pathFalsePotential );
          // At the end of the loop, segmentEndPotential is still the potential
          // at the end of the last straight segment.
          lastSegmentQuadratic = ( ( segmentEndPotential - finalPotential )
                           * inverseOfAuxiliaryStep * inverseOfAuxiliaryStep );
          numberOfNormalSegments = normalSegmentIndex;
          definiteOvershootAuxiliary
          = ( ( numberOfNormalSegments + 2 ) * auxiliaryStep );
        }
      }
      if( trueVacuumLowerThanPathFalseMinimum )
      {
        startOfFinalSegment = ( definiteOvershootAuxiliary - auxiliaryStep );
        if( !foundDefiniteUndershoot )
        {
          // If none of the straight segments crossed into negative potential,
          // then it happened in the last segment, so we find the auxiliary
          // difference from definiteOvershootAuxiliary which brings the
          // potential back up to 0.
          definiteUndershootAuxiliary = ( definiteOvershootAuxiliary
                            - sqrt( -finalPotential / lastSegmentQuadratic ) );
        }
      }
    }
  }

  SplinePotential::~SplinePotential()
  {
    // This does nothing.
  }


  // This returns the value of the potential at auxiliaryValue, by finding the
  // correct segment and then returning its value at that point.
  double
  SplinePotential::operator()( double const auxiliaryValue ) const
  {
    if( auxiliaryValue <= 0.0 )
    {
      return 0.0;
    }
    if( auxiliaryValue < auxiliaryStep )
    {
      return ( auxiliaryValue * auxiliaryValue * firstSegmentQuadratic );
    }
    if( auxiliaryValue >= definiteOvershootAuxiliary )
    {
      return finalPotential;
    }
    if( auxiliaryValue >= startOfFinalSegment )
    {
      return ( finalPotential
               + ( auxiliaryValue * auxiliaryValue * lastSegmentQuadratic ) );
    }
    size_t const auxiliarySteps( auxiliaryValue * inverseOfAuxiliaryStep );
    double const auxiliaryDifference( auxiliaryValue - auxiliarySteps );
    return ( potentialValues[ auxiliarySteps - 1 ]
          + ( auxiliaryDifference * firstDerivatives[ auxiliarySteps - 1 ] ) );
  }

  // This is for debugging.
  std::string SplinePotential::AsDebuggingString() const
  {
    BOL::StringParser doubleFormatter( 6,
                                       '',
                                       8,
                                       2,
                                       "",
                                       "-",
                                       "",
                                       "-",
                                       "*10^");
    std::stringstream returnStream;
    returnStream << "energyBarrierResolved = " << energyBarrierResolved
    << ", trueVacuumLowerThanPathFalseMinimum = "
    << trueVacuumLowerThanPathFalseMinimum
    << ", potential = ( UnitStep[x] * ( "
    << doubleFormatter.doubleToString( firstSegmentQuadratic )
    << " * x^(2) ) * UnitStep["
    << doubleFormatter.doubleToString( auxiliaryStep ) << " - x]";
    double cumulativeAuxiliary( auxiliaryStep );
    for( size_t segmentIndex( 0 );
         segmentIndex < numberOfNormalSegments;
         ++segmentIndex )
    {
      returnStream << " + "
      << "UnitStep[x - "
      << doubleFormatter.doubleToString( cumulativeAuxiliary )
      << "] * ( ("
      << doubleFormatter.doubleToString( potentialValues[ segmentIndex ] )
      << ") + (x-("
      << doubleFormatter.doubleToString( cumulativeAuxiliary )
      << ")) * ("
      << doubleFormatter.doubleToString( firstDerivatives[ segmentIndex ] )
      << ") ) * UnitStep[";
      cumulativeAuxiliary += auxiliaryStep;
      returnStream << doubleFormatter.doubleToString( cumulativeAuxiliary )
      << " - x]";
    }
    returnStream
    << " + UnitStep[x - "
    << doubleFormatter.doubleToString( cumulativeAuxiliary ) << "] * ( ("
    << doubleFormatter.doubleToString( finalPotential )
    << ") + (x-" << definiteOvershootAuxiliary << ")^2 * ("
    << doubleFormatter.doubleToString( lastSegmentQuadratic )
    << ") ) * UnitStep["
    << doubleFormatter.doubleToString( definiteOvershootAuxiliary )
    << " - x] )" << std::endl;
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
