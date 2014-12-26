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
    auxiliaryStep( 1.0 / static_cast< double >( numberOfPotentialSegments ) ),
    inverseOfAuxiliaryStep( numberOfPotentialSegments ),
    numberOfNormalSegments( numberOfPotentialSegments - 2 ),
    potentialValues( numberOfNormalSegments,
                     0.0 ),
    firstDerivatives( numberOfNormalSegments,
                      0.0 ),
    pathFalsePotential( -1.0 ),
    firstSegmentQuadratic( -1.0 ),
    finalPotential( 0.0 ),
    lastSegmentQuadratic( -1.0 ),
    energyBarrierWasResolved( false ),
    tunnelingPossibleOnPath( false )
  {
    // First we have to find the path false minimum. The base constructor
    // already set auxiliaryOfPathPanicVacuum to zero.
    std::vector< double >
    fieldConfiguration( potentialFunction.NumberOfFieldVariables() );
    tunnelPath.PutOnPathAt( fieldConfiguration,
                            auxiliaryOfPathFalseVacuum );
    std::vector< double > const
    pathFalseEndConfiguration( fieldConfiguration );
    double const pathTemperature( tunnelPath.TemperatureValue() );
    pathFalsePotential = potentialFunction( pathFalseEndConfiguration,
                                            pathTemperature );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "Path start: "
    << potentialFunction.FieldConfigurationAsMathematica(
                                                    pathFalseEndConfiguration )
    << ", potential value = " << pathFalsePotential;
    std::cout << std::endl;/**/

    // We start at the assumed path false vacuum, then move along the resolved
    // segment ends to find any that are lower which are not separated by an
    // energy barrier which ends sufficiently far away from the given path
    // false vacuum end. The point where the energy barrier "ends" is where the
    // potential drops below what it is for the path false vacuum, and this
    // point must be separated in field space from the false vacuum end of the
    // path by at least sqrt(minimumSquareDistanceBetweenPathVacua).
    while( !energyBarrierWasResolved )
    {
      potentialValues.resize( 1,
                              -1.0 );
      // The return value of RollForwardToLocalMinimum is false if it moves all
      // the way to the true vacuum end of the path without the potential ever
      // increasing. Here fieldConfiguration is being used as the field
      // configuration of the path false vacuum.
      tunnelingPossibleOnPath = RollToPathLocalMinimum( potentialFunction,
                                                        tunnelPath,
                                                        pathTemperature,
                                                        fieldConfiguration );
      if( !tunnelingPossibleOnPath )
      {
        // There is no point in continuing if tunneling is not possible.
        break;
      }

      // The return value of FindEndOfEnergyBarrier is false if it moves all
      // the way to the true vacuum end of the path without the potential ever
      // dropping below pathFalsePotential. Here fieldConfiguration is being
      // used as the field configuration of the point at the end of the energy
      // barrier.
      tunnelingPossibleOnPath = FindPathPanicVacuum( potentialFunction,
                                                     tunnelPath,
                                                     pathTemperature,
                                                     fieldConfiguration );

      if( !tunnelingPossibleOnPath )
      {
        // There is no point in continuing if tunneling is not possible.
        break;
      }

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
      // potential has dropped below pathFalsePotential by
      // definiteUndershootAuxiliary, and the path panic vacuum is at
      // auxiliaryOfPathPanicVacuum along tunnelPath and is sufficiently
      // separated from the false vacuum end of tunnelPath. Now we have to
      // ensure that the parabolic segments have the proper values.
      firstSegmentQuadratic = ( potentialValues.front()
                                * inverseOfAuxiliaryStep
                                * inverseOfAuxiliaryStep );
      thresholdForNearPathPanic
      = ( auxiliaryOfPathPanicVacuum - auxiliaryStep );
      finalPotential = potentialValues.back();
      size_t const numberOfSegments( potentialValues.size() );
      for( size_t segmentIndex( 0 );
           segmentIndex < ( numberOfSegments - 1 );
           ++segmentIndex )
      {
        firstDerivatives[ segmentIndex ]
        = ( potentialValues[ segmentIndex + 1 ]
            - potentialValues[ segmentIndex ] );
      }

      lastSegmentQuadratic
      = ( ( potentialValues[ numberOfSegments - 2 ] - finalPotential )
          * inverseOfAuxiliaryStep * inverseOfAuxiliaryStep );
    }


    THIS IS WHERE THINGS HAVE TO BE WORKED OUT!


    tunnelPath.PutOnPathAt( pathFalseEndConfiguration,
                            1.0 );
    double const pathEndPotential( potentialFunction( pathFalseEndConfiguration,
                                                      pathTemperature ) );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "Path end: "
    << potentialFunction.FieldConfigurationAsMathematica( pathFalseEndConfiguration )
    << ", potential value = " << pathEndPotential;
    std::cout << std::endl;/**/

    double segmentEndAuxiliary( auxiliaryStep );
    tunnelPath.PutOnPathAt( pathFalseEndConfiguration,
                            segmentEndAuxiliary );
    double segmentEndPotential( potentialFunction( pathFalseEndConfiguration,
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
      tunnelPath.PutOnPathAt( pathFalseEndConfiguration,
                              segmentEndAuxiliary );
      segmentEndPotential = potentialFunction( pathFalseEndConfiguration,
                                               pathTemperature );

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "segmentIndex = " << segmentIndex << ", segmentEndAuxiliary = "
      << segmentEndAuxiliary << ", field configuration at segment end = "
      << potentialFunction.FieldConfigurationAsMathematica(
                                                           pathFalseEndConfiguration )
      << ", potential value = " << pathEndPotential;
      std::cout << std::endl;/**/
    }
    // Now we move on to resolving the energy barrier, looking out for an early
    // path panic minimum. However, this is only done if we resolved the start
    // of an energy barrier.
    if( energyBarrierResolved )
    {
      ++segmentIndex;
      segmentEndAuxiliary += auxiliaryStep;
      // Now segmentEndAuxiliary and segmentIndex are correct for the first
      // straight segment.
      size_t normalSegmentIndex( 0 );
      potentialValues.front() = ( segmentEndPotential - pathFalsePotential );
      tunnelPath.PutOnPathAt( pathFalseEndConfiguration,
                              segmentEndAuxiliary );
      // From this point on, pathFalsePotential is always subtracted from
      // potentials along the path.
      segmentEndPotential
      = ( potentialFunction( pathFalseEndConfiguration,
                             pathTemperature ) - pathFalsePotential );
      firstDerivatives.front() = ( inverseOfAuxiliaryStep *
                           ( segmentEndPotential - potentialValues.front() ) );
      bool foundDefiniteUndershoot( false );

      // We have to check whether the first straight segment goes below the
      // path false minimum.
      if( segmentEndPotential < 0 )
      {
        foundDefiniteUndershoot = true;
        definiteUndershootAuxiliary
        = ( auxiliaryStep * ( 1.0 - ( potentialValues.front()
                                      / firstDerivatives.front() ) ) );
      }

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "First straight segment: segmentIndex = " << segmentIndex
      << ", segmentEndAuxiliary = " << segmentEndAuxiliary
      << ", field configuration at segment end = "
      << potentialFunction.FieldConfigurationAsMathematica(
                                                           pathFalseEndConfiguration )
      << ", potential value = " << ( segmentEndPotential + pathFalsePotential )
      << ", relative value = " << segmentEndPotential;
      std::cout << std::endl;/**/

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "Why do I not note that trueVacuumLowerThanPathFalseMinimum must be"
      << " true when checking for the definite undershoot?" << std::endl
      << "Also, what happens if the given true vacuum is not as deep as the"
      << " point 1 segment back towards the false vacuum?";
      std::cout << std::endl;/**/

      while( segmentIndex < ( numberOfPotentialSegments - 2 ) )
      {
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
        tunnelPath.PutOnPathAt( pathFalseEndConfiguration,
                                segmentEndAuxiliary );
        // From this point on, pathFalsePotential is always subtracted from
        // potentials along the path.
        segmentEndPotential
        = ( potentialFunction( pathFalseEndConfiguration,
                               pathTemperature ) - pathFalsePotential );
        firstDerivatives[ normalSegmentIndex ]
        = ( inverseOfAuxiliaryStep
           * ( segmentEndPotential - potentialValues[ normalSegmentIndex ] ) );

        // We have to check whether this straight segment goes below the path
        // false minimum.
        if( !foundDefiniteUndershoot
            &&
            ( segmentEndPotential < 0 ) )
        {
          foundDefiniteUndershoot = true;
          definiteUndershootAuxiliary
          = ( segmentEndAuxiliary + ( potentialValues[ normalSegmentIndex ]
                                  / firstDerivatives[ normalSegmentIndex ] ) );
          // It is plus because potentialValues[ normalSegmentIndex ] must be
          // positive and firstDerivatives[ normalSegmentIndex ] must be
          // negative for this to happen, and thus definiteUndershootAuxiliary
          // will be less than segmentEndAuxiliary.
        }

        // debugging:
        /**/std::cout << std::endl << "debugging:"
        << std::endl
        << "segmentIndex = " << segmentIndex << ", normalSegmentIndex = "
        << normalSegmentIndex << ", segmentEndAuxiliary = "
        << segmentEndAuxiliary
        << ", field configuration at segment end = "
        << potentialFunction.FieldConfigurationAsMathematica(
                                                           pathFalseEndConfiguration )
        << ", potential value = "
        << ( segmentEndPotential + pathFalsePotential )
        << ", relative value = " << segmentEndPotential
        << ", potentialValues[ " << normalSegmentIndex << " ] = "
        << potentialValues[ normalSegmentIndex ] << ", firstDerivatives[ "
        << normalSegmentIndex << " ] = "
        << firstDerivatives[ normalSegmentIndex ];
        std::cout << std::endl;/**/
      }
      // Now we have to check to see if the end of the path is lower than the
      // path false minimum, if we haven't already found an early path panic
      // minimum.
      if( !trueVacuumLowerThanPathFalseMinimum )
      {
        if( pathEndPotential < pathFalsePotential )
        {
          trueVacuumLowerThanPathFalseMinimum = true;
          finalPotential = ( pathEndPotential - pathFalsePotential );
          // At the end of the loop, segmentEndPotential is still the potential
          // at the end of the last straight segment.
          lastSegmentQuadratic = ( ( segmentEndPotential - finalPotential )
                           * inverseOfAuxiliaryStep * inverseOfAuxiliaryStep );
          numberOfNormalSegments = ( normalSegmentIndex + 1 );
          // After the loop ends without an early path panic minimum,
          // normalSegmentIndex is the index of the last normal segment, and
          // since the index starts at 0, the total number of segments is one
          // larger.
          definiteOvershootAuxiliary
          = ( ( numberOfNormalSegments + 2 ) * auxiliaryStep );
        }
      }
      if( trueVacuumLowerThanPathFalseMinimum )
      {
        thresholdForNearPathPanic
        = ( definiteOvershootAuxiliary - auxiliaryStep );
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
  double SplinePotential::operator()( double auxiliaryValue ) const
  {
    if( auxiliaryValue >= definiteOvershootAuxiliary )
    {
      return finalPotential;
    }
    if( auxiliaryValue >= thresholdForNearPathPanic )
    {
      double const differenceFromMaximumAuxiliary( definiteOvershootAuxiliary
                                                   - auxiliaryValue );
      return ( finalPotential + ( differenceFromMaximumAuxiliary
                  * differenceFromMaximumAuxiliary * lastSegmentQuadratic ) );
    }
    if( auxiliaryValue <= auxiliaryOfPathFalseVacuum )
    {
      return 0.0;
    }
    auxiliaryValue -= auxiliaryOfPathFalseVacuum;
    if( auxiliaryValue < auxiliaryStep )
    {
      return ( auxiliaryValue * auxiliaryValue * firstSegmentQuadratic );
    }
    size_t const auxiliarySteps( auxiliaryValue * inverseOfAuxiliaryStep );
    double const auxiliaryDifference( auxiliaryValue
                                      - ( auxiliarySteps * auxiliaryStep ) );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "SplinePotential::operator()( auxiliaryValue = " << auxiliaryValue
    << " ) called. auxiliarySteps = " << auxiliarySteps
    << ", auxiliaryDifference = " << auxiliaryDifference
    << ", potentialValues[ " << ( auxiliarySteps - 1 ) << " ] = "
    << potentialValues[ auxiliarySteps - 1 ]
    << ", firstDerivatives[ " << ( auxiliarySteps - 1 ) << " ] = "
    << firstDerivatives[ auxiliarySteps - 1 ];
    std::cout << std::endl;*/

    return ( potentialValues[ auxiliarySteps - 1 ]
          + ( auxiliaryDifference * firstDerivatives[ auxiliarySteps - 1 ] ) );
  }

  // This is for debugging.
  std::string SplinePotential::AsDebuggingString() const
  {
    BOL::StringParser doubleFormatter( 6,
                                       ' ',
                                       8,
                                       2,
                                       "",
                                       "-",
                                       "",
                                       "-",
                                       "*10^");
    std::stringstream returnStream;
    returnStream << "numberOfNormalSegments = " << numberOfNormalSegments
    << ", firstSegmentQuadratic = " << firstSegmentQuadratic
    << ", finalPotential = " << finalPotential
    << ", lastSegmentQuadratic = " << lastSegmentQuadratic
    << ", definiteUndershootAuxiliary = " << definiteUndershootAuxiliary
    << ", definiteOvershootAuxiliary = " << definiteOvershootAuxiliary
    << ", thresholdForNearPathPanic = " << thresholdForNearPathPanic
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

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "preparing points for Mathematica.";
    returnStream << "Mathematica points = { ";
    for( size_t sampleIndex( 0 );
         sampleIndex <= ( 3 * ( numberOfNormalSegments + 2 ) );
         ++sampleIndex )
    {
      if( sampleIndex > 0 )
      {
        returnStream << ", ";
      }
      double const sampleAuxiliary( ( sampleIndex * auxiliaryStep ) / 3.0 );
      returnStream << "{ " << sampleAuxiliary << ", "
      << (*this)( sampleAuxiliary ) << " }";
    }
    returnStream << " }" << std::endl;
    std::cout << std::endl;*/

    return returnStream.str();
  }


  // This goes along the segment ends checking for the lowest before the
  // potential starts to rise again, returning the auxiliary value of the
  // segment end with lowest potential, leaving the lowest potential in
  // pathFalsePotential and its field configuration in fieldConfiguration.
  double SplinePotential::FindPathFalseVacuum(
                                    double const pathFalseVacuumAuxiliarySoFar,
                              double& pathFalsePotential,
                              TunnelPath const& tunnelPath,
                              std::vector< double >& fieldConfiguration,
                        double const minimumSquareDistanceBetweenPathVacua );

} /* namespace VevaciousPlusPlus */
