/*
 * SplinePotential.cpp
 *
 *  Created on: Jun 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/TunnelingCalculation.hpp"

namespace VevaciousPlusPlus
{

  SplinePotential::SplinePotential(
                                   double const minimumFalseVacuumConcavity ) :
    auxiliaryValues(),
    potentialValues( 1,
                     0.0 ),
    firstDerivatives(),
    halfSecondDerivatives(),
    finalPotential( NAN ),
    halfFinalSecondDerivative( NAN ),
    finalCubicCoefficient( NAN ),
    minimumFalseVacuumConcavity( minimumFalseVacuumConcavity ),
    definiteUndershootAuxiliary( NAN ),
    definiteOvershootAuxiliary( 1.0 ),
    auxiliaryUpToCurrentSegment( 0.0 ),
    startOfFinalSegment( NAN ),
    sizeOfFinalSegment( NAN )
  {
    // This constructor is just an initialization list.
  }

  SplinePotential::SplinePotential( SplinePotential const& copySource ) :
    auxiliaryValues( copySource.auxiliaryValues ),
    potentialValues( copySource.potentialValues ),
    firstDerivatives( copySource.firstDerivatives ),
    halfSecondDerivatives( copySource.halfSecondDerivatives ),
    finalPotential( copySource.finalPotential ),
    halfFinalSecondDerivative( copySource.halfFinalSecondDerivative ),
    finalCubicCoefficient( copySource.finalCubicCoefficient ),
    minimumFalseVacuumConcavity( copySource.minimumFalseVacuumConcavity ),
    definiteUndershootAuxiliary( copySource.definiteUndershootAuxiliary ),
    definiteOvershootAuxiliary( copySource.definiteOvershootAuxiliary ),
    auxiliaryUpToCurrentSegment( copySource.auxiliaryUpToCurrentSegment ),
    startOfFinalSegment( copySource.startOfFinalSegment ),
    sizeOfFinalSegment( copySource.sizeOfFinalSegment )
  {
    // This constructor is just an initialization list.
  }

  SplinePotential::~SplinePotential()
  {
    // This does nothing.
  }


  // This returns the value of the potential at auxiliaryValue, by finding the
  // correct segment and then returning its value at that point.
  double SplinePotential::operator()( double const auxiliaryValue ) const
  {
    if( auxiliaryValue <= 0.0 )
    {
      return potentialValues.front();
    }
    if( auxiliaryValue >= definiteOvershootAuxiliary )
    {
      return finalPotential;
    }
    size_t segmentIndex( 0 );
    double auxiliaryDifference( auxiliaryValue );
    while( segmentIndex < auxiliaryValues.size() )
    {
      if( auxiliaryDifference < auxiliaryValues[ segmentIndex ] )
      {
        return ( potentialValues[ segmentIndex ]
                 + ( ( firstDerivatives[ segmentIndex ]
                       + ( halfSecondDerivatives[ segmentIndex ]
                           * auxiliaryDifference ) )
                     * auxiliaryDifference ) );
      }
      auxiliaryDifference -= auxiliaryValues[ segmentIndex ];
      ++segmentIndex;
    }
    // If we get to here, we're beyond the last normal segment:
    // auxiliaryValues.back() < auxiliaryValue < definiteOvershootAuxiliary.
    auxiliaryDifference = ( auxiliaryValue - definiteOvershootAuxiliary );
    return ( finalPotential
             + ( ( halfFinalSecondDerivative
                   + ( finalCubicCoefficient * auxiliaryDifference ) )
                 * auxiliaryDifference * auxiliaryDifference ) );
  }

  // This returns the value of the first derivative of the potential at
  // auxiliaryValue, by finding the correct segment and then returning its
  // slope at that point.
  double SplinePotential::FirstDerivative( double const auxiliaryValue ) const
  {
    if( ( auxiliaryValue <= 0.0 )
        ||
        ( auxiliaryValue >= definiteOvershootAuxiliary ) )
    {
      return 0.0;
    }
    size_t segmentIndex( 0 );
    double auxiliaryDifference( auxiliaryValue );
    while( segmentIndex < auxiliaryValues.size() )
    {
      if( auxiliaryDifference < auxiliaryValues[ segmentIndex ] )
      {
        return ( firstDerivatives[ segmentIndex ]
                 + ( 2.0 * halfSecondDerivatives[ segmentIndex ]
                         * auxiliaryDifference ) );
      }
      auxiliaryDifference -= auxiliaryValues[ segmentIndex ];
      ++segmentIndex;
    }
    // If we get to here, we're beyond the last normal segment:
    // auxiliaryValues[ currentGivenSize - 1 ] < auxiliaryValue < 1.0.
    auxiliaryDifference = ( auxiliaryValue - definiteOvershootAuxiliary );
    return ( ( ( 2.0 * halfFinalSecondDerivative )
               + ( 3.0 * finalCubicCoefficient * auxiliaryDifference ) )
                 * auxiliaryDifference );
  }

  // This returns the value of the first derivative of the potential at
  // auxiliaryValue, by finding the correct segment and then returning its
  // slope at that point.
  double SplinePotential::SecondDerivative( double const auxiliaryValue ) const
  {
    if( ( auxiliaryValue <= 0.0 )
        ||
        ( auxiliaryValue >= definiteOvershootAuxiliary ) )
    {
      return 0.0;
    }
    size_t segmentIndex( 0 );
    double auxiliaryDifference( auxiliaryValue );
    while( segmentIndex < auxiliaryValues.size() )
    {
      if( auxiliaryDifference < auxiliaryValues[ segmentIndex ] )
      {
        return ( 2.0 * halfSecondDerivatives[ segmentIndex ] );
      }
      auxiliaryDifference -= auxiliaryValues[ segmentIndex ];
      ++segmentIndex;
    }
    // If we get to here, we're beyond the last normal segment:
    // auxiliaryValues[ currentGivenSize - 1 ] < auxiliaryValue < 1.0.
    auxiliaryDifference = ( auxiliaryValue - definiteOvershootAuxiliary );
    return ( ( 2.0 * halfFinalSecondDerivative )
             + ( 3.0 * finalCubicCoefficient * auxiliaryDifference ) );
  }

  // This sets up the spline based on auxiliaryValues and potentialValues,
  // ensuring that the potential reaches the correct values for the vacua,
  // and that the potential derivative vanishes at the vacua. It also notes
  // the first point where the potential drops below that of the false vacuum
  // in definiteUndershootAuxiliary and the first maximum after that in
  // definiteOvershootAuxiliary, and cuts off the potential at that maximum.
  void
  SplinePotential::SetSpline( double const trueVacuumPotentialDifference )
  {
    firstDerivatives.resize( auxiliaryValues.size() );
    halfSecondDerivatives.resize( auxiliaryValues.size() );
    finalPotential = trueVacuumPotentialDifference;
    firstDerivatives[ 0 ] = 0.0;
    // There is a chance that numerical jitter means that the false vacuum may
    // be slightly ahead or behind the auxiliary variable p = 0. Either way we
    // flatten out the spline a bit (either forced concave by setting both
    // the slope to 0 and the second derivative to minimumFalseVacuumConcavity
    // if the proper false vacuum is at positive p, or just flattening it at
    // p = 0 if it should be at negative p, by setting the slope to 0 at
    // p = 0). It has the side-effect that potentials with energy barriers too
    // thin to be resolved by the path resolution are then forced to have a
    // small barrier at least.
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "potentialValues[ 0 ] = " << potentialValues[ 0 ]
    << ", potentialValues[ 1 ] = " << potentialValues[ 1 ];
    std::cout << std::endl;/**/
    if( potentialValues[ 1 ] < 0.0 )
    {
      halfSecondDerivatives[ 0 ] = minimumFalseVacuumConcavity;
      potentialValues[ 0 ]
      = ( potentialValues[ 1 ]
          - ( minimumFalseVacuumConcavity
              * auxiliaryValues[ 0 ] * auxiliaryValues[ 0 ] ) );
      finalPotential += potentialValues[ 0 ];
      // The true vacuum is also lowered to ensure that tunneling is still
      // possible.
    }
    else
    {
      halfSecondDerivatives[ 0 ] = ( potentialValues[ 1 ]
                           / ( auxiliaryValues[ 0 ] * auxiliaryValues[ 0 ] ) );
    }
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "auxiliaryValues[ 0 ] = " << auxiliaryValues[ 0 ]
    << ", potentialValues[ 0 ] = " << potentialValues[ 0 ]
    << ", firstDerivatives[ 0 ] = " << firstDerivatives[ 0 ]
    << ", halfSecondDerivatives[ 0 ] = " << halfSecondDerivatives[ 0 ];
    std::cout << std::endl;/**/
    bool definiteUndershootFound( false );
    bool definiteOvershootFound( false );
    auxiliaryUpToCurrentSegment = auxiliaryValues[ 0 ];
    for( size_t segmentIndex( 1 );
         segmentIndex < auxiliaryValues.size();
         ++segmentIndex )
    {
      firstDerivatives[ segmentIndex ]
      = ( firstDerivatives[ segmentIndex - 1 ]
          + ( 2.0 * auxiliaryValues[ segmentIndex - 1 ]
                  * halfSecondDerivatives[ segmentIndex - 1 ] ) );
      halfSecondDerivatives[ segmentIndex ]
      = ( ( potentialValues[ segmentIndex + 1 ]
            - potentialValues[ segmentIndex ]
            - ( firstDerivatives[ segmentIndex ]
                * auxiliaryValues[ segmentIndex ] ) )
        / ( auxiliaryValues[ segmentIndex ]
            * auxiliaryValues[ segmentIndex ] ) );

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "auxiliaryUpToCurrentSegment = " << auxiliaryUpToCurrentSegment
      << ", auxiliaryValues[ " << segmentIndex << " ] = "
      << auxiliaryValues[ segmentIndex ]
      << ", potentialValues[ " << segmentIndex << " ] = "
      << potentialValues[ segmentIndex ]
      << ", firstDerivatives[ " << segmentIndex << " ] = "
      << firstDerivatives[ segmentIndex ]
      << ", halfSecondDerivatives[ " << segmentIndex << " ] = "
      << halfSecondDerivatives[ segmentIndex ];
      std::cout << std::endl;/**/

      // Now we check for the potential dropping below potentialValues[ 0 ] in
      // this segment.
      if( !definiteUndershootFound )
      {
        if( ( halfSecondDerivatives[ segmentIndex ] == 0.0 )
            &&
            ( potentialValues[ segmentIndex + 1 ] < potentialValues[ 0 ] ) )
        {
          definiteUndershootAuxiliary
          = ( auxiliaryUpToCurrentSegment
              + ( ( potentialValues[ 0 ] - potentialValues[ segmentIndex ] )
                  / firstDerivatives[ segmentIndex ] ) );
          definiteUndershootFound = true;
        }
        else
        {
          double const
          discriminantValue( ( firstDerivatives[ segmentIndex ]
                               * firstDerivatives[ segmentIndex ] )
                             - ( 4.0 * potentialValues[ segmentIndex ]
                                   * halfSecondDerivatives[ segmentIndex ] ) );
          if( discriminantValue >= 0.0 )
          {
            double const discriminantRoot( sqrt( discriminantValue ) );
            double
            crossingAuxiliary( ( -0.5 * ( firstDerivatives[ segmentIndex ]
                                          + discriminantRoot ) )
                               / ( halfSecondDerivatives[ segmentIndex ] ) );
            // We take the lower value for crossing to lower potential first,
            // but switch to the upper value if the lower value is negative
            // (and thus outside the segment).
            if( crossingAuxiliary < 0.0 )
            {
              crossingAuxiliary += ( discriminantRoot
                                 / ( halfSecondDerivatives[ segmentIndex ] ) );
            }
            // Then we check that crossingAuxiliary is within the segment.
            if( ( crossingAuxiliary >= 0.0 )
                &&
                ( crossingAuxiliary < auxiliaryValues[ segmentIndex ] ) )
            {
              definiteUndershootAuxiliary
              = ( auxiliaryUpToCurrentSegment + crossingAuxiliary );
              definiteUndershootFound = true;
            }
          }
        }
        // debugging:
        /**/std::cout << std::endl << "debugging:"
        << std::endl
        << "after checking, definiteUndershootFound = "
        << definiteUndershootFound << ", definiteUndershootAuxiliary = "
        << definiteUndershootAuxiliary;
        std::cout << std::endl;/**/
      }
      // End of checking for crossing the line where the potential equals its
      // value at the false vacuum.

      // Now we check to see if we are looking for the 1st minimum after the
      // definite undershoot point.
      if( definiteUndershootFound
          &&
          !definiteOvershootFound )
      {
        double const
        extremumAuxiliary( ( -0.5 * firstDerivatives[ segmentIndex ] )
                           / ( halfSecondDerivatives[ segmentIndex ] ) );
        // debugging:
        /**/std::cout << std::endl << "debugging:"
        << std::endl
        << "definiteUndershootFound = " << definiteUndershootFound
        << ", definiteOvershootFound = " << definiteOvershootFound
        << ", extremumAuxiliary = " << extremumAuxiliary
        << ", definiteUndershootAuxiliary = " << definiteUndershootAuxiliary;
        std::cout << std::endl;/**/
        if( ( extremumAuxiliary > std::max( 0.0,
                                            ( definiteUndershootAuxiliary
                                            - auxiliaryUpToCurrentSegment ) ) )
            &&
            ( extremumAuxiliary < auxiliaryValues[ segmentIndex ] ) )
        {
          // If we find a minimum in a normal segment after having crossed to
          // deeper than the false vacuum (since the spline segments are only
          // quadratic, crossing the depth means that the slope is negative
          // and within any segment that starts with a negative slope, if there
          // is an extremum, it is a minimum; the slope also cannot turn
          // positive without having gone through a minimum), we set it up to
          // be the implicit final segment, purely as
          // finalPotential
          // + (definiteOvershootAuxiliary-p)^2 * halfFinalSecondDerivative
          // (no cubic term, and the linear term is absorbed by the shift of
          // the end of the segment to the minimum), and then remove this
          // normal segment and all those after it, and ends the function.
          definiteOvershootAuxiliary
          = ( auxiliaryUpToCurrentSegment + extremumAuxiliary );
          definiteOvershootFound = true;
          finalPotential = ( potentialValues[ segmentIndex ]
                             + ( firstDerivatives[ segmentIndex ]
                                 * extremumAuxiliary )
                             + ( halfSecondDerivatives[ segmentIndex ]
                                 * extremumAuxiliary * extremumAuxiliary ) );
          halfFinalSecondDerivative = halfSecondDerivatives[ segmentIndex ];
          finalCubicCoefficient = 0.0;
          auxiliaryValues.resize( segmentIndex );
          potentialValues.resize( segmentIndex );
          firstDerivatives.resize( segmentIndex );
          halfSecondDerivatives.resize( segmentIndex );
          startOfFinalSegment = auxiliaryUpToCurrentSegment;
          sizeOfFinalSegment = extremumAuxiliary;
          // debugging:
          /**/std::cout << std::endl << "debugging:"
          << std::endl
          << "found early path panic minimum in segment " << segmentIndex
          << ", startOfFinalSegment = " << startOfFinalSegment
          << ", sizeOfFinalSegment = " << sizeOfFinalSegment
          << ", definiteUndershootFound = " << definiteUndershootFound
          << ", definiteUndershootAuxiliary = " << definiteUndershootAuxiliary
          << ", definiteOvershootAuxiliary = " << definiteOvershootAuxiliary;
          std::cout << std::endl;/**/
          return;
        }
      }

      // Now we note the auxiliary value that starts the next segment.
      auxiliaryUpToCurrentSegment += auxiliaryValues[ segmentIndex ];
    }
    startOfFinalSegment = auxiliaryUpToCurrentSegment;

    if( !definiteUndershootFound )
    {
      // If the potential hasn't dropped below the value of the false potential
      // in the other segments, we'll just leave the definite undershoot
      // auxiliary value as the start of the final segment.
      definiteUndershootAuxiliary = auxiliaryUpToCurrentSegment;
    }

    // If we get to here, there is only 1 deeper-than-false-vacuum minimum and
    // it is in the implicit final segment. The last element of potentialValues
    // is already the value of the potential at the end of the last normal
    // segment, but firstDerivatives.back() is only the slope at the start of
    // the last normal segment, so the second derivative must be added,
    // multiplied by the last element of auxiliaryValues.
    sizeOfFinalSegment = ( 1.0 - auxiliaryUpToCurrentSegment );
    double const finalPotentialDrop( potentialValues.back() - finalPotential );
    double const finalStartSlope( firstDerivatives.back()
                                  + ( 2.0 * halfSecondDerivatives.back()
                                          * auxiliaryValues.back() ) );
    // If the slope at the start of the final segment is steep enough that a
    // quadratic spline would have a path panic minimum before p = 1.0, we keep
    // it and just end the spline at that minimum (no cubic term, and the
    // linear term is absorbed by the shift of the end of the segment to the
    // minimum). This is definitely the case if finalPotentialDrop is negative
    // or small enough compared to -finalStartSlope (finalStartSlope must be
    // negative or we would have already found an early path potential
    // minimum).
    if( ( 2.0 * finalPotentialDrop )
        + ( finalStartSlope * sizeOfFinalSegment ) <= 0.0 )
    {
      finalCubicCoefficient = 0.0;
      halfFinalSecondDerivative
      = ( -( finalPotentialDrop + ( finalStartSlope * sizeOfFinalSegment ) )
          / ( sizeOfFinalSegment * sizeOfFinalSegment ) );
      double const
      offsetOfMinimumFromOne( 0.5 * ( sizeOfFinalSegment
                                      - ( finalPotentialDrop
                                          / ( halfFinalSecondDerivative
                                              * sizeOfFinalSegment ) ) ) );
      definiteOvershootAuxiliary = ( 1.0 - offsetOfMinimumFromOne );
      sizeOfFinalSegment -= offsetOfMinimumFromOne;
      finalPotential = ( potentialValues.back()
                         - ( halfFinalSecondDerivative * sizeOfFinalSegment
                                                      * sizeOfFinalSegment ) );
    }
    else
    {
      halfFinalSecondDerivative = ( ( 3.0 * finalPotentialDrop )
                                  + ( finalStartSlope * sizeOfFinalSegment ) );
      finalCubicCoefficient
      = ( ( halfFinalSecondDerivative - finalPotentialDrop )
          / ( sizeOfFinalSegment * sizeOfFinalSegment * sizeOfFinalSegment ) );
      halfFinalSecondDerivative = ( halfFinalSecondDerivative
                               / ( sizeOfFinalSegment * sizeOfFinalSegment ) );
      definiteOvershootAuxiliary = 1.0;
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "Implicit final segment:" << std::endl
    << "startOfFinalSegment = " << startOfFinalSegment
    << ", sizeOfFinalSegment = " << sizeOfFinalSegment
    << ", finalPotential = " << finalPotential
    << ", halfFinalSecondDerivative = " << halfFinalSecondDerivative
    << ", finalCubicCoefficient = " << finalCubicCoefficient
    << ", definiteUndershootAuxiliary = " << definiteUndershootAuxiliary
    << ", definiteOvershootAuxiliary = " << definiteOvershootAuxiliary;
    std::cout << std::endl;/**/
  }

  // This is for debugging.
  std::string SplinePotential::AsDebuggingString() const
  {
    std::stringstream returnStream;
    double cumulativeAuxiliary( 0.0 );
    for( size_t segmentIndex( 0 );
         segmentIndex < auxiliaryValues.size();
         ++segmentIndex )
    {
      if( segmentIndex > 0 )
      {
        returnStream << " + ";
      }
      returnStream
      << "UnitStep[x - " << cumulativeAuxiliary << "] * ( ("
      << potentialValues[ segmentIndex ] << ") + (x-(" << cumulativeAuxiliary
      << ")) * (" << firstDerivatives[ segmentIndex ] << ") + (x-("
      << cumulativeAuxiliary << "))^2 * ("
      << halfSecondDerivatives[ segmentIndex ] << ") ) * UnitStep[";
      cumulativeAuxiliary += auxiliaryValues[ segmentIndex ];
      returnStream << cumulativeAuxiliary << " - x]" << std::endl;
    }
    returnStream
    << " + UnitStep[x - " << cumulativeAuxiliary << "] * ( (" << finalPotential
    << ") + (x-" << definiteOvershootAuxiliary << ")^2 * ("
    << halfFinalSecondDerivative << ") + (x-" << definiteOvershootAuxiliary
    << ")^3 * (" << finalCubicCoefficient << ") ) * UnitStep["
    << definiteOvershootAuxiliary << " - x]";
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
