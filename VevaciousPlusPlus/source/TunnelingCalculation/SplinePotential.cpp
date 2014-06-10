/*
 * SplinePotential.cpp
 *
 *  Created on: Jun 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/TunnelingCalculation.hpp"

namespace VevaciousPlusPlus
{

  SplinePotential::SplinePotential() :
    auxiliaryValues(),
    potentialValues( 1,
                     0.0 ),
    firstDerivatives(),
    halfSecondDerivatives(),
    finalPotential( NAN ),
    halfFinalSecondDerivative( NAN ),
    finalCubicCoefficient( NAN )
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
    finalCubicCoefficient( copySource.finalCubicCoefficient )
  {
    // This constructor is just an initialization list.
  }

  SplinePotential::~SplinePotential()
  {
    // This does nothing.
  }


  // This returns the value of the potential at auxiliaryValue, by finding
  // the correct spline and then returning its value at that point.
  double SplinePotential::operator()( double const auxiliaryValue ) const
  {
    if( auxiliaryValue <= 0.0 )
    {
      return 0.0;
    }
    if( auxiliaryValue >= 1.0 )
    {
      return finalPotential;
    }
    size_t splineIndex( 0 );
    double auxiliaryDifference( auxiliaryValue );
    while( splineIndex < auxiliaryValues.size() )
    {
      if( auxiliaryDifference < auxiliaryValues[ splineIndex ] )
      {
        return ( potentialValues[ splineIndex ]
                 + ( ( firstDerivatives[ splineIndex ]
                       + ( halfSecondDerivatives[ splineIndex ]
                           * auxiliaryDifference ) )
                     * auxiliaryDifference ) );
      }
      auxiliaryDifference -= auxiliaryValues[ splineIndex ];
      ++splineIndex;
    }
    // If we get to here, we're beyond the last normal spline:
    // auxiliaryValues[ currentGivenSize - 1 ] < auxiliaryValue < 1.0.
    auxiliaryDifference = ( auxiliaryValue - 1.0 );
    return ( finalPotential
             + ( ( halfFinalSecondDerivative
                   + ( finalCubicCoefficient * auxiliaryDifference ) )
                 * ( auxiliaryDifference * auxiliaryDifference ) ) );
  }

  // This returns the value of the first derivative of the potential at
  // auxiliaryValue, by finding the correct spline and then returning its slope
  // at that point.
  double SplinePotential::FirstDerivative( double const auxiliaryValue ) const
  {
    if( ( auxiliaryValue <= 0.0 )
        ||
        ( auxiliaryValue >= 1.0 ) )
    {
      return 0.0;
    }
    size_t splineIndex( 0 );
    double auxiliaryDifference( auxiliaryValue );
    while( splineIndex < auxiliaryValues.size() )
    {
      if( auxiliaryDifference < auxiliaryValues[ splineIndex ] )
      {
        return ( firstDerivatives[ splineIndex ]
                 + ( 2.0 * halfSecondDerivatives[ splineIndex ]
                         * auxiliaryDifference ) );
      }
      auxiliaryDifference -= auxiliaryValues[ splineIndex ];
      ++splineIndex;
    }
    // If we get to here, we're beyond the last normal spline:
    // auxiliaryValues[ currentGivenSize - 1 ] < auxiliaryValue < 1.0.
    auxiliaryDifference = ( 1.0 - auxiliaryValue );
    return ( ( ( 2.0 * halfFinalSecondDerivative )
               + ( 3.0 * finalCubicCoefficient * auxiliaryDifference ) )
             * auxiliaryDifference );
  }

  // This sets up the splines based on auxiliaryValues and potentialValues,
  // ensuring that the potential reaches the correct values for the vacua,
  // and that the potential derivative vanishes at the vacua.
  void
  SplinePotential::SetSplines( double const trueVacuumPotentialDifference )
  {
    firstDerivatives.resize( auxiliaryValues.size() );
    halfSecondDerivatives.resize( auxiliaryValues.size() );
    finalPotential = trueVacuumPotentialDifference;
    firstDerivatives[ 0 ] = 0.0;
    halfSecondDerivatives[ 0 ] = ( potentialValues[ 1 ]
                           / ( auxiliaryValues[ 0 ] * auxiliaryValues[ 0 ] ) );
    double finalDifference( 1.0 - auxiliaryValues[ 0 ] );
    for( size_t splineIndex( 1 );
         splineIndex < auxiliaryValues.size();
         ++splineIndex )
    {
      firstDerivatives[ splineIndex ]
      = ( firstDerivatives[ splineIndex - 1 ]
          + ( 2.0 * auxiliaryValues[ splineIndex - 1 ]
                  * halfSecondDerivatives[ splineIndex - 1 ] ) );
      halfSecondDerivatives[ splineIndex ]
      = ( ( potentialValues[ splineIndex + 1 ]
            - potentialValues[ splineIndex ]
            - ( firstDerivatives[ splineIndex ]
                * auxiliaryValues[ splineIndex ] ) )
        / ( auxiliaryValues[ splineIndex ]
            * auxiliaryValues[ splineIndex ] ) );
      finalDifference -= auxiliaryValues[ splineIndex ];
    }
    // The last element of potentialValues is already the value of the
    // potential at the end of the last normal spline, but
    // firstDerivatives.back() is only the slope at the start of the last
    // normal spline, so the second derivative must be added, multiplied by the
    // last element of auxiliaryValues.
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "finalDifference = " << finalDifference
    << ", finalPotential = " << finalPotential;
    std::cout << std::endl;/**/
    halfFinalSecondDerivative
    = ( ( 3.0 * ( potentialValues.back() - finalPotential ) )
        - ( ( firstDerivatives.back()
              + ( 2.0 * halfSecondDerivatives.back()
                      * auxiliaryValues.back() ) )
            * finalDifference ) );
    finalCubicCoefficient
    = ( ( potentialValues.back() - finalPotential
          - halfFinalSecondDerivative )
        / ( finalDifference * finalDifference * finalDifference ) );
    halfFinalSecondDerivative = ( halfFinalSecondDerivative
                                   / ( finalDifference * finalDifference ) );
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "THIS IS STILL WRONG! Something is wrong with the final coefficients.";
    std::cout << std::endl;/**/
  }

  // This is for debugging.
  std::string SplinePotential::AsDebuggingString() const
  {
    std::stringstream returnStream;
    double cumulativeAuxiliary( 0.0 );
    for( size_t splineIndex( 0 );
         splineIndex < auxiliaryValues.size();
         ++splineIndex )
    {
      if( splineIndex > 0 )
      {
        returnStream << " + ";
      }
      returnStream
      << "UnitStep[x - " << cumulativeAuxiliary << "] * ( ("
      << potentialValues[ splineIndex ] << ") + (x-(" << cumulativeAuxiliary
      << ")) * (" << firstDerivatives[ splineIndex ] << ") + (x-("
      << cumulativeAuxiliary << "))^2 * ("
      << halfSecondDerivatives[ splineIndex ] << ") ) * UnitStep[";
      cumulativeAuxiliary += auxiliaryValues[ splineIndex ];
      returnStream << cumulativeAuxiliary << " - x]" << std::endl;
    }
    returnStream
    << " + UnitStep[x - " << cumulativeAuxiliary << "] * ( (" << finalPotential
    << ") + (x-1)^2 * (" << halfFinalSecondDerivative << ") + (x-1)^3 * ("
    << finalCubicCoefficient << ") ) * UnitStep[1 - x]";
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
