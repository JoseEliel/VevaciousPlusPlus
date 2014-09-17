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
                                    size_t const numberOfPotentialSegments ) :
    auxiliaryStep( 1.0 / static_cast< double >( numberOfPotentialSegments ) ),
    inverseOfAuxiliaryStep( numberOfPotentialSegments ),
    potentialValues( numberOfPotentialSegments - 2 ),
    firstDerivatives( numberOfPotentialSegments - 2 ),
    firstSegmentQuadratic( NAN ),
    finalPotential( NAN ),
    lastSegmentQuadratic( NAN ),
    definiteUndershootAuxiliary( NAN ),
    definiteOvershootAuxiliary( 1.0 ),
    startOfFinalSegment( NAN )
  {
    std::vector< double >
    fieldConfiguration( potentialFunction.NumberOfFieldVariables() );
    tunnelPath.PutOnPathAt( fieldConfiguration,
                            0.0 );
    double const pathTemperature( tunnelPath.TemperatureValue() );
    double falseVacuumPotential( potentialFunction( fieldConfiguration,
                                                    pathTemperature ) );
    tunnelPath.PutOnPathAt( fieldConfiguration,
                            auxiliaryStep );
    bool stillLookingForFalseVacuum( true );
    size_t skippedFalseSideSegments( 0 );
    bool foundPathPanicVacuum( false );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "SplinePotential::SplinePotential( ..., tunnelPath =" << std::endl;
    std::cout << tunnelPath.AsDebuggingString() << std::endl;
    std::cout << ", numberOfPotentialSegments = " << numberOfPotentialSegments
    << " called ). falseVacuumPotential = " << falseVacuumPotential;
    std::cout << std::endl;/**/

    double segmentEndPotential( NAN );
    for( size_t nodeIndex( 0 );
         nodeIndex < ( numberOfPotentialSegments - 3 );
         ++nodeIndex )
    {
      tunnelPath.PutOnPathAt( fieldConfiguration,
                              ( ( nodeIndex + 1 ) * auxiliaryStep ) );
      segmentEndPotential = ( potentialFunction( fieldConfiguration,
                                                 pathTemperature )
                              - falseVacuumPotential );

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "nodeIndex = " << nodeIndex << ", falseVacuumPotential = "
      << falseVacuumPotential << ", segmentEndPotential = "
      << segmentEndPotential;
      std::cout << std::endl;/**/

      // If we have not found the start of an energy barrier yet...
      if( stillLookingForFalseVacuum )
      {
        // ... we check to see if we have found an energy barrier.
        if( segmentEndPotential > 0.0 )
        {
          skippedFalseSideSegments = nodeIndex;
          firstSegmentQuadratic = ( segmentEndPotential
                           * inverseOfAuxiliaryStep * inverseOfAuxiliaryStep );
          potentialValues.front() = segmentEndPotential;
          stillLookingForFalseVacuum = false;
        }
        else
        {
          falseVacuumPotential = segmentEndPotential;
        }
      }
      else
      {
        // ... otherwise we start mapping out the potential.
        size_t const
        currentSegmentIndex( nodeIndex - skippedFalseSideSegments );
        potentialValues[ currentSegmentIndex + 1 ]
        = segmentEndPotential;
        firstDerivatives[ currentSegmentIndex ]
        = ( ( segmentEndPotential - potentialValues[ currentSegmentIndex ] )
            * inverseOfAuxiliaryStep );
        // Now we finish early if we've found the path panic vacuum.
        if( ( potentialValues[ currentSegmentIndex ] < 0.0 )
            &&
            ( firstDerivatives[ currentSegmentIndex ] > 0.0 ) )
        {
          finalPotential = potentialValues[ currentSegmentIndex ];
          lastSegmentQuadratic = ( potentialValues[ currentSegmentIndex - 1 ]
                           * inverseOfAuxiliaryStep * inverseOfAuxiliaryStep );
          potentialValues.resize( currentSegmentIndex - 1 );
          firstDerivatives.resize( currentSegmentIndex - 1 );
          foundPathPanicVacuum = true;
          break;
        }
      }
    }
    if( !foundPathPanicVacuum )
    {
      // If we did not break early because of an early path panic minimum, we
      // still have to check if the (implicit) final segment for an early panic
      // minimum.
      tunnelPath.PutOnPathAt( fieldConfiguration,
                              1.0 );
      finalPotential = ( potentialFunction( fieldConfiguration,
                                            pathTemperature )
                        - falseVacuumPotential );
      if( finalPotential < segmentEndPotential )
      {
        lastSegmentQuadratic = ( ( segmentEndPotential - finalPotential )
                           * inverseOfAuxiliaryStep * inverseOfAuxiliaryStep );
      }
      else
      {
        finalPotential = segmentEndPotential;
        lastSegmentQuadratic = ( ( potentialValues.back() - finalPotential )
                           * inverseOfAuxiliaryStep * inverseOfAuxiliaryStep );
        potentialValues.resize( potentialValues.size() - 1 );
        firstDerivatives.resize( firstDerivatives.size() - 1 );
      }
    }
    definiteOvershootAuxiliary = ( ( potentialValues.size() + 2 )
                                   * auxiliaryStep );
    startOfFinalSegment = ( definiteOvershootAuxiliary - auxiliaryStep );
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
    returnStream << "UnitStep[x] * ( "
    << doubleFormatter.doubleToString( firstSegmentQuadratic )
    << " * x^(2) ) * UnitStep["
    << doubleFormatter.doubleToString( auxiliaryStep ) << " - x]";
    double cumulativeAuxiliary( auxiliaryStep );
    for( size_t segmentIndex( 0 );
         segmentIndex < potentialValues.size();
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
      << " - x]" << std::endl;
    }
    returnStream
    << " + UnitStep[x - "
    << doubleFormatter.doubleToString( cumulativeAuxiliary ) << "] * ( ("
    << doubleFormatter.doubleToString( finalPotential )
    << ") + (x-" << definiteOvershootAuxiliary << ")^2 * ("
    << doubleFormatter.doubleToString( lastSegmentQuadratic )
    << ") ) * UnitStep["
    << doubleFormatter.doubleToString( definiteOvershootAuxiliary ) << " - x]";
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
