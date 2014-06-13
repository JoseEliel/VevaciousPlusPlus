/*
 * OdeintBubbleDerivatives.cpp
 *
 *  Created on: Jun 5, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/TunnelingCalculation.hpp"

namespace VevaciousPlusPlus
{

  OdeintBubbleDerivatives::OdeintBubbleDerivatives(
                       PathFieldsAndPotential const& pathFieldsAndPotential ) :
    potentialSpline( pathFieldsAndPotential.PotentialApproximation() ),
    firstDerivatives( pathFieldsAndPotential.FieldDerivatives() ),
    numberOfFields( pathFieldsAndPotential.FieldDerivatives().size() ),
    secondDerivatives( firstDerivatives.size() ),
    dampingFactor( pathFieldsAndPotential.DampingFactor() )
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "OdeintBubbleDerivatives::OdeintBubbleDerivatives("
    << " pathFieldsAndPotential ="
    << pathFieldsAndPotential.AsDebuggingString() << " ) called.";
    std::cout << std::endl << "potentialSpline =" << std::endl
    << potentialSpline.AsDebuggingString();
    std::cout << std::endl;/**/

    for( unsigned int fieldIndex( 0 );
         fieldIndex < firstDerivatives.size();
         ++fieldIndex )
    {
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "firstDerivatives[ " << fieldIndex << " ] = "
      << firstDerivatives[ fieldIndex ].AsDebuggingString();
      std::cout << std::endl;/**/
      secondDerivatives[ fieldIndex ].BecomeFirstDerivativeOf(
                                              firstDerivatives[ fieldIndex ] );
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "secondDerivatives[ " << fieldIndex << " ] = "
      << secondDerivatives[ fieldIndex ].AsDebuggingString();
      std::cout << std::endl;/**/
    }
  }

  OdeintBubbleDerivatives::~OdeintBubbleDerivatives()
  {
    // This does nothing.
  }


  // This is in the form required for the Boost odeint package.
  void OdeintBubbleDerivatives::operator()(
                      std::vector< double > const& auxiliaryAndFirstDerivative,
                              std::vector< double >& firstAndSecondDerivatives,
                                            double const radialValue )
  {
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "OdeintBubbleDerivatives::operator( auxiliaryAndFirstDerivative = { ";
    for( std::vector< double >::const_iterator
         inputData( auxiliaryAndFirstDerivative.begin() );
         inputData < auxiliaryAndFirstDerivative.end();
         ++inputData )
    {
      if( inputData > auxiliaryAndFirstDerivative.begin() )
      {
        std::cout << ", ";
      }
      std::cout << *inputData;
    }
    std::cout << " }, firstAndSecondDerivatives = { ";
    for( std::vector< double >::const_iterator
         inputData( firstAndSecondDerivatives.begin() );
         inputData < firstAndSecondDerivatives.end();
         ++inputData )
    {
      if( inputData > firstAndSecondDerivatives.begin() )
      {
        std::cout << ", ";
      }
      std::cout << *inputData;
    }
    std::cout << " }, radialValue = " << radialValue << " ) called.";
    std::cout << std::endl;*/
    double const auxiliaryValue( auxiliaryAndFirstDerivative[ 0 ] );
    // This cheats if there has already been an overshoot, to try to avoid the
    // integration going to small step sizes to resolve the oscillations of the
    // field going off to infinity.
    if( auxiliaryValue < 0.0 )
    {
      firstAndSecondDerivatives[ 0 ] = 0.0;
      firstAndSecondDerivatives[ 1 ] = 0.0;
      return;
    }
    double const auxiliaryDerivative( auxiliaryAndFirstDerivative[ 1 ] );
    double firstDerivativeValue( 0.0 );
    double fieldDerivativeSquared( 0.0 );
    double fieldFirstDotSecondDerivatives( 0.0 );
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      firstDerivativeValue = firstDerivatives[ fieldIndex ]( auxiliaryValue );
      fieldDerivativeSquared
      += ( firstDerivativeValue * firstDerivativeValue );
      fieldFirstDotSecondDerivatives
      += ( firstDerivativeValue
           * secondDerivatives[ fieldIndex ]( auxiliaryValue ) );
    }
    firstAndSecondDerivatives[ 0 ] = auxiliaryDerivative;
    firstAndSecondDerivatives[ 1 ]
     = ( ( ( potentialSpline.FirstDerivative( auxiliaryValue )
             - ( fieldFirstDotSecondDerivatives
                 * auxiliaryDerivative * auxiliaryDerivative ) )
           / fieldDerivativeSquared )
         - ( ( dampingFactor * auxiliaryDerivative ) / radialValue ) );
  }

} /* namespace VevaciousPlusPlus */
