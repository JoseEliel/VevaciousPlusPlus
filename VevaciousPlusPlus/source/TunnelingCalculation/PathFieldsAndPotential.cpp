/*
 * PathFieldsAndPotential.cpp
 *
 *  Created on: Jun 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/TunnelingCalculation.hpp"

namespace VevaciousPlusPlus
{

  PathFieldsAndPotential::PathFieldsAndPotential(
                                       Eigen::MatrixXd const& pathCoefficients,
                         std::vector< double > const& falseVacuumConfiguration,
                                                double const falseVacuumDepth,
                                                 double const trueVacuumDepth,
                                              double const givenTemperature ) :
    potentialApproximation(),
    numberOfFields( pathCoefficients.cols() ),
    fieldPath( numberOfFields,
               SimplePolynomial( pathCoefficients.rows() + 1 ) ),
    pathTangent( fieldPath ),
    fieldConfiguration( fieldPath.size(),
                        0.0 ),
    falseVacuumDepth( falseVacuumDepth ),
    trueVacuumDepth( trueVacuumDepth ),
    nonZeroTemperature( givenTemperature > 0.0 ),
    givenTemperature( givenTemperature )
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "PathFieldsAndPotential::PathFieldsAndPotential( ... ) called.";
    std::cout << std::endl;
    std::cout << "pathCoefficients =" << std::endl << pathCoefficients;
    std::cout << std::endl;/**/
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      std::vector< double >&
      coefficientVector( fieldPath[ fieldIndex ].CoefficientVector() );
      coefficientVector[ 0 ] = falseVacuumConfiguration[ fieldIndex ];
      for( size_t coefficientIndex( 0 );
           coefficientIndex < pathCoefficients.rows();
           ++coefficientIndex )
      {
        coefficientVector[ coefficientIndex + 1 ]
        = pathCoefficients( coefficientIndex,
                            fieldIndex );
      }
      pathTangent[ fieldIndex ].BecomeFirstDerivativeOf(
                                                     fieldPath[ fieldIndex ] );
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "PathFieldsAndPotential::PathFieldsAndPotential(...) finishing;"
    << " fieldPath now is"
    << std::endl;
    for( std::vector< SimplePolynomial >::const_iterator
         fieldAsPolynomial( fieldPath.begin() );
         fieldAsPolynomial < fieldPath.end();
         ++fieldAsPolynomial )
    {
      std::cout << fieldAsPolynomial->AsDebuggingString() << std::endl;
    }
    std::cout << "pathTangent now is"
    << std::endl;
    for( std::vector< SimplePolynomial >::const_iterator
         fieldAsPolynomial( pathTangent.begin() );
         fieldAsPolynomial < pathTangent.end();
         ++fieldAsPolynomial )
    {
      std::cout << fieldAsPolynomial->AsDebuggingString() << std::endl;
    }
    std::cout << "(2 sets of " << pathTangent.size()
    << " polynomials)";
    std::cout << std::endl;/**/
  }

  PathFieldsAndPotential::PathFieldsAndPotential(
                                   PathFieldsAndPotential const& copySource ) :
    potentialApproximation( copySource.potentialApproximation ),
    numberOfFields( copySource.numberOfFields ),
    fieldPath( copySource.fieldPath ),
    pathTangent( copySource.pathTangent ),
    fieldConfiguration( copySource.fieldConfiguration ),
    falseVacuumDepth( copySource.falseVacuumDepth ),
    trueVacuumDepth( copySource.trueVacuumDepth ),
    nonZeroTemperature( copySource.nonZeroTemperature ),
    givenTemperature( copySource.givenTemperature )
  {
    // This constructor is just an initialization list.
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "PathFieldsAndPotential copy constructor called.";
    std::cout << std::endl;/**/
  }

  PathFieldsAndPotential::~PathFieldsAndPotential()
  {
    // This does nothing.
  }


  // This is for debugging.
  std::string PathFieldsAndPotential::AsDebuggingString() const
  {
    std::stringstream returnStream;
    returnStream << "potentialApproximation = "
    << potentialApproximation.AsDebuggingString() << std::endl
    << "numberOfFields = " << numberOfFields << std::endl << "fieldPath = { ";
    for( std::vector< SimplePolynomial >::const_iterator
         fieldAsPolynomial( fieldPath.begin() );
         fieldAsPolynomial < fieldPath.end();
         ++fieldAsPolynomial )
    {
      if( fieldAsPolynomial != fieldPath.begin() )
      {
        returnStream << "," << std::endl;
      }
      returnStream << fieldAsPolynomial->AsDebuggingString();
    }
    returnStream << " }" << std::endl << "pathTangent = {";
    for( std::vector< SimplePolynomial >::const_iterator
         fieldAsPolynomial( pathTangent.begin() );
         fieldAsPolynomial < pathTangent.end();
         ++fieldAsPolynomial )
    {
      if( fieldAsPolynomial != pathTangent.begin() )
      {
        returnStream << "," << std::endl;
      }
      returnStream << fieldAsPolynomial->AsDebuggingString();
    }
    returnStream << " }" << std::endl
    << "falseVacuumDepth = " << falseVacuumDepth << ", trueVacuumDepth = "
    << trueVacuumDepth << ", nonZeroTemperature = " << nonZeroTemperature
    << ", givenTemperature = " << givenTemperature;
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
