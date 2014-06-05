/*
 * PathFieldsAndPotential.cpp
 *
 *  Created on: Jun 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/TunnelingCalculation/PathFieldsAndPotential.hpp"

namespace VevaciousPlusPlus
{

  PathFieldsAndPotential::PathFieldsAndPotential() :
    fieldConfiguration( fieldPath.size() )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "STILL NEEDS WORK!";
    std::cout << std::endl;/**/
    fieldsAsPolynomials.resize( numberOfFields,
                            SimplePolynomial( numberOfVaryingPathNodes + 2 ) );
    fieldDerivativesAsPolynomials = fieldsAsPolynomials;
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      std::vector< double >& coefficientVector(
                       fieldsAsPolynomials[ fieldIndex ].CoefficientVector() );
      coefficientVector[ 0 ] = falseVacuumConfiguration[ fieldIndex ];
      for( unsigned int coefficientIndex( 0 );
           coefficientIndex <= numberOfVaryingPathNodes;
           ++coefficientIndex )
      {
        coefficientVector[ coefficientIndex + 1 ]
        = pathCoefficients( coefficientIndex,
                            fieldIndex );
      }
      fieldDerivativesAsPolynomials[ fieldIndex ]
      = fieldsAsPolynomials[ fieldIndex ].FirstDerivative();
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "PathFromNodes::operator() finishing; fieldsAsPolynomials now is"
    << std::endl;
    for( std::vector< SimplePolynomial >::const_iterator
         fieldAsPolynomial( fieldsAsPolynomials.begin() );
         fieldAsPolynomial < fieldsAsPolynomials.end();
         ++fieldAsPolynomial )
    {
      std::cout << fieldAsPolynomial->AsDebuggingString() << std::endl;
    }
    std::cout << "fieldDerivativesAsPolynomials now is"
    << std::endl;
    for( std::vector< SimplePolynomial >::const_iterator
         fieldAsPolynomial( fieldDerivativesAsPolynomials.begin() );
         fieldAsPolynomial < fieldDerivativesAsPolynomials.end();
         ++fieldAsPolynomial )
    {
      std::cout << fieldAsPolynomial->AsDebuggingString() << std::endl;
    }
    std::cout << "(2 sets of " << fieldDerivativesAsPolynomials.size()
    << " polynomials)";
    std::cout << std::endl;/**/
  }

  PathFieldsAndPotential::~PathFieldsAndPotential()
  {
    // TODO Auto-generated destructor stub
  }

} /* namespace VevaciousPlusPlus */
