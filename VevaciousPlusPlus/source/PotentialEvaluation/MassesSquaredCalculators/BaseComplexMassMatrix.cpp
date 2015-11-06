/*
 * BaseComplexMassMatrix.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/MassesSquaredCalculators/BaseComplexMassMatrix.hpp"

namespace VevaciousPlusPlus
{

  BaseComplexMassMatrix::BaseComplexMassMatrix( size_t const numberOfRows,
                   std::map< std::string, std::string > const& attributeMap ) :
    MassesSquaredFromMatrix< std::complex< double > >( numberOfRows,
                                                       attributeMap ),
    matrixElements( ( numberOfRows * numberOfRows ),
         ComplexParametersAndFieldsProductSum( ParametersAndFieldsProductSum(),
                                            ParametersAndFieldsProductSum() ) )
  {
    // This constructor is just an initialization list.
  }

  BaseComplexMassMatrix::BaseComplexMassMatrix(
                                    BaseComplexMassMatrix const& copySource ) :
    MassesSquaredFromMatrix< std::complex< double > >( copySource ),
    matrixElements( copySource.matrixElements )
  {
    // This constructor is just an initialization list.
  }

  BaseComplexMassMatrix::BaseComplexMassMatrix() :
    MassesSquaredFromMatrix< std::complex< double > >(),
    matrixElements()
  {
    // This constructor is just an initialization list.
  }

  BaseComplexMassMatrix::~BaseComplexMassMatrix()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
