/*
 * RealMassesSquaredMatrix.cpp
 *
 *  Created on: Oct 09, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/MassesSquaredCalculators/RealMassesSquaredMatrix.hpp"

namespace VevaciousPlusPlus
{

  RealMassesSquaredMatrix::RealMassesSquaredMatrix( size_t const numberOfRows,
                   std::map< std::string, std::string > const& attributeMap ) :
    MassesSquaredFromMatrix< double >( numberOfRows,
                                       attributeMap ),
    matrixElements( ( numberOfRows * numberOfRows ),
                    ParametersAndFieldsProductSum() )
  {
    // This constructor is just an initialization list.
  }

  RealMassesSquaredMatrix::RealMassesSquaredMatrix(
                                  RealMassesSquaredMatrix const& copySource ) :
    MassesSquaredFromMatrix< double >( copySource ),
    matrixElements( copySource.matrixElements )
  {
    // This constructor is just an initialization list.
  }


  RealMassesSquaredMatrix::RealMassesSquaredMatrix() :
    MassesSquaredFromMatrix< double >(),
    matrixElements()
  {
    // This constructor is just an initialization list.
  }

  RealMassesSquaredMatrix::~RealMassesSquaredMatrix()
  {
    // This does nothing.
  }


  // This returns a matrix of the values of the elements for a field
  // configuration given by fieldConfiguration, using the values for the
  // Lagrangian parameters found in parameterValues.
  Eigen::MatrixXd RealMassesSquaredMatrix::CurrentValues(
                                  std::vector< double > const& parameterValues,
                        std::vector< double > const& fieldConfiguration ) const
  {
    size_t rowsTimesLength( 0 );
    Eigen::MatrixXd valuesMatrix( numberOfRows,
                                  numberOfRows );
    for( size_t rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex )
      = matrixElements[ rowsTimesLength + rowIndex ]( parameterValues,
                                                      fieldConfiguration );
      for( size_t columnIndex( rowIndex + 1 );
           columnIndex < numberOfRows;
           ++columnIndex )
      {
        valuesMatrix.coeffRef( rowIndex,
                               columnIndex )
        = matrixElements[ rowsTimesLength + columnIndex ]( parameterValues,
                                                          fieldConfiguration );
        valuesMatrix.coeffRef( columnIndex,
                               rowIndex ) = valuesMatrix.coeff( rowIndex,
                                                                columnIndex );
      }
      rowsTimesLength += numberOfRows;
    }
    return valuesMatrix;
  }

  // This returns a matrix of the values of the elements for a field
  // configuration given by fieldConfiguration, using the values for the
  // Lagrangian parameters from the last call of UpdateForFixedScale.
  Eigen::MatrixXd RealMassesSquaredMatrix::CurrentValues(
                        std::vector< double > const& fieldConfiguration ) const
  {
    size_t rowsTimesLength( 0 );
    Eigen::MatrixXd valuesMatrix( numberOfRows,
                                  numberOfRows );
    for( size_t rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex )
      = matrixElements[ rowsTimesLength + rowIndex ]( fieldConfiguration );
      for( size_t columnIndex( rowIndex + 1 );
           columnIndex < numberOfRows;
           ++columnIndex )
      {
        valuesMatrix.coeffRef( rowIndex,
                               columnIndex )
        = matrixElements[ rowsTimesLength + columnIndex ](
                                                          fieldConfiguration );
        valuesMatrix.coeffRef( columnIndex,
                               rowIndex ) = valuesMatrix.coeff( rowIndex,
                                                                columnIndex );
      }
      rowsTimesLength += numberOfRows;
    }
    return valuesMatrix;
  }

} /* namespace VevaciousPlusPlus */
