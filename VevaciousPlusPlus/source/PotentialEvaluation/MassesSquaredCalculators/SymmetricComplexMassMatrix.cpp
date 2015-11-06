/*
 * SymmetricComplexMassMatrix.cpp
 *
 *  Created on: Oct 09, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/MassesSquaredCalculators/SymmetricComplexMassMatrix.hpp"

namespace VevaciousPlusPlus
{

  SymmetricComplexMassMatrix::SymmetricComplexMassMatrix(
                                                     size_t const numberOfRows,
                   std::map< std::string, std::string > const& attributeMap ) :
    BaseComplexMassMatrix( numberOfRows,
                           attributeMap )
  {
    // This constructor is just an initialization list.
  }

  SymmetricComplexMassMatrix::SymmetricComplexMassMatrix(
                               SymmetricComplexMassMatrix const& copySource ) :
    BaseComplexMassMatrix( copySource )
  {
    // This constructor is just an initialization list.
  }

  SymmetricComplexMassMatrix::SymmetricComplexMassMatrix() :
    BaseComplexMassMatrix()
  {
    // This constructor is just an initialization list.
  }

  SymmetricComplexMassMatrix::~SymmetricComplexMassMatrix()
  {
    // This does nothing.
  }


  // This returns a matrix of the values of the elements for a field
  // configuration given by fieldConfiguration, using the values for the
  // Lagrangian parameters found in parameterValues.
  Eigen::MatrixXcd SymmetricComplexMassMatrix::MatrixToSquare(
                                  std::vector< double > const& parameterValues,
                        std::vector< double > const& fieldConfiguration ) const
  {
    size_t rowsTimesLength( 0 );
    Eigen::MatrixXcd valuesMatrix( numberOfRows,
                                   numberOfRows );
    for( size_t rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      for( size_t columnIndex( 0 );
           columnIndex < rowIndex;
           ++columnIndex )
      {
        valuesMatrix.coeffRef( rowIndex,
                               columnIndex ).real()
        = matrixElements[ rowsTimesLength + columnIndex ].first(
                                                               parameterValues,
                                                          fieldConfiguration );
        valuesMatrix.coeffRef( rowIndex,
                               columnIndex ).imag()
        = matrixElements[ rowsTimesLength + columnIndex ].second(
                                                               parameterValues,
                                                          fieldConfiguration );
        // We use the fact that the matrix is symmetric.
        valuesMatrix.coeffRef( columnIndex,
                               rowIndex )
        = valuesMatrix.coeff( rowIndex,
                              columnIndex );
      }
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).real()
      = matrixElements[ rowsTimesLength + rowIndex ].first( parameterValues,
                                                          fieldConfiguration );
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).imag()
      = matrixElements[ rowsTimesLength + rowIndex ].second( parameterValues,
                                                          fieldConfiguration );
      rowsTimesLength += numberOfRows;
    }
    return valuesMatrix;
  }

  // This returns a matrix of the values of the elements for a field
  // configuration given by fieldConfiguration, using the values for the
  // Lagrangian parameters from the last call of UpdateForFixedScale.
  Eigen::MatrixXcd SymmetricComplexMassMatrix::MatrixToSquare(
                        std::vector< double > const& fieldConfiguration ) const
  {
    size_t rowsTimesLength( 0 );
    Eigen::MatrixXcd valuesMatrix( numberOfRows,
                                   numberOfRows );
    for( size_t rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      for( size_t columnIndex( 0 );
           columnIndex < rowIndex;
           ++columnIndex )
      {
        valuesMatrix.coeffRef( rowIndex,
                               columnIndex ).real()
        = matrixElements[ rowsTimesLength + columnIndex ].first(
                                                          fieldConfiguration );
        valuesMatrix.coeffRef( rowIndex,
                               columnIndex ).imag()
        = matrixElements[ rowsTimesLength + columnIndex ].second(
                                                          fieldConfiguration );
        // We use the fact that the matrix is symmetric.
        valuesMatrix.coeffRef( columnIndex,
                               rowIndex )
        = valuesMatrix.coeff( rowIndex,
                              columnIndex );
      }
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).real()
      = matrixElements[ rowsTimesLength + rowIndex ].first(
                                                          fieldConfiguration );
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).imag()
      = matrixElements[ rowsTimesLength + rowIndex ].second(
                                                          fieldConfiguration );
      rowsTimesLength += numberOfRows;
    }
    return valuesMatrix;
  }

  // This returns a matrix that is the lower-triangular part (only column
  // index <= row index) of the square of matrixToSquare.
  Eigen::MatrixXcd SymmetricComplexMassMatrix::LowerTriangleOfSquareMatrix(
                                 Eigen::MatrixXcd const& matrixToSquare ) const
  {
    Eigen::MatrixXcd valuesSquaredMatrix( numberOfRows,
                                          numberOfRows );
    for( size_t rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      for( size_t columnIndex( 0 );
           columnIndex <= rowIndex;
           ++columnIndex )
      {
        valuesSquaredMatrix.coeffRef( rowIndex,
                                      columnIndex ).real() = 0.0;
        valuesSquaredMatrix.coeffRef( rowIndex,
                                      columnIndex ).imag() = 0.0;
        for( size_t sumIndex( 0 );
             sumIndex < numberOfRows;
             ++sumIndex )
        {
          valuesSquaredMatrix.coeffRef( rowIndex,
                                        columnIndex ).real()
          += ( ( matrixToSquare.coeff( sumIndex,
                                     rowIndex ).real()
                 * matrixToSquare.coeff( sumIndex,
                                       columnIndex ).real() )
               + ( matrixToSquare.coeff( sumIndex,
                                       rowIndex ).imag()
                   * matrixToSquare.coeff( sumIndex,
                                         columnIndex ).imag() ) );
          valuesSquaredMatrix.coeffRef( rowIndex,
                                        columnIndex ).imag()
          += ( ( matrixToSquare.coeff( sumIndex,
                                     rowIndex ).real()
                 * matrixToSquare.coeff( sumIndex,
                                       columnIndex ).imag() )
               - ( matrixToSquare.coeff( sumIndex,
                                       rowIndex ).imag()
                   * matrixToSquare.coeff( sumIndex,
                                         columnIndex ).real() ) );
          // The Eigen routines don't bother looking at elements of
          // valuesSquaredMatrix where columnIndex > rowIndex, so we don't even
          // bother filling them with the conjugates of the transpose.
        }
      }
    }
    return valuesSquaredMatrix;
  }

} /* namespace VevaciousPlusPlus */
