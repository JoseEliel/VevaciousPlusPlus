/*
 * SymmetricComplexMassMatrix.cpp
 *
 *  Created on: Mar 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/MassesSquaredCalculators/OldSymmetricComplexMassMatrix.hpp"

namespace VevaciousPlusPlus
{

  OldSymmetricComplexMassMatrix::OldSymmetricComplexMassMatrix(
                                                     size_t const numberOfRows,
                   std::map< std::string, std::string > const& attributeMap ) :
    OldMassesSquaredFromMatrix< std::complex< double > >( numberOfRows,
                                                          attributeMap ),
    matrixElements( ( numberOfRows * numberOfRows ),
                    std::pair< PolynomialSum, PolynomialSum >( PolynomialSum(),
                                                            PolynomialSum() ) )
  {
    // This constructor is just an initialization list.
  }

  OldSymmetricComplexMassMatrix::OldSymmetricComplexMassMatrix(
                               OldSymmetricComplexMassMatrix const& copySource ) :
    OldMassesSquaredFromMatrix( copySource ),
    matrixElements( copySource.matrixElements )
  {
    // This constructor is just an initialization list.
  }

  OldSymmetricComplexMassMatrix::OldSymmetricComplexMassMatrix() :
    OldMassesSquaredFromMatrix(),
    matrixElements()
  {
    // This constructor is just an initialization list.
  }

  OldSymmetricComplexMassMatrix::~OldSymmetricComplexMassMatrix()
  {
    // This does nothing.
  }


  // This returns a matrix of the values of the elements for a field
  // configuration given by fieldConfiguration, with all functionoids
  // evaluated at the last scale which was used to update them.
  Eigen::MatrixXcd OldSymmetricComplexMassMatrix::MatrixToSquare(
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

  // This returns a matrix of the values of the elements for a field
  // configuration given by fieldConfiguration, with all functionoids
  // evaluated at the natural exponent of logarithmOfScale.
  Eigen::MatrixXcd OldSymmetricComplexMassMatrix::MatrixToSquare(
                               std::vector< double > const& fieldConfiguration,
                                          double const logarithmOfScale ) const
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
                                                            fieldConfiguration,
                                                            logarithmOfScale );
        valuesMatrix.coeffRef( rowIndex,
                               columnIndex ).imag()
        = matrixElements[ rowsTimesLength + columnIndex ].second(
                                                            fieldConfiguration,
                                                            logarithmOfScale );
        // We use the fact that the matrix is symmetric.
        valuesMatrix.coeffRef( columnIndex,
                               rowIndex )
        = valuesMatrix.coeff( rowIndex,
                              columnIndex );
      }
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).real()
      = matrixElements[ rowsTimesLength + rowIndex ].first( fieldConfiguration,
                                                            logarithmOfScale );
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).imag()
      = matrixElements[ rowsTimesLength + rowIndex ].second(
                                                            fieldConfiguration,
                                                            logarithmOfScale );
      rowsTimesLength += numberOfRows;
    }
    return valuesMatrix;
  }

  // This returns a matrix that is the lower-triangular part (only column
  // index <= row index) of the square of matrixToSquare.
  Eigen::MatrixXcd OldSymmetricComplexMassMatrix::LowerTriangleOfSquareMatrix(
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
