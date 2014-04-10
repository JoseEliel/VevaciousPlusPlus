/*
 * ComplexMassSquaredMatrix.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  ComplexMassSquaredMatrix::ComplexMassSquaredMatrix(
                                               unsigned int const numberOfRows,
                   std::map< std::string, std::string > const& attributeMap ) :
    MassesSquaredFromMatrix< std::complex< double > >( numberOfRows,
                                                       attributeMap ),
    matrixElements( ( numberOfRows * numberOfRows ),
                    std::pair< PolynomialSum, PolynomialSum >( PolynomialSum(),
                                                           PolynomialSum() ) )
  {
    // This constructor is just an initialization list.
  }

  ComplexMassSquaredMatrix::ComplexMassSquaredMatrix(
                                 ComplexMassSquaredMatrix const& copySource ) :
    MassesSquaredFromMatrix< std::complex< double > >( copySource ),
    matrixElements( copySource.matrixElements )
  {
    // This constructor is just an initialization list.
  }

  ComplexMassSquaredMatrix::ComplexMassSquaredMatrix() :
    MassesSquaredFromMatrix< std::complex< double > >(),
    matrixElements()
  {
    // This constructor is just an initialization list.
  }

  ComplexMassSquaredMatrix::~ComplexMassSquaredMatrix()
  {
    // This does nothing.
  }


  // This returns a matrix of the values of the elements for a field
  // configuration given by fieldConfiguration, with all functionoids
  // evaluated at the last scale which was used to update them.
  Eigen::MatrixXcd ComplexMassSquaredMatrix::CurrentValues(
                        std::vector< double > const& fieldConfiguration ) const
  {
    unsigned int rowsTimesLength( 0 );
    Eigen::MatrixXcd valuesMatrix( numberOfRows,
                                   numberOfRows );
    for( unsigned int rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      for( unsigned int columnIndex( 0 );
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
        // The Eigen routines don't bother looking at elements of valuesMatrix
        // where columnIndex > rowIndex, so we don't even bother filling them
        // with the conjugates of the transpose.
      }
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).real()
      = matrixElements[ rowsTimesLength + rowIndex ].first(
                                                          fieldConfiguration );
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).imag() = 0.0;
      rowsTimesLength += numberOfRows;
    }
    return valuesMatrix;
  }

  // This returns a matrix of the values of the elements for a field
  // configuration given by fieldConfiguration, with all functionoids
  // evaluated at the natural exponent of logarithmOfScale.
  Eigen::MatrixXcd ComplexMassSquaredMatrix::CurrentValues(
                               std::vector< double > const& fieldConfiguration,
                                          double const logarithmOfScale ) const
  {
    unsigned int rowsTimesLength( 0 );
    Eigen::MatrixXcd valuesMatrix( numberOfRows,
                                   numberOfRows );
    for( unsigned int rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      for( unsigned int columnIndex( 0 );
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
        // The Eigen routines don't bother looking at elements of valuesMatrix
        // where columnIndex > rowIndex, so we don't even bother filling them
        // with the conjugates of the transpose.
      }
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).real()
      = matrixElements[ rowsTimesLength + rowIndex ].first( fieldConfiguration,
                                                            logarithmOfScale );
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).imag() = 0.0;
      rowsTimesLength += numberOfRows;
    }
    return valuesMatrix;
  }

} /* namespace VevaciousPlusPlus */
