/*
 * RealMassesSquaredMatrix.cpp
 *
 *  Created on: Mar 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  RealMassesSquaredMatrix::RealMassesSquaredMatrix(
                                               unsigned int const numberOfRows,
                   std::map< std::string, std::string > const& attributeMap ) :
    MassesSquaredFromMatrix< double >( numberOfRows,
                                       attributeMap ),
    matrixElements( ( numberOfRows * numberOfRows ),
                    PolynomialSum() )
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
  // configuration given by fieldConfiguration, with all functionoids
  // evaluated at the last scale which was used to update them.
  Eigen::MatrixXd RealMassesSquaredMatrix::CurrentValues(
                        std::vector< double > const& fieldConfiguration ) const
  {
    unsigned int rowsTimesLength( 0 );
    Eigen::MatrixXd valuesMatrix( numberOfRows,
                                  numberOfRows );
    for( unsigned int rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex )
      = matrixElements[ rowsTimesLength + rowIndex ]( fieldConfiguration );
      for( unsigned int columnIndex( rowIndex + 1 );
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

  // This returns a matrix of the values of the elements for a field
  // configuration given by fieldConfiguration, with all functionoids
  // evaluated at the natural exponent of logarithmOfScale.
  Eigen::MatrixXd RealMassesSquaredMatrix::CurrentValues(
                               std::vector< double > const& fieldConfiguration,
                                          double const logarithmOfScale ) const
  {
    unsigned int rowsTimesLength( 0 );
    Eigen::MatrixXd valuesMatrix( numberOfRows,
                                  numberOfRows );
    for( unsigned int rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex )
      = matrixElements[ rowsTimesLength + rowIndex ]( fieldConfiguration,
                                                      logarithmOfScale );
      for( unsigned int columnIndex( rowIndex + 1 );
           columnIndex < numberOfRows;
           ++columnIndex )
      {
        valuesMatrix.coeffRef( rowIndex,
                               columnIndex )
        = matrixElements[ rowsTimesLength + columnIndex ]( fieldConfiguration,
                                                           logarithmOfScale );
        valuesMatrix.coeffRef( columnIndex,
                               rowIndex ) = valuesMatrix.coeff( rowIndex,
                                                                columnIndex );
      }
      rowsTimesLength += numberOfRows;
    }
    return valuesMatrix;
  }

} /* namespace VevaciousPlusPlus */
