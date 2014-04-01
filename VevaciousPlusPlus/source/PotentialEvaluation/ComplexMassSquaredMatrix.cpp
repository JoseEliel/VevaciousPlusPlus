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
    MassesSquaredFromPolynomials( attributeMap ),
    numberOfRows( numberOfRows ),
    matrixElements( ( numberOfRows * numberOfRows ),
                    std::pair< PolynomialSum, PolynomialSum >( PolynomialSum(),
                                                           PolynomialSum() ) ),
    valuesMatrix( numberOfRows,
                  numberOfRows ),
    eigenvalueFinder( numberOfRows )
  {
    massesSquared.assign( numberOfRows,
                          NAN );
  }

  ComplexMassSquaredMatrix::ComplexMassSquaredMatrix(
                                 ComplexMassSquaredMatrix const& copySource ) :
    MassesSquaredFromPolynomials( copySource ),
    numberOfRows( copySource.numberOfRows ),
    matrixElements( copySource.matrixElements ),
    valuesMatrix( copySource.valuesMatrix ),
    eigenvalueFinder( copySource.eigenvalueFinder )
  {
    // This constructor is just an initialization list.
  }

  ComplexMassSquaredMatrix::ComplexMassSquaredMatrix() :
    MassesSquaredFromPolynomials(),
    numberOfRows( 0 ),
    matrixElements(),
    valuesMatrix(),
    eigenvalueFinder()
  {
    // This constructor is just an initialization list.
  }

  ComplexMassSquaredMatrix::~ComplexMassSquaredMatrix()
  {
    // This does nothing.
  }


  // This returns the eigenvalues of the square of the matrix.
  std::vector< double > const&
  ComplexMassSquaredMatrix::MassesSquared(
                              std::vector< double > const& fieldConfiguration )
  {
    unsigned int rowsTimesLength( 0 );
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
        /*
        valuesMatrix.coeffRef( columnIndex,
                               rowIndex ).real()
        = valuesMatrix.coeff( rowIndex,
                              columnIndex ).real();
        valuesMatrix.coeffRef( columnIndex,
                               rowIndex ).imag()
        = -(valuesMatrix.coeff( rowIndex,
                                columnIndex ).imag());
        */
      }
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).real()
      = matrixElements[ rowsTimesLength + rowIndex ].first(
                                                          fieldConfiguration );
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).imag() = 0.0;
      rowsTimesLength += numberOfRows;
    }

    eigenvalueFinder.compute( valuesMatrix,
                              Eigen::EigenvaluesOnly );
    for( unsigned int whichIndex( 0 );
         whichIndex < numberOfRows;
         ++whichIndex )
    {
      massesSquared[ whichIndex ]
      = eigenvalueFinder.eigenvalues()( whichIndex );
    }

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "making extra matrices!";
    std::cout << std::endl;
    Eigen::MatrixXcd zeroRowLessThanColumn( valuesMatrix );
    Eigen::MatrixXcd zeroColumnLessThanRow( valuesMatrix );
    for( unsigned int rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      for( unsigned int columnIndex( rowIndex + 1 );
           columnIndex < numberOfRows;
           ++columnIndex )
      {
        zeroRowLessThanColumn.coeffRef( rowIndex,
                                        columnIndex ).real() = 0.0;
        zeroRowLessThanColumn.coeffRef( rowIndex,
                                        columnIndex ).imag() = 0.0;
        zeroColumnLessThanRow.coeffRef( columnIndex,
                                        rowIndex ).real() = 0.0;
        zeroColumnLessThanRow.coeffRef( columnIndex,
                                        rowIndex ).imag() = 0.0;
      }
    }
    std::cout << std::endl;
    std::cout
    << "valuesMatrix = " << std::endl
    << valuesMatrix << std::endl
    << "zeroRowLessThanColumn = " << std::endl
    << zeroRowLessThanColumn << std::endl
    << "zeroColumnLessThanRow = " << std::endl
    << zeroColumnLessThanRow << std::endl;
    std::cout << std::endl;
    std::cout
    << "valuesMatrix = " << std::endl
    << valuesMatrix << std::endl
    << "valuesMatrix eigenvalues = " << std::endl
    << eigenvalueFinder.eigenvalues() << std::endl;
    eigenvalueFinder.compute( zeroRowLessThanColumn,
                              Eigen::EigenvaluesOnly );
    std::cout << std::endl;
    std::cout
    << "zeroRowLessThanColumn eigenvalues = " << std::endl
    << eigenvalueFinder.eigenvalues() << std::endl;
    eigenvalueFinder.compute( zeroColumnLessThanRow,
                              Eigen::EigenvaluesOnly );
    std::cout << std::endl;
    std::cout
    << "zeroColumnLessThanRow eigenvalues = " << std::endl
    << eigenvalueFinder.eigenvalues() << std::endl;
    std::cout << std::endl;*/

    return massesSquared;
  }

} /* namespace VevaciousPlusPlus */
