/*
 * SymmetricComplexMassMatrix.cpp
 *
 *  Created on: Mar 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  SymmetricComplexMassMatrix::SymmetricComplexMassMatrix(
                                               unsigned int const numberOfRows,
                   std::map< std::string, std::string > const& attributeMap ) :
    MassesSquaredFromPolynomials( attributeMap ),
    numberOfRows( numberOfRows ),
    matrixElements( ( numberOfRows * numberOfRows ),
                    std::pair< PolynomialSum, PolynomialSum >( PolynomialSum(),
                                                           PolynomialSum() ) ),
    valuesMatrix( numberOfRows,
                  numberOfRows ),
    valuesSquaredMatrix( numberOfRows,
                         numberOfRows ),
    eigenvalueFinder( numberOfRows )
  {
    massesSquared.assign( numberOfRows,
                          NAN );
  }

  SymmetricComplexMassMatrix::SymmetricComplexMassMatrix(
                               SymmetricComplexMassMatrix const& copySource ) :
    MassesSquaredFromPolynomials( copySource ),
    numberOfRows( copySource.numberOfRows ),
    matrixElements( copySource.matrixElements ),
    valuesMatrix( copySource.valuesMatrix ),
    valuesSquaredMatrix( copySource.valuesSquaredMatrix ),
    eigenvalueFinder( copySource.eigenvalueFinder )
  {
    // This constructor is just an initialization list.
  }

  SymmetricComplexMassMatrix::SymmetricComplexMassMatrix() :
    MassesSquaredFromPolynomials(),
    numberOfRows( 0 ),
    matrixElements(),
    valuesMatrix(),
    valuesSquaredMatrix(),
    eigenvalueFinder()
  {
    // This constructor is just an initialization list.
  }

  SymmetricComplexMassMatrix::~SymmetricComplexMassMatrix()
  {
    // This does nothing.
  }


  // This returns the eigenvalues of the square of the matrix.
  std::vector< double > const&
  SymmetricComplexMassMatrix::MassesSquared(
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
    for( unsigned int rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      for( unsigned int columnIndex( 0 );
           columnIndex <= rowIndex;
           ++columnIndex )
      {
        valuesSquaredMatrix.coeffRef( rowIndex,
                                      columnIndex ).real() = 0.0;
        valuesSquaredMatrix.coeffRef( rowIndex,
                                      columnIndex ).imag() = 0.0;
        for( unsigned int sumIndex( 0 );
             sumIndex < numberOfRows;
             ++sumIndex )
        {
          valuesSquaredMatrix.coeffRef( rowIndex,
                                        columnIndex ).real()
          += ( ( valuesMatrix.coeff( sumIndex,
                                     rowIndex ).real()
                 * valuesMatrix.coeff( sumIndex,
                                       columnIndex ).real() )
               + ( valuesMatrix.coeff( sumIndex,
                                       rowIndex ).imag()
                   * valuesMatrix.coeff( sumIndex,
                                         columnIndex ).imag() ) );
          valuesSquaredMatrix.coeffRef( rowIndex,
                                        columnIndex ).imag()
          += ( ( valuesMatrix.coeff( sumIndex,
                                     rowIndex ).real()
                 * valuesMatrix.coeff( sumIndex,
                                       columnIndex ).imag() )
               - ( valuesMatrix.coeff( sumIndex,
                                       rowIndex ).imag()
                   * valuesMatrix.coeff( sumIndex,
                                         columnIndex ).real() ) );
          // The Eigen routines don't bother looking at elements of
          // valuesSquaredMatrix where columnIndex > rowIndex, so we don't even
          // bother filling them with the conjugates of the transpose.
        }
      }
    }

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "valuesMatrix = " << std::endl
    << valuesMatrix << std::endl
    << "valuesMatrix.adjoint() = " << std::endl
    << valuesMatrix.adjoint() << std::endl
    << "valuesSquaredMatrix = " << std::endl
    << valuesSquaredMatrix << std::endl
    << "valuesMatrix.adjoint() * valuesMatrix = " << std::endl;
    Eigen::MatrixXcd productMatrix( valuesMatrix.adjoint() * valuesMatrix );
    std::cout << productMatrix << std::endl;
    std::cout << std::endl;*/

    eigenvalueFinder.compute( valuesSquaredMatrix,
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
    << "massesSquared = {" << std::endl;
    for( unsigned int whichIndex( 0 );
         whichIndex < numberOfRows;
         ++whichIndex )
    {
      std::cout << " " << massesSquared[ whichIndex ];
    }
    std::cout
    << " }" << std::endl << "eigenvalueFinder.eigenvalues() = " << std::endl
    << eigenvalueFinder.eigenvalues() << std::endl;
    eigenvalueFinder.compute( productMatrix,
                              Eigen::EigenvaluesOnly );
    std::cout
    << "for productMatrix, eigenvalueFinder.eigenvalues() = " << std::endl
    << eigenvalueFinder.eigenvalues() << std::endl;
    std::cout << std::endl;*/

    return massesSquared;
  }

} /* namespace VevaciousPlusPlus */
