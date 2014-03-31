/*
 * ComplexMassMatrix.cpp
 *
 *  Created on: Mar 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  ComplexMassMatrix::ComplexMassMatrix( unsigned int const numberOfRows,
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

  ComplexMassMatrix::ComplexMassMatrix( ComplexMassMatrix const& copySource ) :
    MassesSquaredFromPolynomials( copySource ),
    numberOfRows( copySource.numberOfRows ),
    matrixElements( copySource.matrixElements ),
    valuesMatrix( copySource.valuesMatrix ),
    valuesSquaredMatrix( copySource.valuesSquaredMatrix ),
    eigenvalueFinder( copySource.eigenvalueFinder )
  {
    // This constructor is just an initialization list.
  }

  ComplexMassMatrix::ComplexMassMatrix() :
    MassesSquaredFromPolynomials(),
    numberOfRows( 0 ),
    matrixElements(),
    valuesMatrix(),
    valuesSquaredMatrix(),
    eigenvalueFinder()
  {
    // This constructor is just an initialization list.
  }

  ComplexMassMatrix::~ComplexMassMatrix()
  {
    // This does nothing.
  }


  // This returns the eigenvalues of the square of the matrix.
  std::vector< double > const&
  ComplexMassMatrix::MassesSquared(
                              std::vector< double > const& fieldConfiguration )
  {
    unsigned int rowsTimesLength( 0 );
    for( unsigned int rowIndex( 0 );
         rowIndex < numberOfRows;
         ++rowIndex )
    {
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).real()
      = matrixElements[ rowsTimesLength + rowIndex ].first(
                                                          fieldConfiguration );
      valuesMatrix.coeffRef( rowIndex,
                             rowIndex ).imag() = 0.0;
      valuesSquaredMatrix.coeffRef( rowIndex,
                                 rowIndex ) = valuesMatrix.coeff( rowIndex,
                                                                  rowIndex );
      for( unsigned int columnIndex( rowIndex + 1 );
           columnIndex < numberOfRows;
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
        valuesMatrix.coeffRef( columnIndex,
                               rowIndex ).real()
        = valuesMatrix.coeff( rowIndex,
                              columnIndex ).real();
        valuesMatrix.coeffRef( columnIndex,
                               rowIndex ).imag()
        = -(valuesMatrix.coeff( rowIndex,
                                columnIndex ).imag());
        valuesSquaredMatrix.coeffRef( columnIndex,
                                      rowIndex )
        = valuesMatrix.coeff( columnIndex,
                              rowIndex );
        valuesSquaredMatrix.coeffRef( rowIndex,
                                      columnIndex )
        = valuesMatrix.coeff( rowIndex,
                              columnIndex );
      }
      rowsTimesLength += numberOfRows;
    }
    valuesSquaredMatrix *= valuesMatrix;

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "valuesMatrix = " << std::endl
    << valuesMatrix << std::endl
    << "valuesSquaredMatrix = " << std::endl
    << valuesSquaredMatrix;
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
    std::cout << " }" << std::endl;*/

    return massesSquared;
  }

} /* namespace VevaciousPlusPlus */
