/*
 * RealMassesSquaredMatrix.cpp
 *
 *  Created on: Mar 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  RealMassesSquaredMatrix::RealMassesSquaredMatrix(
                                               unsigned int const numberOfRows,
                   std::map< std::string, std::string > const& attributeMap ) :
    MassesSquaredFromPolynomials( attributeMap ),
    numberOfRows( numberOfRows ),
    matrixElements( ( numberOfRows * numberOfRows ),
                    PolynomialSum() ),
    valuesMatrix( numberOfRows,
                  numberOfRows ),
    eigenvalueFinder( numberOfRows )
  {
    massesSquared.assign( numberOfRows,
                          NAN );
  }

  RealMassesSquaredMatrix::RealMassesSquaredMatrix(
                                  RealMassesSquaredMatrix const& copySource ) :
    MassesSquaredFromPolynomials( copySource ),
    numberOfRows( copySource.numberOfRows ),
    matrixElements( copySource.matrixElements ),
    valuesMatrix( copySource.valuesMatrix ),
    eigenvalueFinder( copySource.eigenvalueFinder )
  {
    // This constructor is just an initialization list.
  }


  RealMassesSquaredMatrix::RealMassesSquaredMatrix() :
    MassesSquaredFromPolynomials(),
    numberOfRows( 0 ),
    matrixElements(),
    valuesMatrix(),
    eigenvalueFinder()
  {
    // This constructor is just an initialization list.
  }

  RealMassesSquaredMatrix::~RealMassesSquaredMatrix()
  {
    // This does nothing.
  }


  // This returns the eigenvalues of the matrix.
  std::vector< double > const&
  RealMassesSquaredMatrix::MassesSquared(
                              std::vector< double > const& fieldConfiguration )
  {
    unsigned int rowsTimesLength( 0 );
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

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "valuesMatrix = " << std::endl
    << valuesMatrix << std::endl;
    std::cout << std::endl;/**/

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
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "massesSquared = {" << std::endl;
    for( unsigned int whichIndex( 0 );
         whichIndex < numberOfRows;
         ++whichIndex )
    {
      std::cout << " " << massesSquared[ whichIndex ];
    }
    std::cout << " }" << std::endl;/**/

    return massesSquared;
  }

} /* namespace VevaciousPlusPlus */
