/*
 * RealMassesSquaredMatrix.hpp
 *
 *  Created on: Mar 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef REALMASSESSQUAREDMATRIX_HPP_
#define REALMASSESSQUAREDMATRIX_HPP_

#include "../StandardIncludes.hpp"
#include "MassesSquaredFromPolynomials.hpp"
#include "Eigen/Dense"

namespace VevaciousPlusPlus
{

  class RealMassesSquaredMatrix : public MassesSquaredFromPolynomials
  {
  public:
    RealMassesSquaredMatrix( unsigned int const numberOfRows,
                    std::map< std::string, std::string > const& attributeMap );
    RealMassesSquaredMatrix( RealMassesSquaredMatrix const& copySource );
    RealMassesSquaredMatrix();
    virtual
    ~RealMassesSquaredMatrix();


    // This returns the eigenvalues of the matrix.
    virtual std::vector< double > const&
    MassesSquared( std::vector< double > const& fieldConfiguration );

    // This allows access to the polynomial sum for a given index.
    PolynomialSum& ElementAt( unsigned int const elementIndex )
    { return matrixElements[ elementIndex ]; }

    // This is mainly for debugging:
    std::string AsString() const;

  protected:
    unsigned int numberOfRows;
    std::vector< PolynomialSum > matrixElements;
    Eigen::MatrixXd valuesMatrix;
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > eigenvalueFinder;
  };




  // This is mainly for debugging:
  inline std::string RealMassesSquaredMatrix::AsString() const
  {
    std::stringstream returnStream;
    returnStream
    << "numberOfRows = " << numberOfRows << ", matrixElements = ";
    for( std::vector< PolynomialSum >::const_iterator
         whichElement( matrixElements.begin() );
         whichElement < matrixElements.end();
         ++whichElement )
    {
      returnStream
      << std::endl << whichElement->AsString()
      << std::endl;
    }
    returnStream
    << "valuesMatrix = " << std::endl << valuesMatrix << std::endl
    << "eigenvalueFinder = ...";
    return std::string( returnStream.str() );
  }

} /* namespace VevaciousPlusPlus */
#endif /* REALMASSESSQUAREDMATRIX_HPP_ */
