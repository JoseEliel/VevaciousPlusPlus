/*
 * ComplexMassSquaredMatrix.hpp
 *
 *  Created on: Mar 31, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef COMPLEXMASSSQUAREDMATRIX_HPP_
#define COMPLEXMASSSQUAREDMATRIX_HPP_

#include "../StandardIncludes.hpp"
#include "MassesSquaredFromPolynomials.hpp"
#include "PolynomialSum.hpp"
#include <complex>
#include "Eigen/Dense"

namespace VevaciousPlusPlus
{

  class ComplexMassSquaredMatrix : public MassesSquaredFromPolynomials
  {
  public:
    ComplexMassSquaredMatrix( unsigned int const numberOfElements,
                    std::map< std::string, std::string > const& attributeMap );
    ComplexMassSquaredMatrix( ComplexMassSquaredMatrix const& copySource );
    ComplexMassSquaredMatrix();
    virtual
    ~ComplexMassSquaredMatrix();


    // This returns the eigenvalues of the square of the matrix.
    virtual std::vector< double > const&
    MassesSquared( std::vector< double > const& fieldConfiguration );

    // This allows access to the pair of polynomial sums for a given index.
    std::pair< PolynomialSum, PolynomialSum >&
    ElementAt( unsigned int const elementIndex )
    { return matrixElements[ elementIndex ]; }

    // This is mainly for debugging:
    std::string AsString() const;

  protected:
    unsigned int numberOfRows;
    std::vector< std::pair< PolynomialSum, PolynomialSum > > matrixElements;
    Eigen::MatrixXcd valuesMatrix;
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXcd > eigenvalueFinder;
  };




  // This is mainly for debugging:
  inline std::string ComplexMassSquaredMatrix::AsString() const
  {
    std::stringstream returnStream;
    returnStream
    << "numberOfRows = " << numberOfRows << ", matrixElements = ";
    for( std::vector< std::pair< PolynomialSum,
                                 PolynomialSum > >::const_iterator
         whichPair( matrixElements.begin() );
         whichPair < matrixElements.end();
         ++whichPair )
    {
      returnStream
      << std::endl << "real:" << std::endl << whichPair->first.AsString()
      << std::endl << "imag:" << std::endl << whichPair->second.AsString()
      << std::endl;
    }
    returnStream
    << "valuesMatrix = " << std::endl << valuesMatrix
    << std::endl << "eigenvalueFinder = ...";
    return std::string( returnStream.str() );
  }

} /* namespace VevaciousPlusPlus */
#endif /* COMPLEXMASSSQUAREDMATRIX_HPP_ */
