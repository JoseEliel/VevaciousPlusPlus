/*
 * SymmetricComplexMassMatrix.hpp
 *
 *  Created on: Mar 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef COMPLEXMASSMATRIX_HPP_
#define COMPLEXMASSMATRIX_HPP_

#include "CommonIncludes.hpp"
#include "MassesSquaredFromMatrix.hpp"
#include "BasicFunctions/PolynomialSum.hpp"
#include "Eigen/Dense"

namespace VevaciousPlusPlus
{

  class SymmetricComplexMassMatrix :
                       public MassesSquaredFromMatrix< std::complex< double > >
  {
  public:
    SymmetricComplexMassMatrix( size_t const numberOfElements,
                    std::map< std::string, std::string > const& attributeMap );
    SymmetricComplexMassMatrix( SymmetricComplexMassMatrix const& copySource );
    SymmetricComplexMassMatrix();
    virtual
    ~SymmetricComplexMassMatrix();


    // This allows access to the pair of polynomial sums for a given index.
    std::pair< PolynomialSum, PolynomialSum >&
    ElementAt( size_t const elementIndex )
    { return matrixElements[ elementIndex ]; }

    std::vector< std::pair< PolynomialSum, PolynomialSum > > const&
    MatrixElements() const{ return matrixElements; }

    // This is mainly for debugging:
    std::string AsString() const;


  protected:
    std::vector< std::pair< PolynomialSum, PolynomialSum > > matrixElements;

    // This returns a matrix of the values of the square of the matrix of
    // elements for a field configuration given by fieldConfiguration, with all
    // functionoids evaluated at the last scale which was used to update them.
    virtual Eigen::MatrixXcd
    CurrentValues( std::vector< double > const& fieldConfiguration ) const
    { return
      LowerTriangleOfSquareMatrix( MatrixToSquare( fieldConfiguration ) ); }

    // This returns a matrix of the values of the square of the matrix of
    // elements for a field configuration given by fieldConfiguration, with all
    // functionoids evaluated at the natural exponent of logarithmOfScale.
    virtual Eigen::MatrixXcd
    CurrentValues( std::vector< double > const& fieldConfiguration,
                   double const logarithmOfScale ) const
    { return LowerTriangleOfSquareMatrix( MatrixToSquare( fieldConfiguration,
                                                        logarithmOfScale ) ); }

    // This returns a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, with all functionoids
    // evaluated at the last scale which was used to update them.
    Eigen::MatrixXcd
    MatrixToSquare( std::vector< double > const& fieldConfiguration ) const;

    // This returns a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, with all functionoids
    // evaluated at the natural exponent of logarithmOfScale.
    Eigen::MatrixXcd
    MatrixToSquare( std::vector< double > const& fieldConfiguration,
                    double const logarithmOfScale ) const;

    // This returns a matrix that is the lower-triangular part (only column
    // index <= row index) of the square of matrixToSquare.
    Eigen::MatrixXcd LowerTriangleOfSquareMatrix(
                                Eigen::MatrixXcd const& matrixToSquare ) const;
  };




  // This is mainly for debugging:
  inline std::string SymmetricComplexMassMatrix::AsString() const
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
      << std::endl << "real:" << std::endl
      << whichPair->first.AsDebuggingString()
      << std::endl << "imag:" << std::endl
      << whichPair->second.AsDebuggingString()
      << std::endl;
    }
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* COMPLEXMASSMATRIX_HPP_ */
