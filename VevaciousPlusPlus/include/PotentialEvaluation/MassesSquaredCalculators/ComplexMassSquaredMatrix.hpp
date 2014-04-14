/*
 * ComplexMassSquaredMatrix.hpp
 *
 *  Created on: Mar 31, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef COMPLEXMASSSQUAREDMATRIX_HPP_
#define COMPLEXMASSSQUAREDMATRIX_HPP_

#include "../../StandardIncludes.hpp"
#include "MassesSquaredFromMatrix.hpp"
#include "../PolynomialSum.hpp"
#include "Eigen/Dense"

namespace VevaciousPlusPlus
{

  class ComplexMassSquaredMatrix :
                       public MassesSquaredFromMatrix< std::complex< double > >
  {
  public:
    ComplexMassSquaredMatrix( unsigned int const numberOfElements,
                    std::map< std::string, std::string > const& attributeMap );
    ComplexMassSquaredMatrix( ComplexMassSquaredMatrix const& copySource );
    ComplexMassSquaredMatrix();
    virtual
    ~ComplexMassSquaredMatrix();


    // This allows access to the pair of polynomial sums for a given index.
    std::pair< PolynomialSum, PolynomialSum >&
    ElementAt( unsigned int const elementIndex )
    { return matrixElements[ elementIndex ]; }

    // This is mainly for debugging:
    std::string AsString() const;


  protected:
    std::vector< std::pair< PolynomialSum, PolynomialSum > > matrixElements;

    // This returns a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, with all functionoids
    // evaluated at the last scale which was used to update them.
    virtual Eigen::MatrixXcd
    CurrentValues( std::vector< double > const& fieldConfiguration ) const;

    // This returns a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, with all functionoids
    // evaluated at the natural exponent of logarithmOfScale.
    virtual Eigen::MatrixXcd
    CurrentValues( std::vector< double > const& fieldConfiguration,
                   double const logarithmOfScale ) const;
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
      << std::endl << "real:" << std::endl
      << whichPair->first.AsDebuggingString()
      << std::endl << "imag:" << std::endl
      << whichPair->second.AsDebuggingString()
      << std::endl;
    }
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* COMPLEXMASSSQUAREDMATRIX_HPP_ */
