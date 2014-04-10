/*
 * RealMassesSquaredMatrix.hpp
 *
 *  Created on: Mar 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef REALMASSESSQUAREDMATRIX_HPP_
#define REALMASSESSQUAREDMATRIX_HPP_

#include "../StandardIncludes.hpp"
#include "MassesSquaredFromMatrix.hpp"
#include "PolynomialSum.hpp"
#include "Eigen/Dense"

namespace VevaciousPlusPlus
{

  class RealMassesSquaredMatrix : public MassesSquaredFromMatrix< double >
  {
  public:
    RealMassesSquaredMatrix( unsigned int const numberOfRows,
                    std::map< std::string, std::string > const& attributeMap );
    RealMassesSquaredMatrix( RealMassesSquaredMatrix const& copySource );
    RealMassesSquaredMatrix();
    virtual
    ~RealMassesSquaredMatrix();


    // This allows access to the polynomial sum for a given index.
    PolynomialSum& ElementAt( unsigned int const elementIndex )
    { return matrixElements[ elementIndex ]; }

    // This is mainly for debugging:
    std::string AsString() const;

  protected:
    std::vector< PolynomialSum > matrixElements;

    // This returns a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, with all functionoids
    // evaluated at the last scale which was used to update them.
    virtual Eigen::MatrixXd
    CurrentValues( std::vector< double > const& fieldConfiguration ) const;

    // This returns a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, with all functionoids
    // evaluated at the natural exponent of logarithmOfScale.
    virtual Eigen::MatrixXd
    CurrentValues( std::vector< double > const& fieldConfiguration,
                   double const logarithmOfScale ) const;
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
      << std::endl << whichElement->AsDebuggingString()
      << std::endl;
    }
    return std::string( returnStream.str() );
  }

} /* namespace VevaciousPlusPlus */
#endif /* REALMASSESSQUAREDMATRIX_HPP_ */
