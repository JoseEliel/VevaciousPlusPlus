/*
 * RealMassesSquaredMatrix.hpp
 *
 *  Created on: Mar 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef OLDREALMASSESSQUAREDMATRIX_HPP_
#define OLDREALMASSESSQUAREDMATRIX_HPP_

#include "CommonIncludes.hpp"
#include "BasicFunctions/PolynomialSum.hpp"
#include "Eigen/Dense"
#include "OldMassesSquaredFromMatrix.hpp"

namespace VevaciousPlusPlus
{

  class OldRealMassesSquaredMatrix :
                                    public OldMassesSquaredFromMatrix< double >
  {
  public:
    OldRealMassesSquaredMatrix( size_t const numberOfRows,
                    std::map< std::string, std::string > const& attributeMap );
    OldRealMassesSquaredMatrix( OldRealMassesSquaredMatrix const& copySource );
    OldRealMassesSquaredMatrix();
    virtual ~OldRealMassesSquaredMatrix();


    // This allows access to the polynomial sum for a given index.
    PolynomialSum& ElementAt( size_t const elementIndex )
    { return matrixElements[ elementIndex ]; }

    std::vector< PolynomialSum > const& MatrixElements() const
    { return matrixElements; }


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
  inline std::string OldRealMassesSquaredMatrix::AsString() const
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
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* OLDREALMASSESSQUAREDMATRIX_HPP_ */
