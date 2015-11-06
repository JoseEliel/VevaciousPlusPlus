/*
 * RealMassesSquaredMatrix.hpp
 *
 *  Created on: Oct 09, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef REALMASSESSQUAREDMATRIX_HPP_
#define REALMASSESSQUAREDMATRIX_HPP_

#include "CommonIncludes.hpp"
#include "BasicFunctions/ParametersAndFieldsProductSum.hpp"
#include "Eigen/Dense"
#include "MassesSquaredFromMatrix.hpp"

namespace VevaciousPlusPlus
{

  class RealMassesSquaredMatrix : public MassesSquaredFromMatrix< double >
  {
  public:
    RealMassesSquaredMatrix( size_t const numberOfRows,
                    std::map< std::string, std::string > const& attributeMap );
    RealMassesSquaredMatrix( RealMassesSquaredMatrix const& copySource );
    RealMassesSquaredMatrix();
    virtual ~RealMassesSquaredMatrix();


    // This calls UpdateForFixedScale on each element of matrixElements.
    virtual void
    UpdateForFixedScale( std::vector< double > const& parameterValues );

    // This allows access to the polynomial sum for a given index.
    ParametersAndFieldsProductSum& ElementAt( size_t const elementIndex )
    { return matrixElements[ elementIndex ]; }

    std::vector< ParametersAndFieldsProductSum > const& MatrixElements() const
    { return matrixElements; }


    // This is mainly for debugging:
    std::string AsString() const;


  protected:
    std::vector< ParametersAndFieldsProductSum > matrixElements;

    // This returns a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, using the values for the
    // Lagrangian parameters found in parameterValues.
    virtual Eigen::MatrixXd
    CurrentValues( std::vector< double > const& parameterValues,
                   std::vector< double > const& fieldConfiguration ) const;

    // This returns a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, using the values for the
    // Lagrangian parameters from the last call of UpdateForFixedScale.
    virtual Eigen::MatrixXd
    CurrentValues( std::vector< double > const& fieldConfiguration ) const;
  };





  // This calls UpdateForFixedScale on each element of matrixElements.
  inline void RealMassesSquaredMatrix::UpdateForFixedScale(
                                 std::vector< double > const& parameterValues )
  {
    for( std::vector< ParametersAndFieldsProductSum >::iterator
         parametersAndFieldsProduct( matrixElements.begin() );
         parametersAndFieldsProduct < matrixElements.end();
         ++parametersAndFieldsProduct )
    {
      parametersAndFieldsProduct->UpdateForFixedScale( parameterValues );
    }
  }

  // This is mainly for debugging:
  inline std::string RealMassesSquaredMatrix::AsString() const
  {
    std::stringstream returnStream;
    returnStream
    << "numberOfRows = " << numberOfRows << ", matrixElements = ";
    for( std::vector< ParametersAndFieldsProductSum >::const_iterator
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
#endif /* REALMASSESSQUAREDMATRIX_HPP_ */
