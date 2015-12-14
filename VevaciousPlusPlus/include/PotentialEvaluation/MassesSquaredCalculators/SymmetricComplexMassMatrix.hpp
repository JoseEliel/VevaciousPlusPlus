/*
 * SymmetricComplexMassMatrix.hpp
 *
 *  Created on: Oct 09, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef COMPLEXMASSMATRIX_HPP_
#define COMPLEXMASSMATRIX_HPP_

#include "BaseComplexMassMatrix.hpp"
#include <cstddef>
#include <map>
#include <string>
#include "Eigen/Dense"
#include <vector>

namespace VevaciousPlusPlus
{

  class SymmetricComplexMassMatrix : public BaseComplexMassMatrix
  {
  public:
    SymmetricComplexMassMatrix( size_t const numberOfElements,
                    std::map< std::string, std::string > const& attributeMap );
    SymmetricComplexMassMatrix( SymmetricComplexMassMatrix const& copySource );
    SymmetricComplexMassMatrix();
    virtual ~SymmetricComplexMassMatrix();


  protected:
    // This returns a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, using the values for the
    // Lagrangian parameters found in parameterValues.
    virtual Eigen::MatrixXcd
    CurrentValues( std::vector< double > const& parameterValues,
                   std::vector< double > const& fieldConfiguration ) const
    { return LowerTriangleOfSquareMatrix( MatrixToSquare( parameterValues,
                                                      fieldConfiguration ) ); }

    // This returns a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, using the values for the
    // Lagrangian parameters from the last call of UpdateForFixedScale.
    virtual Eigen::MatrixXcd
    CurrentValues( std::vector< double > const& fieldConfiguration ) const
    { return
      LowerTriangleOfSquareMatrix( MatrixToSquare( fieldConfiguration ) ); }

    // This returns a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, using the values for the
    // Lagrangian parameters found in parameterValues.
    Eigen::MatrixXcd
    MatrixToSquare( std::vector< double > const& parameterValues,
                    std::vector< double > const& fieldConfiguration ) const;

    // This returns a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, using the values for the
    // Lagrangian parameters from the last call of UpdateForFixedScale.
    Eigen::MatrixXcd
    MatrixToSquare( std::vector< double > const& fieldConfiguration ) const;

    // This returns a matrix that is the lower-triangular part (only column
    // index <= row index) of the square of matrixToSquare.
    Eigen::MatrixXcd LowerTriangleOfSquareMatrix(
                                Eigen::MatrixXcd const& matrixToSquare ) const;
  };

} /* namespace VevaciousPlusPlus */
#endif /* COMPLEXMASSMATRIX_HPP_ */
