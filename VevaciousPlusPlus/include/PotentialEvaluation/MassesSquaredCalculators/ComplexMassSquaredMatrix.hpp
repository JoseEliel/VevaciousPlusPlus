/*
 * ComplexMassSquaredMatrix.hpp
 *
 *  Created on: Oct 09, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef COMPLEXMASSSQUAREDMATRIX_HPP_
#define COMPLEXMASSSQUAREDMATRIX_HPP_

#include "../../BasicFunctions/ParametersAndFieldsProductTermSum.hpp"
#include "BaseComplexMassMatrix.hpp"
#include "CommonIncludes.hpp"
#include "Eigen/Dense"
#include "MassesSquaredFromMatrix.hpp"

namespace VevaciousPlusPlus
{

  class ComplexMassSquaredMatrix : public BaseComplexMassMatrix
  {
  public:
    ComplexMassSquaredMatrix( size_t const numberOfElements,
                    std::map< std::string, std::string > const& attributeMap );
    ComplexMassSquaredMatrix( ComplexMassSquaredMatrix const& copySource );
    ComplexMassSquaredMatrix();
    virtual ~ComplexMassSquaredMatrix();


  protected:
    // This returns a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, using the values for the
    // Lagrangian parameters found in parameterValues.
    virtual Eigen::MatrixXcd
    CurrentValues( std::vector< double > const& parameterValues,
                   std::vector< double > const& fieldConfiguration ) const;

    // This returns a matrix of the values of the elements for a field
    // configuration given by fieldConfiguration, using the values for the
    // Lagrangian parameters from the last call of UpdateForFixedScale.
    virtual Eigen::MatrixXcd
    CurrentValues( std::vector< double > const& fieldConfiguration ) const;
  };

} /* namespace VevaciousPlusPlus */
#endif /* COMPLEXMASSSQUAREDMATRIX_HPP_ */
