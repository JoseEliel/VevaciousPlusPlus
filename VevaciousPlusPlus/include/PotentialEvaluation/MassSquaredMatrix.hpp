/*
 * MassSquaredMatrix.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MASSSQUAREDMATRIX_HPP_
#define MASSSQUAREDMATRIX_HPP_

#include "../StandardIncludes.hpp"
#include "Eigen/Dense"

namespace VevaciousPlusPlus
{
  class MassSquaredMatrix
  {
  public:
    MassSquaredMatrix(
                    std::map< std::string, std::string > const& attributeMap );
    virtual
    ~MassSquaredMatrix();


    // This adds a new element and returns a reference to it.
    PolynomialSum& AddNewElement();


  protected:
    std::vector< PolynomialSum > matrixElements;
    Eigen::MatrixXd eigenMatrix;
  };




  inline PolynomialSum& MassSquaredMatrix::AddNewElement()
  {
    matrixElements.push_back( PolynomialSum() );
    return matrixElements.back();
  }

} /* namespace VevaciousPlusPlus */
#endif /* MASSSQUAREDMATRIX_HPP_ */
