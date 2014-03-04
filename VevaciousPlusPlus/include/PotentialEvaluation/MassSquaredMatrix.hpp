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
#include "PolynomialSum.hpp"


namespace VevaciousPlusPlus
{
  class MassSquaredMatrix
  {
  public:
    enum SpinType
    {
      scalarBoson,
      weylFermion,
      gaugeBoson,
      notSet
    };

    MassSquaredMatrix(
                    std::map< std::string, std::string > const& attributeMap );
    virtual
    ~MassSquaredMatrix();


    // This returns the eigenvalues of the matrix.
    std::vector< double > const&
    MassesSquared( std::vector< double > const& fieldConfiguration );

    // This returns the number of identical copies of this mass-squared matrix
    // that the model has.
    double MultiplicityFactor() const{ return multiplicityFactor; }

    SpinType GetSpinType() const{ return spinType; }

    // This adds a new element and returns a reference to it.
    PolynomialSum& AddNewElement();


  protected:
    std::vector< PolynomialSum > matrixElements;
    Eigen::MatrixXd eigenMatrix;
    std::vector< double > massesSquared;
    double multiplicityFactor;
    SpinType spinType;
  };




  inline PolynomialSum& MassSquaredMatrix::AddNewElement()
  {
    matrixElements.push_back( PolynomialSum() );
    return matrixElements.back();
  }

} /* namespace VevaciousPlusPlus */
#endif /* MASSSQUAREDMATRIX_HPP_ */
