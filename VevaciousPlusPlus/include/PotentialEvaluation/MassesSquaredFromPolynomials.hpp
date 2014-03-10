/*
 * MassesSquaredFromPolynomials.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MASSSQUAREDMATRIX_HPP_
#define MASSSQUAREDMATRIX_HPP_

#include "../StandardIncludes.hpp"
#include "PolynomialSum.hpp"


namespace VevaciousPlusPlus
{
  class MassesSquaredFromPolynomials
  {
  public:
    enum SpinType
    {
      scalarBoson,
      weylFermion,
      gaugeBoson,
      notSet
    };

    MassesSquaredFromPolynomials(
                    std::map< std::string, std::string > const& attributeMap );
    MassesSquaredFromPolynomials(
                              MassesSquaredFromPolynomials const& copySource );
    MassesSquaredFromPolynomials();
    virtual
    ~MassesSquaredFromPolynomials();


    // This should return the eigenvalues of the matrix.
    virtual std::vector< double > const&
    MassesSquared( std::vector< double > const& fieldConfiguration ) = 0;

    // This returns the number of identical copies of this mass-squared matrix
    // that the model has.
    double MultiplicityFactor() const{ return multiplicityFactor; }

    SpinType GetSpinType() const{ return spinType; }


  protected:
    std::vector< double > massesSquared;
    double multiplicityFactor;
    SpinType spinType;
  };

} /* namespace VevaciousPlusPlus */
#endif /* MASSSQUAREDMATRIX_HPP_ */
