/*
 * PolynomialTerm.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALTERM_HPP_
#define POLYNOMIALTERM_HPP_

#include "../StandardIncludes.hpp"
#include "ParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class PolynomialTerm
  {
  public:
    PolynomialTerm();
    PolynomialTerm( PolynomialTerm const& copySource );
    virtual
    ~PolynomialTerm();


    // This multiplies the relevant field values with the coefficient and the
    // values from the functionoids.
    double operator()( std::vector< double > const& fieldConfiguration ) const;

    // This returns a flag that can be false if the coefficient got multiplied
    // by a NULL functionoid.
    bool IsValid() const { return isValid; }

    // This raises the power of the field given by fieldIndex by the number
    // given by powerInt.
    void RaiseFieldPower( unsigned int const fieldIndex,
                          unsigned int const powerInt );

    // This multiplies coefficientConstant by multiplicationFactor.
    void MultiplyBy( double const multiplicationFactor );

    // This adds runningParameter to the set of functionoids which multiply
    // coefficientConstant to form the scale-dependent coefficient.
    void MultiplyBy( ParameterFunctionoid* const runningParameter,
                     unsigned int const powerInt );



  protected:
    bool isValid;
    double coefficientConstant;
    std::vector< unsigned int > fieldProductByIndex;
    std::vector< ParameterFunctionoid* > functionoidProduct;
  };




  // This multiplies the relevant field values with the coefficient and the
  // values from the functionoids.
  inline double PolynomialTerm::operator()(
                        std::vector< double > const& fieldConfiguration ) const
  {
    double returnValue( coefficientConstant );
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "PolynomialTerm::operator()(";
    for( std::vector< double >::const_iterator
         whichField( fieldConfiguration.begin() );
         whichField < fieldConfiguration.end();
         ++whichField )
    {
      std::cout << ' ' << *whichField;
    }
    std::cout
    << " ) called. returnValue = " << returnValue << ", coefficientConstant = "
    << coefficientConstant << std::endl;/**/

    for( std::vector< unsigned int >::const_iterator
         whichField( fieldProductByIndex.begin() );
         whichField < fieldProductByIndex.end();
         ++whichField )
    {
      returnValue *= fieldConfiguration[ *whichField ];

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "fieldConfiguration[ *whichField ] = "
      << fieldConfiguration[ *whichField ] << ", returnValue = "
      << returnValue;
      std::cout << std::endl;/**/
    }
    for( std::vector< ParameterFunctionoid* >::const_iterator
         whichFunctionoid( functionoidProduct.begin() );
         whichFunctionoid < functionoidProduct.end();
         ++whichFunctionoid )
    {
      returnValue *= (*(*whichFunctionoid))();

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "(*(*whichFunctionoid))() = " << (*(*whichFunctionoid))()
      << ", returnValue = " << returnValue;
      std::cout << std::endl;/**/
    }
    return returnValue;
  }

  // This multiplies coefficientConstant by multiplicationFactor.
  inline void PolynomialTerm::MultiplyBy( double const multiplicationFactor )
  {
    coefficientConstant *= multiplicationFactor;
  }

  // This raises the power of the field given by fieldIndex by the number
  // given by powerInt.
  inline void PolynomialTerm::RaiseFieldPower( unsigned int const fieldIndex,
                                               unsigned int const powerInt )
  {
    for( unsigned int powerCount( 0 );
         powerCount < powerInt;
         ++powerCount )
    {
      fieldProductByIndex.push_back( fieldIndex );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALTERM_HPP_ */
