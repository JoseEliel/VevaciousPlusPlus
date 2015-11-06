/*
 * ParametersAndFieldsProductSum.hpp
 *
 *  Created on: Oct 9, 2015
 *      Author: bol
 */

#ifndef PARAMETERSANDFIELDSPRODUCTSUM_HPP_
#define PARAMETERSANDFIELDSPRODUCTSUM_HPP_

#include "CommonIncludes.hpp"
#include "ParametersAndFieldsProduct.hpp"

namespace VevaciousPlusPlus
{

  class ParametersAndFieldsProductSum
  {
  public:
    ParametersAndFieldsProductSum();
    ParametersAndFieldsProductSum(
                             ParametersAndFieldsProductSum const& copySource );
    virtual ~ParametersAndFieldsProductSum();


    // This calls UpdateForFixedScale on each element of
    // parametersAndFieldsProducts.
    void UpdateForFixedScale( std::vector< double > const& parameterValues );

    // This returns the sum of operator() for each element of
    // parametersAndFieldsProducts.
    double operator()( std::vector< double > const& parameterValues,
                       std::vector< double > const& fieldConfiguration ) const;

    // This returns the sum of operator() for each element of
    // parametersAndFieldsProducts.
    double operator()( std::vector< double > const& fieldConfiguration ) const;

    std::vector< ParametersAndFieldsProduct > const&
    ParametersAndFieldsProducts() const
    { return parametersAndFieldsProducts; }

    std::vector< ParametersAndFieldsProduct >& ParametersAndFieldsProducts()
    { return parametersAndFieldsProducts; }

    // This returns the highest sum of field powers of all the terms in
    // parametersAndFieldsProducts.
    unsigned int HighestFieldPower() const;

    // This is mainly for debugging:
    std::string AsDebuggingString() const;

    // This returns a string that should be valid Python assuming that the
    // field configuration is given as an array called "fv" and that the
    // Lagrangian parameters are in an array called "lp".
    std::string AsPython() const;


  protected:
    std::vector< ParametersAndFieldsProduct > parametersAndFieldsProducts;
  };





  // This calls UpdateForFixedScale on each element of
  // parametersAndFieldsProducts.
  inline void ParametersAndFieldsProductSum::UpdateForFixedScale(
                                 std::vector< double > const& parameterValues )
  {
    for( std::vector< ParametersAndFieldsProduct >::iterator
         parametersAndFieldsProduct( parametersAndFieldsProducts.begin() );
         parametersAndFieldsProduct < parametersAndFieldsProducts.end();
         ++parametersAndFieldsProduct )
    {
      parametersAndFieldsProduct->UpdateForFixedScale( parameterValues );
    }
  }

  // This returns the sum of operator() for each element of
  // parametersAndFieldsProducts.
  inline double ParametersAndFieldsProductSum::operator()(
                                  std::vector< double > const& parameterValues,
                        std::vector< double > const& fieldConfiguration ) const
  {
    double returnSum( 0.0 );
    for( std::vector< ParametersAndFieldsProduct >::const_iterator
         parametersAndFieldsProduct( parametersAndFieldsProducts.begin() );
         parametersAndFieldsProduct < parametersAndFieldsProducts.end();
         ++parametersAndFieldsProduct )
    {
      returnSum += (*parametersAndFieldsProduct)( parameterValues,
                                                  fieldConfiguration );
    }
    return returnSum;
  }

  // This returns the sum of operator() for each element of
  // parametersAndFieldsProducts.
  inline double ParametersAndFieldsProductSum::operator()(
                        std::vector< double > const& fieldConfiguration ) const
  {
    double returnSum( 0.0 );
    for( std::vector< ParametersAndFieldsProduct >::const_iterator
         parametersAndFieldsProduct( parametersAndFieldsProducts.begin() );
         parametersAndFieldsProduct < parametersAndFieldsProducts.end();
         ++parametersAndFieldsProduct )
    {
      returnSum += (*parametersAndFieldsProduct)( fieldConfiguration );
    }
    return returnSum;
  }

  // This returns the highest sum of field powers of all the terms in
  // parametersAndFieldsProducts.
  inline unsigned int ParametersAndFieldsProductSum::HighestFieldPower() const
  {
    unsigned int highestPower( 0 );
    for( std::vector< ParametersAndFieldsProduct >::const_iterator
         parametersAndFieldsProduct( parametersAndFieldsProducts.begin() );
         parametersAndFieldsProduct < parametersAndFieldsProducts.end();
         ++parametersAndFieldsProduct )
    {
      if( parametersAndFieldsProduct->FieldPower() > highestPower )
      {
        highestPower = parametersAndFieldsProduct->FieldPower();
      }
    }
    return highestPower;
  }

  // This is mainly for debugging:
  inline std::string ParametersAndFieldsProductSum::AsDebuggingString() const
  {
    std::stringstream returnStream;
    returnStream << "ParametersAndFieldsProductSum =" << std::endl;
    for( std::vector< ParametersAndFieldsProduct >::const_iterator
         parametersAndFieldsProduct( parametersAndFieldsProducts.begin() );
         parametersAndFieldsProduct < parametersAndFieldsProducts.end();
         ++parametersAndFieldsProduct )
    {
      returnStream
      << parametersAndFieldsProduct->AsDebuggingString() << std::endl;
    }
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */

#endif /* PARAMETERSANDFIELDSPRODUCTSUM_HPP_ */
