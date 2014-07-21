/*
 * PolynomialsFromCoefficientsFactory.hpp
 *
 *  Created on: Jul 21, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALSFROMCOEFFICIENTSFACTORY_HPP_
#define POLYNOMIALSFROMCOEFFICIENTSFACTORY_HPP_

#include "CommonIncludes.hpp"
#include "TunnelPathFactory.hpp"
#include "TunnelPath.hpp"
#include "PolynomialPathFromCoefficients.hpp"

namespace VevaciousPlusPlus
{

  class PolynomialsFromCoefficientsFactory : public TunnelPathFactory
  {
  public:
    PolynomialsFromCoefficientsFactory( size_t const numberOfFields,
                            size_t const numberOfVaryingCoefficientsPerField );
    virtual ~PolynomialsFromCoefficientsFactory();


    // This should reset the TunnelPathFactory so that it will produce
    // TunnelPath*s that parameterize the path between the given vacua.
    virtual void SetVacua( PotentialMinimum const& falseVacuum,
                           PotentialMinimum const& trueVacuum );

    // This returns a PolynomialPathFromCoefficients* as a TunnelPath*.
    // The calling code is responsible for memory management! (It'd be great to
    // return a std::unique_Ptr< TunnelPath >, but we're stubbornly sticking to
    // the requirements not including a C++11-compliant compiler.)
    virtual TunnelPath*
    operator()( std::vector< double > const& pathParameterization,
                double const pathTemperature = 0.0 ) const
    { return new PolynomialPathFromCoefficients( numberOfFields,
                                                 falseVacuum,
                                                 vacuumDifference,
                                                 pathParameterization,
                                                 pathTemperature,
                                                 referenceFieldIndex,
                                                 coefficientsBeyondLinear ); }


  protected:
    std::vector< double > falseVacuum;
    std::vector< double > vacuumDifference;
    size_t referenceFieldIndex;
    size_t const coefficientsBeyondLinear;
  };




  // This should reset the TunnelPathFactory so that it will produce
  // TunnelPath*s that parameterize the path between the given vacua.
  inline void PolynomialsFromCoefficientsFactory::SetVacua(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    this->falseVacuum = falseVacuum.FieldConfiguration();
    referenceFieldIndex = 0;
    double largestFieldDifference( trueVacuum.FieldConfiguration().front()
                                  - falseVacuum.FieldConfiguration().front() );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      vacuumDifference[ fieldIndex ]
      = ( trueVacuum.FieldConfiguration()[ fieldIndex ]
          - falseVacuum.FieldConfiguration()[ fieldIndex ] );
      if( vacuumDifference[ fieldIndex ] > largestFieldDifference )
      {
        referenceFieldIndex = fieldIndex;
        largestFieldDifference = vacuumDifference[ fieldIndex ];
      }
    }
    initialStepSizes.assign( zeroParameterization.size(),
                             ( 0.5 * largestFieldDifference ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALSFROMCOEFFICIENTSFACTORY_HPP_ */
