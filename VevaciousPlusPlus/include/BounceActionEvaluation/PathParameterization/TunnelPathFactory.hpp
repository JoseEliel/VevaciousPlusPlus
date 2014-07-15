/*
 * TunnelPathFactory.hpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef TUNNELPATHFACTORY_HPP_
#define TUNNELPATHFACTORY_HPP_

#include "CommonIncludes.hpp"
#include "TunnelPath.hpp"

namespace VevaciousPlusPlus
{

  class TunnelPathFactory
  {
  public:
    TunnelPathFactory( size_t const numberOfFields,
                       std::vector< double > const& zeroParameterization );
    virtual ~TunnelPathFactory();


    // This should return a pointer to an appropriate derived class instance.
    // The calling code is responsible for memory management! (It'd be great to
    // return a std::unique_Ptr< TunnelPath >, but we're stubbornly sticking to
    // the requirements not including a C++11-compliant compiler.)
    virtual TunnelPath*
    operator()( std::vector< double > const& pathParameterization,
                double const pathTemperature = 0.0 ) const = 0;


    virtual std::vector< double > const&
    ZeroParameterization() const{ return zeroParameterization; }

    virtual std::vector< double > const&
    InitialStepSizes() const{ return initialStepSizes; }


  protected:
    size_t const numberOfFields;
    std::vector< double > zeroParameterization;
    std::vector< double > initialStepSizes;
  };

} /* namespace VevaciousPlusPlus */
#endif /* TUNNELPATHFACTORY_HPP_ */
