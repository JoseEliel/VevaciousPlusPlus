/*
 * MinuitPathPotentialMinimizer.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITPATHPOTENTIALMINIMIZER_HPP_
#define MINUITPATHPOTENTIALMINIMIZER_HPP_

#include "CommonIncludes.hpp"
#include "FullPathVaryingMinuit.hpp"
#include "../PathParameterization/TunnelPathFactory.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"

namespace VevaciousPlusPlus
{

  class MinuitPathPotentialMinimizer : public FullPathVaryingMinuit
  {
  public:
    MinuitPathPotentialMinimizer( TunnelPathFactory* const pathFactory,
                                  PotentialFunction const& potentialFunction,
                                  size_t const numberOfPotentialSamplePoints,
                                  size_t const movesPerImprovement = 100,
                                  unsigned int const minuitStrategy = 1,
                                  double const minuitToleranceFraction = 0.5 );
    virtual ~MinuitPathPotentialMinimizer();


    // It may seem unwise to have this object call Minuit on itself, but really
    // it's just a handy way of keeping the minimization function within the
    // class that ends up finding its minimum. The temperature is set by
    // SetInitialPath.
    virtual double
    operator()( std::vector< double > const& pathParameterization ) const;


  protected:
    PotentialFunction const& potentialFunction;
    size_t const numberOfPotentialSamplePoints;
    size_t const numberOfFields;
    double const pathSegmentSize;
  };




  // It may seem unwise to have this object call Minuit on itself, but really
  // it's just a handy way of keeping the minimization function within the
  // class that ends up finding its minimum. The temperature is set by
  // SetInitialPath.
  inline double MinuitPathPotentialMinimizer::operator()(
                      std::vector< double > const& pathParameterization ) const
  {
    double potentialSum( 0.0 );
    TunnelPath* tunnelPath( (*pathFactory)( pathParameterization,
                                            pathTemperature ) );
    std::vector< double > fieldConfiguration( numberOfFields );
    for( size_t sampleIndex( 0 );
         sampleIndex < numberOfPotentialSamplePoints;
         ++sampleIndex )
    {
      tunnelPath->PutOnPathAt( fieldConfiguration,
                               (double)( sampleIndex + 1 ) * pathSegmentSize );
      potentialSum += potentialFunction( fieldConfiguration,
                                         pathTemperature );
    }
    delete tunnelPath;
    return potentialSum;
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINUITPATHPOTENTIALMINIMIZER_HPP_ */
