/*
 * MinuitPathBounceMinimizer.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITPATHBOUNCEMINIMIZER_HPP_
#define MINUITPATHBOUNCEMINIMIZER_HPP_

#include "CommonIncludes.hpp"
#include "FullPathVaryingMinuit.hpp"
#include "../PathParameterization/TunnelPathFactory.hpp"
#include "../BounceActionCalculator.hpp"

namespace VevaciousPlusPlus
{

  class MinuitPathBounceMinimizer : public FullPathVaryingMinuit
  {
  public:
    MinuitPathBounceMinimizer( TunnelPathFactory* const pathFactory,
                          BounceActionCalculator* const bounceActionCalculator,
                               size_t const movesPerImprovement = 100,
                               unsigned int const minuitStrategy = 1,
                               double const minuitToleranceFraction = 0.5 );
    virtual ~MinuitPathBounceMinimizer();


    // It may seem unwise to have this object call Minuit on itself, but really
    // it's just a handy way of keeping the minimization function within the
    // class that ends up finding its minimum. The temperature is set by
    // SetInitialPath.
    virtual double
    operator()( std::vector< double > const& pathParameterization ) const;


  protected:
    BounceActionCalculator* const bounceActionCalculator;
  };




  // It may seem unwise to have this object call Minuit on itself, but really
  // it's just a handy way of keeping the minimization function within the
  // class that ends up finding its minimum. The temperature is set by
  // SetInitialPath.
  inline double MinuitPathBounceMinimizer::operator()(
                      std::vector< double > const& pathParameterization ) const
  {
    TunnelPath* tunnelPath( (*pathFactory)( pathParameterization,
                                            pathTemperature ) );
    double const bounceAction( (*bounceActionCalculator)( *tunnelPath ) );
    delete tunnelPath;
    return bounceAction;
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINUITPATHBOUNCEMINIMIZER_HPP_ */
