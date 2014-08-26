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


    // This resets the BouncePathFinder so that it sets up currentPath as its
    // initial path between the given vacua. It also resets pathCanBeImproved
    // and sets pathTemperature appropriately.
    virtual void SetInitialPath( PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 TunnelPath const* startingPath = NULL,
                                 double const pathTemperature = 0.0 );


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
    if( NanParameterFromMinuit( pathParameterization ) )
    {
      return functionValueForNanInput;
    }
    TunnelPath* tunnelPath( (*pathFactory)( pathParameterization,
                                            pathTemperature ) );
    double const bounceAction( (*bounceActionCalculator)( *tunnelPath ) );
    delete tunnelPath;
    return bounceAction;
  }

  // This resets the BouncePathFinder so that it sets up currentPath as its
  // initial path between the given vacua. It also resets pathCanBeImproved
  // and sets pathTemperature appropriately.
  inline void MinuitPathBounceMinimizer::SetInitialPath(
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                                TunnelPath const* startingPath,
                                                 double const pathTemperature )
  {
    bounceActionCalculator->ResetVacua( falseVacuum,
                                        trueVacuum );
    SetUpPathFactoryAndMinuit( falseVacuum,
                               trueVacuum,
                               startingPath,
                               pathTemperature );
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINUITPATHBOUNCEMINIMIZER_HPP_ */
