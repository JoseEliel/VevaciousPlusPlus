/*
 * MinuitPathFinder.hpp
 *
 *  Created on: Jul 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITPATHFINDER_HPP_
#define MINUITPATHFINDER_HPP_

#include <limits>
#include "CommonIncludes.hpp"
#include "../BouncePathFinder.hpp"
#include "Minuit2/FCNBase.h"

namespace VevaciousPlusPlus
{

  class MinuitPathFinder : public BouncePathFinder,
                           public ROOT::Minuit2::FCNBase
  {
  public:
    MinuitPathFinder( size_t const movesPerImprovement = 100,
                      unsigned int const minuitStrategy = 1,
                      double const minuitToleranceFraction = 0.5 );
    virtual ~MinuitPathFinder();


    // This class implements Up() inherited from ROOT::Minuit2::FCNBase, and
    // its derived classes should implement operator(). There is no deadly
    // diamond of doom because BouncePathFinder does not have an operator() or
    // Up() at all.

    // This implements Up() for ROOT::Minuit2::FCNBase just to stick to a basic
    // value.
    virtual double Up() const { return 1.0; }


  protected:
    static double const functionValueForNanInput;


    size_t movesPerImprovement;
    unsigned int minuitStrategy;
    double minuitToleranceFraction;
    // double currentMinuitTolerance;


    // This returns true if any element of minuitParameterization appears to be
    // a NAN by simultaneously being not >= 0.0 and not < 0.0, false otherwise.
    bool NanParameterFromMinuit(
                   std::vector< double > const& minuitParameterization ) const;
  };




  // This returns true if any element of minuitParameterization appears to be
  // a NAN by simultaneously being not >= 0.0 and not < 0.0, false otherwise.
  inline bool MinuitPathFinder::NanParameterFromMinuit(
                    std::vector< double > const& minuitParameterization ) const
  {
    for( std::vector< double >::const_iterator
         minuitParameter( minuitParameterization.begin() );
         minuitParameter < minuitParameterization.end();
         ++minuitParameter )
    {
      if( !( *minuitParameter >= 0.0 )
          &&
          !( *minuitParameter < 0.0 ) )
      {
        return true;
      }
    }
    return false;
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINUITPATHFINDER_HPP_ */
