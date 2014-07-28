/*
 * MinuitPathFinder.hpp
 *
 *  Created on: Jul 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITPATHFINDER_HPP_
#define MINUITPATHFINDER_HPP_

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
    size_t movesPerImprovement;
    unsigned int minuitStrategy;
    double minuitToleranceFraction;
    double currentMinuitTolerance;
  };

} /* namespace VevaciousPlusPlus */
#endif /* MINUITPATHFINDER_HPP_ */
