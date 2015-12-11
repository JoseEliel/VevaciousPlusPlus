/*
 * BouncePathFinder.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BOUNCEPATHFINDER_HPP_
#define BOUNCEPATHFINDER_HPP_

#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include "BubbleProfile.hpp"
#include "PathParameterization/TunnelPath.hpp"

namespace VevaciousPlusPlus
{

  class BouncePathFinder
  {
  public:
    BouncePathFinder() : pathTemperature( 0.0 ) {}

    virtual ~BouncePathFinder() {}


    // This should reset the BouncePathFinder and set pathTemperature
    // appropriately.
    virtual void SetPotentialAndVacuaAndTemperature(
                                    PotentialFunction const& potentialFunction,
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                      double const pathTemperature = 0.0 ) = 0;

    // This should prompt the derived class instance to assess whether the path
    // can be improved based on the bubble profile which resulted from the last
    // call of TryToImprovePath (or from the path given by a previous
    // BouncePathFinder instance), and return true if it might be improved.
    // Internal variables in the derived class should take care of noticing if
    // this object had not yet had its TryToImprovePath called or what the
    // bounce action previous to that given by bubbleFromLastPath is, if
    // relevant to the decision.
    virtual bool
    PathCanBeImproved( BubbleProfile const& bubbleFromLastPath ) = 0;

    // This should adjust the internal representation of the path towards
    // extremizing the bounce action, ideally in an incremental manner so that
    // the process can be broken off early if it is going to be the case that
    // the survival probability will be lower than a threshold, and return a
    // TunnelPath const* for this path (if we allowed ourselves to require a
    // C++11-compliant compiler, this function would return a const
    // unique_ptr< TunnelPath >). The overshoot/undershoot algorithm along the
    // path should ensure that the bounce action calculated on the path is an
    // upper bound on the bounce action at the desired extremum. The bubble
    // profile from the last path is provided in addition for strategies that
    // rely on the radial dependence of the path auxiliary being approximated
    // reasonably by that of the previous path, similarly to the strategy used
    // by CosmoTransitions.
    virtual TunnelPath const* TryToImprovePath( TunnelPath const& lastPath,
                                 BubbleProfile const& bubbleFromLastPath ) = 0;


  protected:
    double pathTemperature;
  };

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEPATHFINDER_HPP_ */
