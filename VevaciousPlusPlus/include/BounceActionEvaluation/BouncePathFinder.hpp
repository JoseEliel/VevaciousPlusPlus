/*
 * BouncePathFinder.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BOUNCEPATHFINDER_HPP_
#define BOUNCEPATHFINDER_HPP_

#include "CommonIncludes.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include "PathParameterization/TunnelPath.hpp"

namespace VevaciousPlusPlus
{

  class BouncePathFinder
  {
  public:
    BouncePathFinder();
    virtual ~BouncePathFinder();


    // This should reset the BouncePathFinder and return a pointer to an
    // initial path between the given vacua, and the management of this memory
    // is entirely up to the calling code (if we allowed ourselves to require a
    // C++11-compliant compiler, this function would return a
    // unique_ptr< TunnelPath >). It should also reset
    // pathCanBeImproved and set pathTemperature appropriately.
    virtual void SetVacuaAndTemperature( PotentialMinimum const& falseVacuum,
                                         PotentialMinimum const& trueVacuum,
                                      double const pathTemperature = 0.0 ) = 0;

    // This returns true if the last reset or path improvement did not meet
    // the criterion the derived class has that the path cannot be improved any
    // more.
    bool PathCanBeImproved() const{ return pathCanBeImproved; }

    // This should adjust the internal representation of the path towards
    // extremizing the bounce action, ideally in an incremental manner so that
    // the process can be broken off early if it is going to be the case that
    // the survival probability will be lower than a threshold, and return a
    // TunnelPath* for this path (if we allowed ourselves to require a
    // C++11-compliant compiler, this function would return a
    // unique_ptr< TunnelPath >). The overshoot/undershoot algorithm along the
    // path should ensure that the bounce action calculated on the path is an
    // upper bound on the bounce action at the desired extremum. Whether or not
    // the last call managed to lower the bounce action or not is given by
    // lastImprovementWorked, which might be used to decide whether to change
    // path improvement strategy internally in derived classes.
    virtual TunnelPath const*
    TryToImprovePath( bool const lastImprovementWorked = true ) = 0;


  protected:
    double pathTemperature;
    bool pathCanBeImproved;
  };

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEPATHFINDER_HPP_ */
