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


    // This should reset the BouncePathFinder so that it sets up currentPath as
    // its initial path between the given vacua. It should also reset
    // pathCanBeImproved and set pathTemperature appropriately.
    virtual void SetInitialPath( PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 TunnelPath const* startingPath = NULL,
                                 double const pathTemperature = 0.0 ) = 0;

    // This returns what should be the best path so far.
    virtual TunnelPath const& CurrentPath() const{ return *currentPath; }

    // This returns true if the last reset or path improvement did not meet
    // the criterion the derived class has that the path cannot be improved any
    // more.
    bool PathCanBeImproved() const{ return pathCanBeImproved; }

    // This should adjust the path towards extremizing the bounce action,
    // ideally in an incremental manner so that the process can be broken off
    // early if it is going to be the case that the survival probability will
    // be lower than a threshold. The overshoot/undershoot algorithm along the
    // path should ensure that the bounce action calculated on the path is an
    // upper bound on the bounce action at the desired extremum.
    virtual void ImprovePath() = 0;


  protected:
    double pathTemperature;
    bool pathCanBeImproved;


    // This would be unnecessary if we were to require a C++11-compliant
    // compiler, because then currentPath could be a
    // std::unique_ptr< TunnelPath > and memory management would be automatic.
    void SetCurrentPathPointer( TunnelPath* currentPath )
    { delete this->currentPath; this->currentPath = currentPath; }


  private:
    TunnelPath* currentPath;
  };

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEPATHFINDER_HPP_ */
