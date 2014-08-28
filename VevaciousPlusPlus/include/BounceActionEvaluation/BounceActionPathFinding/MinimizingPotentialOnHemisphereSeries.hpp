/*
 * MinimizingPotentialOnHemisphereSeries.hpp
 *
 *  Created on: Aug 27, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINIMIZINGPOTENTIALONHEMISPHERESERIES_HPP_
#define MINIMIZINGPOTENTIALONHEMISPHERESERIES_HPP_

#include "CommonIncludes.hpp"
#include "MinuitPathFinder.hpp"
#include "Minuit2/MnMigrad.h"
#include "PotentialMinimization/GradientBasedMinimization/MinuitMinimum.hpp"
#include "../PathParameterization/NodesOnBisectingPlanes.hpp"
#include "../PathParameterization/LinearSplineThroughNodesFactory.hpp"

namespace VevaciousPlusPlus
{

  class MinimizingPotentialOnHemisphereSeries : public MinuitPathFinder
  {
  public:
    MinimizingPotentialOnHemisphereSeries();
    virtual ~MinimizingPotentialOnHemisphereSeries();


    // This resets the BouncePathFinder so that it sets up currentPath as its
    // initial path between the given vacua. It also resets pathCanBeImproved
    // and sets pathTemperature appropriately.
    virtual void SetInitialPath( PotentialMinimum const& falseVacuum,
                                 PotentialMinimum const& trueVacuum,
                                 TunnelPath const* startingPath = NULL,
                                 double const pathTemperature = 0.0 )
    { SetUpPathFactoryAndMinuit( falseVacuum,
                                 trueVacuum,
                                 startingPath,
                                 pathTemperature ); }

    // This allows Minuit2 to adjust the full path a set number of times to try
    // to minimize the sum of potentials at a set of nodes or bounce action
    // along the adjusted path, and then sets the path.
    virtual void ImprovePath();


  protected:
    LinearSplineThroughNodesFactory pathFactory;
    std::vector< std::vector< double > > pathNodes;


    // This sets up pathFactory, pathTemperature, currentMinuitTolerance, and
    // currentMinuitResult to be appropriate for the initial path.
    void SetUpPathFactoryAndMinuit( PotentialMinimum const& falseVacuum,
                                    PotentialMinimum const& trueVacuum,
                                    TunnelPath const* startingPath = NULL,
                                    double const pathTemperature = 0.0 );
  };

} /* namespace VevaciousPlusPlus */
#endif /* MINIMIZINGPOTENTIALONHEMISPHERESERIES_HPP_ */
