/*
 * MinimizingPotentialOnParallelPlanes.hpp
 *
 *  Created on: Oct 27, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINIMIZINGPOTENTIALONPARALLELPLANES_HPP_
#define MINIMIZINGPOTENTIALONPARALLELPLANES_HPP_

#include "CommonIncludes.hpp"
#include "MinimizingPotentialOnHypersurfaces.hpp"
#include "Eigen/Dense"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include "../PathParameterization/LinearSplineThroughNodes.hpp"

namespace VevaciousPlusPlus
{

  class MinimizingPotentialOnParallelPlanes :
                                      public MinimizingPotentialOnHypersurfaces
  {
  public:
    MinimizingPotentialOnParallelPlanes(
                                    PotentialFunction const& potentialFunction,
                                         size_t const numberOfPathSegments,
                                         unsigned int const minuitStrategy = 1,
                                  double const minuitToleranceFraction = 0.5 );
    virtual ~MinimizingPotentialOnParallelPlanes();

    // This should minimize on a series of hyperplanes all perpendicular to
    // the vector difference of the vacua, but maybe to avoid "zig-zag-iness"
    // from Minuit2 coming to rest on a jittery line reasonably far from the
    // starting path, each plane should have its starting point be (previous
    // node + vector difference of previous node from the node before it), with
    // special cases for the first and last varying nodes, of course.


    // It may seem unwise to have this object call Minuit on itself, but really
    // it's just a handy way of keeping the minimization function within the
    // class that ends up finding its minimum. In this case,
    // nodeParameterization is just the parameterization of the node.
    virtual double
    operator()( std::vector< double > const& nodeParameterization ) const;

    // This just returns true if it has not yet provided a path, or false if it
    // already has, as this class does all it can in a single call of
    // TryToImprovePath.
    virtual bool
    PathCanBeImproved( BubbleProfile const& bubbleFromLastPath ) const
    { return notYetProvidedPath; }

    // This minimizes the potential on a series of hyperplanes, all
    // perpendicular to the vector difference of the vacua. It tries to avoid
    // "zig-zag-iness" from Minuit2 coming to rest on a jittery line reasonably
    // far from the starting path by setting the starting point on each plane
    // to be the previous node plus the vector difference of previous node from
    // the node before it, with special cases for the first and last varying
    // nodes, of course. It ignores both arguments, and also sets
    // notYetProvidedPath to false.
    virtual TunnelPath const* TryToImprovePath( TunnelPath const& lastPath,
                                     BubbleProfile const& bubbleFromLastPath );


  protected:
    bool notYetProvidedPath;
    double const planeDifferenceFraction;
  };

} /* namespace VevaciousPlusPlus */
#endif /* MINIMIZINGPOTENTIALONPARALLELPLANES_HPP_ */
