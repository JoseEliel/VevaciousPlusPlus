/*
 * MinuitOnPotentialOnParallelPlanes.hpp
 *
 *  Created on: Oct 27, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITONPOTENTIALONPARALLELPLANES_HPP_
#define MINUITONPOTENTIALONPARALLELPLANES_HPP_

#include "CommonIncludes.hpp"
#include "MinuitOnHypersurfaces.hpp"
#include "Minuit2/MnMigrad.h"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include "../PathParameterization/LinearSplineThroughNodes.hpp"

namespace VevaciousPlusPlus
{

  class MinuitOnPotentialOnParallelPlanes : public MinuitOnHypersurfaces
  {
  public:
    MinuitOnPotentialOnParallelPlanes( size_t const numberOfPathSegments,
                                       unsigned int const minuitStrategy = 1,
                                  double const minuitToleranceFraction = 0.5 );
    virtual ~MinuitOnPotentialOnParallelPlanes();


    // This just returns false, as this class does all it can in a single call
    // of TryToImprovePath, which is called before PathCanBeImproved under
    // the only circumstances assumed.
    virtual bool PathCanBeImproved( BubbleProfile const& bubbleFromLastPath )
    { return false; }

    // This minimizes the potential on a series of hyperplanes, all
    // perpendicular to the vector difference of the vacua. It tries to avoid
    // "zig-zag-iness" from Minuit2 coming to rest on a jittery line reasonably
    // far from the starting path by setting the starting point on each plane
    // to be the previous node plus the vector difference of previous node from
    // the node before it, with a special case for the first varying node. It
    // ignores both arguments, and also sets notYetProvidedPath to false.
    virtual TunnelPath const* TryToImprovePath( TunnelPath const& lastPath,
                                     BubbleProfile const& bubbleFromLastPath );
  };

} /* namespace VevaciousPlusPlus */
#endif /* MINUITONPOTENTIALONPARALLELPLANES_HPP_ */
