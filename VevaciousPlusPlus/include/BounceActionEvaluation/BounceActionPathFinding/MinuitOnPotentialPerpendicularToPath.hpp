/*
 * MinuitOnPotentialPerpendicularToPath.hpp
 *
 *  Created on: Nov 5, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITONPOTENTIALPERPENDICULARTOPATH_HPP_
#define MINUITONPOTENTIALPERPENDICULARTOPATH_HPP_

#include "CommonIncludes.hpp"
#include "MinuitOnHypersurfaces.hpp"
#include "Minuit2/MnMigrad.h"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include "../PathParameterization/LinearSplineThroughNodes.hpp"

namespace VevaciousPlusPlus
{

  class MinuitOnPotentialPerpendicularToPath : public MinuitOnHypersurfaces
  {
  public:
    MinuitOnPotentialPerpendicularToPath(
                                    PotentialFunction const& potentialFunction,
                                          size_t const numberOfPathSegments,
                                       int const numberOfAllowedWorsenings = 3,
                                         unsigned int const minuitStrategy = 1,
                                  double const minuitToleranceFraction = 0.5 );
    virtual ~MinuitOnPotentialPerpendicularToPath();


    // This takes the bounce action from bubbleFromLastPath and updates
    // numberOfAllowedWorsenings based on whether it was an improvement on the
    // previous path. Then it returns false if too many paths have been tried
    // which just made the action bigger, true otherwise.
    virtual bool PathCanBeImproved( BubbleProfile const& bubbleFromLastPath );

    // This minimizes the potential on hyperplanes which are as close to
    // perpendicular to the path given by lastPath, at numberOfVaryingNodes
    // equally-space points along that path, and the minima on those
    // hyperplanes are then used to construct the returned path. (The bubble
    // profile from the last path is ignored.)
    virtual TunnelPath const* TryToImprovePath( TunnelPath const& lastPath,
                                     BubbleProfile const& bubbleFromLastPath );


  protected:
    int numberOfAllowedWorsenings;
    double bounceBeforeLastPath;
    std::vector< double > lastPathFalseSideNode;
    std::vector< double > lastPathTrueSideNode;


    // This is an empty hook that can be over-ridden to account for the
    // bubble profile from the last path.
    virtual void AccountForBubbleProfileAroundNode( size_t nodeIndex,
                                    BubbleProfile const& bubbleFromLastPath ){}
  };

} /* namespace VevaciousPlusPlus */
#endif /* MINUITONPOTENTIALPERPENDICULARTOPATH_HPP_ */
