/*
 * MinuitOnPotentialPerpendicularToPath.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinuitOnPotentialPerpendicularToPath.hpp"

namespace VevaciousPlusPlus
{

  MinuitOnPotentialPerpendicularToPath::MinuitOnPotentialPerpendicularToPath(
                                    PotentialFunction const& potentialFunction,
                                             size_t const numberOfPathSegments,
                                           int const numberOfAllowedWorsenings,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
     MinuitOnHypersurfaces( potentialFunction,
                            numberOfPathSegments,
                            minuitStrategy,
                            minuitToleranceFraction ),
     numberOfAllowedWorsenings(
                             static_cast< int >( numberOfAllowedWorsenings ) ),
     bounceBeforeLastPath( NAN )
  {
    // This constructor is just an initialization list.
  }

  MinuitOnPotentialPerpendicularToPath::~MinuitOnPotentialPerpendicularToPath()
  {
    // This does nothing.
  }

  // This takes the bounce action from bubbleFromLastPath and updates
  // numberOfAllowedWorsenings based on whether it was an improvement on the
  // previous path. Then it returns false if too many paths have been tried
  // which just made the action bigger, true otherwise.
  bool MinuitOnPotentialPerpendicularToPath::PathCanBeImproved(
                                      BubbleProfile const& bubbleFromLastPath )
  {
    if( bounceBeforeLastPath <= bubbleFromLastPath.BounceAction() )
    {
      --numberOfAllowedWorsenings;
    }
    bounceBeforeLastPath = bubbleFromLastPath.BounceAction();
    return ( numberOfAllowedWorsenings > 0 );
  }

  // This minimizes the potential on hyperplanes which are as close to
  // perpendicular to the path given by lastPath, at numberOfVaryingNodes
  // equally-space points along that path, and the minima on those
  // hyperplanes are then used to construct the returned path. (The bubble
  // profile from the last path is ignored.)
  TunnelPath const* MinuitOnPotentialPerpendicularToPath::TryToImprovePath(
                                                    TunnelPath const& lastPath,
                                      BubbleProfile const& bubbleFromLastPath )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitOnPotentialPerpendicularToPath::TryToImprovePath(...) returning"
    << " NULL";
    std::cout << std::endl;
    return NULL;/**/
  }

} /* namespace VevaciousPlusPlus */
