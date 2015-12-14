/*
 * MinuitOnPotentialPerpendicularToPath.hpp
 *
 *  Created on: Nov 5, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITONPOTENTIALPERPENDICULARTOPATH_HPP_
#define MINUITONPOTENTIALPERPENDICULARTOPATH_HPP_

#include "MinuitOnHypersurfaces.hpp"
#include <vector>
#include "BounceActionEvaluation/BubbleProfile.hpp"
#include "BounceActionEvaluation/PathParameterization/TunnelPath.hpp"
#include "Eigen/Dense"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include <cstddef>
#include "BounceActionEvaluation/PathParameterization/LinearSplineThroughNodes.hpp"

namespace VevaciousPlusPlus
{

  class MinuitOnPotentialPerpendicularToPath : public MinuitOnHypersurfaces
  {
  public:
    MinuitOnPotentialPerpendicularToPath(
                                       unsigned int const numberOfPathSegments,
                                  unsigned int const numberOfAllowedWorsenings,
                                    double const nodeMovementThresholdFraction,
                                          double const dampingFactor,
                       std::vector< double > const neighborDisplacementWeights,
                                          unsigned int const minuitStrategy,
                                        double const minuitToleranceFraction );
    virtual ~MinuitOnPotentialPerpendicularToPath();


    // This returns false if the last call of TryToImprovePath set
    // nodesConverged to be true because the nodes did not move very far.
    // Otherwise, it takes the bounce action from bubbleFromLastPath and
    // updates numberOfAllowedWorsenings based on whether it was an improvement
    // on the previous path. Then it returns false if too many paths have been
    // tried which just made the action bigger, true otherwise.
    virtual bool PathCanBeImproved( BubbleProfile const& bubbleFromLastPath );

    // This minimizes the potential on hyperplanes which are as close to
    // perpendicular to the path given by lastPath, at numberOfVaryingNodes
    // equally-space points along that path, and the minima on those
    // hyperplanes are then used to construct the returned path. If the minima
    // are all sufficiently close to their starting points from lastPath,
    // nodesConverged is set to true. (The bubble profile from the last path is
    // ignored, but there is an empty hook in the loop to allow derived classes
    // to use it.)
    virtual TunnelPath const* TryToImprovePath( TunnelPath const& lastPath,
                                     BubbleProfile const& bubbleFromLastPath );


  protected:
    unsigned int const numberOfAllowedWorsenings;
    unsigned int numberOfWorseningsSoFar;
    double const nodeMovementThresholdFractionSquared;
    double const dampingFactor;
    std::vector< double > const neighborDisplacementWeights;
    std::vector< Eigen::VectorXd > lastPathNodes;
    std::vector< Eigen::VectorXd > nodeDisplacements;
    double bounceBeforeLastPath;
    bool nodesConverged;


    // This sets up nodeDisplacements and lastPathNodes to have the correct
    // sizes, sets returnPathNodes.front() and lastPathNodes.front() to be
    // falseVacuum.FieldConfiguration(), and returnPathNodes.back()
    // and lastPathNodes.back() to be trueVacuum.FieldConfiguration(),
    // assuming that every path given is between the vacua. It also resets
    // numberOfWorseningsSoFar.
    virtual void SetNodesForInitialPath( PotentialMinimum const& falseVacuum,
                                         PotentialMinimum const& trueVacuum );

    // This is an empty hook that can be over-ridden to account for the
    // bubble profile from the last path.
    virtual void AccountForBubbleProfileAroundNode( size_t nodeIndex,
                                   BubbleProfile const& bubbleFromLastPath ) {}
  };




  // This returns false if the last call of TryToImprovePath set nodesConverged
  // to be true because the nodes did not move very far. Otherwise, it takes
  // the bounce action from bubbleFromLastPath and updates
  // numberOfAllowedWorsenings based on whether it was an improvement on the
  // previous path. Then it returns false if too many paths have been tried
  // which just made the action bigger, true otherwise.
  inline bool MinuitOnPotentialPerpendicularToPath::PathCanBeImproved(
                                      BubbleProfile const& bubbleFromLastPath )
  {
    if( nodesConverged )
    {
      return false;
    }
    if( bubbleFromLastPath.BounceAction() >= bounceBeforeLastPath )
    {
      ++numberOfWorseningsSoFar;
    }
    bounceBeforeLastPath = bubbleFromLastPath.BounceAction();
    return ( numberOfAllowedWorsenings > numberOfWorseningsSoFar );
  }

  // This sets up nodeDisplacements and lastPathNodes to have the correct
  // sizes, sets returnPathNodes.front() and lastPathNodes.front() to be
  // falseVacuum.FieldConfiguration(), and returnPathNodes.back()
  // and lastPathNodes.back() to be trueVacuum.FieldConfiguration(),
  // assuming that every path given is between the vacua. It also resets
  // numberOfWorseningsSoFar.
  inline void MinuitOnPotentialPerpendicularToPath::SetNodesForInitialPath(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    nodeDisplacements = lastPathNodes
    = std::vector< Eigen::VectorXd >( ( numberOfVaryingNodes + 2 ),
                                     Eigen::VectorXd::Zero( numberOfFields ) );
    std::vector< double > const&
    falseConfiguration( falseVacuum.FieldConfiguration() );
    std::vector< double > const&
    trueConfiguration( trueVacuum.FieldConfiguration() );
    returnPathNodes.front() = falseConfiguration;
    returnPathNodes.back() = trueConfiguration;
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      lastPathNodes.front()( fieldIndex ) = falseConfiguration[ fieldIndex ];
      lastPathNodes.back()( fieldIndex ) = trueConfiguration[ fieldIndex ];
    }
    // We reset the number of allowed worsening steps here.
    numberOfWorseningsSoFar = 0;
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINUITONPOTENTIALPERPENDICULARTOPATH_HPP_ */
