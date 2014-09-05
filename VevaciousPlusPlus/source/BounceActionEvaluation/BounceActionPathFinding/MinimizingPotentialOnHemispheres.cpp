/*
 * MinimizingPotentialOnHemispheres.cpp
 *
 *  Created on: Sep 2, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinimizingPotentialOnHemispheres.hpp"

namespace VevaciousPlusPlus
{

  MinimizingPotentialOnHemispheres::MinimizingPotentialOnHemispheres(
                                    PotentialFunction const& potentialFunction,
                                               MinuitBetweenPaths* pathRefiner,
                                             size_t const minimumNumberOfNodes,
                                              size_t const movesPerImprovement,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
    MinimizingPotentialOnHypersurfaces( potentialFunction,
                                        pathRefiner,
                                        movesPerImprovement,
                                        minuitStrategy,
                                        minuitToleranceFraction ),
    minimumNumberOfNodes( minimumNumberOfNodes ),
    stepSize( NAN ),
    currentAngleScaling( NAN )
  {
    // This constructor is just an initialization list.
  }

  MinimizingPotentialOnHemispheres::~MinimizingPotentialOnHemispheres()
  {
    // This does nothing.
  }


  // This takes the parameterization as the components of a vector
  // perpendicular to axis 0 which is then added to a unit vector along axis
  // 0. This then could be scaled to give a unit vector within the hemisphere
  // of positive field 0, but first the angle it makes with axis 0 is scaled
  // so that the next node will be closer to the true vacuum than the current
  // last node before the true vacuum. This vector is then scaled to be of
  // length stepSize.
  Eigen::VectorXd MinimizingPotentialOnHemispheres::StepVector(
                      std::vector< double > const& nodeParameterization ) const
  {
    double tangentSquared( 0.0 );
    for( std::vector< double >::const_iterator
         nodeParameter( nodeParameterization.begin() );
         nodeParameter < nodeParameterization.end();
         ++nodeParameter )
    {
      tangentSquared += ( (*nodeParameter) * (*nodeParameter) );
    }
    double const parameterizationRadius( sqrt( tangentSquared ) );
    double const
    scaledAngle( currentAngleScaling * atan( parameterizationRadius ) );
    double const radialScaling( ( stepSize * sin( scaledAngle ) )
                                / parameterizationRadius );
    Eigen::VectorXd untransformedNode( numberOfFields );
    untransformedNode( 0 ) = ( stepSize * cos( scaledAngle ) );
    for( size_t fieldIndex( 1 );
        fieldIndex < numberOfFields;
        ++fieldIndex )
    {
      untransformedNode( fieldIndex )
      = ( radialScaling * nodeParameterization[ fieldIndex - 1 ] );
    }
    return ( reflectionMatrix * untransformedNode );
  }

} /* namespace VevaciousPlusPlus */
