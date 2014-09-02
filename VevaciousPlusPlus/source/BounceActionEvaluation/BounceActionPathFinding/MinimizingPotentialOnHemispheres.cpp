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
                                             size_t const minimumNumberOfNodes,
                                              size_t const movesPerImprovement,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
    MinimizingPotentialOnHypersurfaces( potentialFunction,
                                        movesPerImprovement,
                                        minuitStrategy,
                                        minuitToleranceFraction ),
    minimumNumberOfNodes( minimumNumberOfNodes ),
    currentAngleScaling( NAN )
  {
    // This constructor is just an initialization list.
  }

  MinimizingPotentialOnHemispheres::~MinimizingPotentialOnHemispheres()
  {
    // This does nothing.
  }


  // This allows Minuit2 to adjust the full path a set number of times to try
  // to minimize the sum of potentials at a set of nodes or bounce action
  // along the adjusted path, and then sets the path.
  TunnelPath const* MinimizingPotentialOnHemispheres::ImprovePath()
  {

  }

  // It may seem unwise to have this object call Minuit on itself, but really
  // it's just a handy way of keeping the minimization function within the
  // class that ends up finding its minimum. In this case,
  // nodeParameterization is just the parameterization of the node.
  double MinimizingPotentialOnHemispheres::operator()(
                      std::vector< double > const& nodeParameterization ) const
  {
    // The parameterization is taken as the components of a vector
    // perpendicular to axis 0 which is then added to a unit vector along axis
    // 0. This then could be scaled to give a unit vector within the hemisphere
    // of positive field 0, but first the angle it makes with axis 0 is scaled
    // so that the next node will be closer to the true vacuum than the current
    // last node before the true vacuum.
    double tangentSquared( 0.0 );
    for( std::vector< double >::const_iterator
         nodeParameter( nodeParameterization.begin() );
         nodeParameter < nodeParameterization.end();
         ++nodeParameter )
    {
      tangentSquared += ( (*nodeParameter) * (*nodeParameter) );
    }
    double const
    scaledAngle( currentAngleScaling * atan( sqrt( tangentSquared ) ) );

    Eigen::VectorXd const
    untransformedNode( UntransformedNode( nodeParameterization ) );
  }

  // This takes the parameterization as the components of a vector
  // perpendicular to axis 0 which is then added to a unit vector along axis
  // 0. This then could be scaled to give a unit vector within the hemisphere
  // of positive field 0, but first the angle it makes with axis 0 is scaled
  // so that the next node will be closer to the true vacuum than the current
  // last node before the true vacuum. This vector is then scaled to be of
  // length stepSize.
  std::vector< double > MinimizingPotentialOnHemispheres::StepVector(
      std::vector< double > const& nodeParameterization ) const
  {

  }

} /* namespace VevaciousPlusPlus */
