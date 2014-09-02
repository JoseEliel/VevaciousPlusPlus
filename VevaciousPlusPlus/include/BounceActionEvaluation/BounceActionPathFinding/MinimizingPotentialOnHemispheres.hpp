/*
 * MinimizingPotentialOnHemispheres.hpp
 *
 *  Created on: Sep 2, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINIMIZINGPOTENTIALONHEMISPHERES_HPP_
#define MINIMIZINGPOTENTIALONHEMISPHERES_HPP_

#include "CommonIncludes.hpp"
#include "MinimizingPotentialOnHypersurfaces.hpp"

namespace VevaciousPlusPlus
{

  class MinimizingPotentialOnHemispheres :
                                      public MinimizingPotentialOnHypersurfaces
  {
  public:
    MinimizingPotentialOnHemispheres(
                                    PotentialFunction const& potentialFunction,
                                      size_t const minimumNumberOfNodes,
                                      size_t const movesPerImprovement = 100,
                                      unsigned int const minuitStrategy = 1,
                                  double const minuitToleranceFraction = 0.5 );
    virtual ~MinimizingPotentialOnHemispheres();


    // This allows Minuit2 to adjust the full path a set number of times to try
    // to minimize the sum of potentials at a set of nodes or bounce action
    // along the adjusted path, and then sets the path.
    virtual TunnelPath const* ImprovePath();

    // It may seem unwise to have this object call Minuit on itself, but really
    // it's just a handy way of keeping the minimization function within the
    // class that ends up finding its minimum. In this case,
    // nodeParameterization is just the parameterization of the node.
    virtual double
    operator()( std::vector< double > const& nodeParameterization ) const;


  protected:
    size_t const minimumNumberOfNodes;
    double stepSize;
    double currentAngleScaling;


    // This sets pathNodes to just be the false vacuum node and the true vacuum
    // node, and also sets stepSize to be the distance between the vacua
    // divided by ( minimumNumberOfNodes + 1 )
    virtual void SetNodesForInitialPath( PotentialMinimum const& falseVacuum,
                                         PotentialMinimum const& trueVacuum );

    // This takes the parameterization as the components of a vector
    // perpendicular to axis 0 which is then added to a unit vector along axis
    // 0. This then could be scaled to give a unit vector within the hemisphere
    // of positive field 0, but first the angle it makes with axis 0 is scaled
    // so that the next node will be closer to the true vacuum than the current
    // last node before the true vacuum. This vector is then scaled to be of
    // length stepSize.
    std::vector< double >
    StepVector( std::vector< double > const& nodeParameterization ) const;
  };




  // This sets pathNodes to just be the false vacuum node and the true vacuum
  // node, and also sets stepSize to be the distance between the vacua
  // divided by ( minimumNumberOfNodes + 1 )
  inline void MinimizingPotentialOnHemispheres::SetNodesForInitialPath(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    pathNodes.resize( 2 );
    pathNodes.front() = falseVacuum.FieldConfiguration();
    pathNodes.back() = trueVacuum.FieldConfiguration();
    stepSize = ( sqrt( trueVacuum.SquareDistanceTo( falseVacuum ) )
                 / static_cast< double >( minimumNumberOfNodes + 1 ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINIMIZINGPOTENTIALONHEMISPHERES_HPP_ */
