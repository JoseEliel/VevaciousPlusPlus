/*
 * MinimizingPotentialOnHemispheres.hpp
 *
 *  Created on: Sep 2, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINIMIZINGPOTENTIALONHEMISPHERES_HPP_
#define MINIMIZINGPOTENTIALONHEMISPHERES_HPP_

#include <limits>
#include "CommonIncludes.hpp"
#include "boost/math/constants/constants.hpp"
#include "MinimizingPotentialOnHypersurfaces.hpp"

namespace VevaciousPlusPlus
{

  class MinimizingPotentialOnHemispheres :
                                      public MinimizingPotentialOnHypersurfaces
  {
  public:
    MinimizingPotentialOnHemispheres(
                                    PotentialFunction const& potentialFunction,
                                      MinuitBetweenPaths* pathRefiner,
                                      size_t const minimumNumberOfNodes,
                                      size_t const movesPerImprovement = 100,
                                      unsigned int const minuitStrategy = 1,
                                  double const minuitToleranceFraction = 0.5 );
    virtual ~MinimizingPotentialOnHemispheres();

    // It may seem unwise to have this object call Minuit on itself, but really
    // it's just a handy way of keeping the minimization function within the
    // class that ends up finding its minimum. In this case,
    // nodeParameterization is just the parameterization of the node.
    virtual double
    operator()( std::vector< double > const& nodeParameterization ) const;


  protected:
    static double const pathForceAlignmentForZeroForce;

    size_t const minimumNumberOfNodes;
    double stepSize;
    std::vector< double > const* currentHemisphereNode;
    double currentAngleScaling;
    double currentDistanceToTrueVacuum;
    bool headingToSaddlePoint;
    double lastPotential;


    // This sets up the nodes by minimizing the potential on a hypersector of a
    // hypersphere of radius equal to ( 0.5 * stepSize ) from the last set node
    // beginning with the false vacuum, facing the true vacuum, and then makes
    // a step of length stepSize through this point on the hemisphere from the
    // node to create a new node, which is then used as the center of the
    // "hemisphere" for the next node. This continues until the direct distance
    // to the true vacuum is less than stepSize, or until Minuit2 has exceeded
    // movesPerImprovement function calls in total with this call of
    // ImproveNodes(), and nodesConverged is set appropriately.
    virtual void ImproveNodes();

    // This sets pathNodes to just be the false vacuum node and the true vacuum
    // node, and also sets stepSize to be the distance between the vacua
    // divided by ( minimumNumberOfNodes + 1 ).
    virtual void SetNodesForInitialPath( PotentialMinimum const& falseVacuum,
                                         PotentialMinimum const& trueVacuum );

    // This takes the parameterization as the components of a vector
    // perpendicular to axis 0 which is then added to a unit vector along axis
    // 0. This then could be scaled to give a unit vector within the hemisphere
    // of positive field 0, but first the angle it makes with axis 0 is scaled
    // so that the next node will be closer to the true vacuum than the current
    // last node before the true vacuum. This vector is then scaled to be of
    // length stepSize.
    Eigen::VectorXd
    StepVector( std::vector< double > const& nodeParameterization ) const;

    // This sets currentAngleScaling to be the ratio of the maximum angle
    // relative to currentParallelComponent that a step can take is asymptotic
    // to pi/2 (so that every step must put its node closer to the true vacuum
    // than the last node).
    virtual void SetUpCurrentAngleScaling()
    { currentAngleScaling = ( acos( stepSize / currentDistanceToTrueVacuum )
                              / ( boost::math::double_constants::half_pi ) ); }

    // This sets currentNode to be the vector sum of centerNode with stepVector
    // and also updates currentDistanceToTrueVacuum.
    void UpdateNode( std::vector< double >& currentNode,
                     std::vector< double > const& centerNode,
                     Eigen::VectorXd const& stepVector );

    // This returns -(F.S)^2/(F.F), where F is the gradient of the
    // potential at the pointOnHemisphere and S is stepVector. Essentially this
    // is returning -1/(S.S) times the square of the cosine of the angle
    // between the force and the path, and since stepVector is a fixed length,
    // minimizing -(F.S)^2/(F.F) is equivalent to maximizing the alignment of
    // the force with the path.
    double NegativeOfForcePathAlignment(
                                std::vector< double > const& pointOnHemisphere,
                                     Eigen::VectorXd const& stepVector ) const;
  };




  // It may seem unwise to have this object call Minuit on itself, but really
  // it's just a handy way of keeping the minimization function within the
  // class that ends up finding its minimum. In this case,
  // nodeParameterization is just the parameterization of the node.
  inline double MinimizingPotentialOnHemispheres::operator()(
                      std::vector< double > const& nodeParameterization ) const
  {
    Eigen::VectorXd const stepVector( StepVector( nodeParameterization ) );
    std::vector< double > pointOnHemisphere( *currentHemisphereNode );
    // Now pointOnHemisphere is the last node on the path before the true
    // vacuum.
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      pointOnHemisphere[ fieldIndex ] += ( 0.5 * stepVector( fieldIndex ) );
    }
    if( headingToSaddlePoint )
    {
      return NegativeOfForcePathAlignment( pointOnHemisphere,
                                           stepVector );
    }
    else
    {
      return potentialFunction( pointOnHemisphere,
                                pathTemperature );
    }
  }

  // This sets pathNodes to just be the false vacuum node and the true vacuum
  // node, and also sets stepSize to be the distance between the vacua
  // divided by ( minimumNumberOfNodes + 1 ).
  inline void MinimizingPotentialOnHemispheres::SetNodesForInitialPath(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    pathNodes.resize( 2 );
    pathNodes.front() = falseVacuum.FieldConfiguration();
    pathNodes.back() = trueVacuum.FieldConfiguration();
    currentDistanceToTrueVacuum
    = sqrt( trueVacuum.SquareDistanceTo( falseVacuum ) );
    stepSize = ( currentDistanceToTrueVacuum
                 / static_cast< double >( minimumNumberOfNodes + 1 ) );
    lastPotential = potentialFunction( falseVacuum.VariableValues() );
    headingToSaddlePoint = true;
    nodesConverged = false;
  }

  // This sets currentNode to be the vector sum of centerNode with stepVector
  // and also updates currentDistanceToTrueVacuum.
  inline void MinimizingPotentialOnHemispheres::UpdateNode(
                                            std::vector< double >& currentNode,
                                       std::vector< double > const& centerNode,
                                            Eigen::VectorXd const& stepVector )
  {
    double currentSquareDistanceToTrueVacuum( 0.0 );
    std::vector< double > const& trueVacuum( pathNodes.back() );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      currentNode[ fieldIndex ] = ( centerNode[ fieldIndex ]
                                    + stepVector( fieldIndex ) );
      double const differenceFromTrueVacuum( trueVacuum[ fieldIndex ]
                                             - currentNode[ fieldIndex ] );
      currentSquareDistanceToTrueVacuum += ( differenceFromTrueVacuum
                                             * differenceFromTrueVacuum );
    }
    currentDistanceToTrueVacuum = sqrt( currentSquareDistanceToTrueVacuum );
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINIMIZINGPOTENTIALONHEMISPHERES_HPP_ */
