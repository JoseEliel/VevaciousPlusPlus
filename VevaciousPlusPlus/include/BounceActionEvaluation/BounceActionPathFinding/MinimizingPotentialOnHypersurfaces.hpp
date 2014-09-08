/*
 * MinimizingPotentialOnHypersurfaces.hpp
 *
 *  Created on: Sep 2, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINIMIZINGPOTENTIALONHYPERSURFACES_HPP_
#define MINIMIZINGPOTENTIALONHYPERSURFACES_HPP_

#include "CommonIncludes.hpp"
#include "MinuitPathFinder.hpp"
#include "Eigen/Dense"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include "../PathParameterization/LinearSplineThroughNodes.hpp"
#include "MinuitBetweenPaths.hpp"

namespace VevaciousPlusPlus
{

  class MinimizingPotentialOnHypersurfaces : public MinuitPathFinder
  {
  public:
    MinimizingPotentialOnHypersurfaces(
                                    PotentialFunction const& potentialFunction,
                                        MinuitBetweenPaths* pathRefiner = NULL,
                                        size_t const movesPerImprovement = 100,
                                        unsigned int const minuitStrategy = 1,
                                  double const minuitToleranceFraction = 0.5 );
    virtual ~MinimizingPotentialOnHypersurfaces();


    // This sets the vacua to be those given, and resets the nodes to describe
    // a straight path between the new vacua, as well as setting
    // pathTemperature and currentMinuitTolerance appropriately.
    virtual TunnelPath const*
    SetInitialPath( PotentialMinimum const& falseVacuum,
                    PotentialMinimum const& trueVacuum,
                    double const pathTemperature = 0.0 );


    // This allows Minuit2 to adjust the full path a set number of times to try
    // to minimize the sum of potentials at a set of nodes or bounce action
    // along the adjusted path, and then sets the path.
    virtual TunnelPath const* ImprovePath();


  protected:
    PotentialFunction const& potentialFunction;
    size_t const numberOfFields;
    std::vector< std::vector< double > > pathNodes;
    std::vector< double > currentParallelComponent;
    // This is the vector between the current 2 reference nodes, scaled to be
    // appropriate for the step from the false-side node under a
    // parameterization that is just a set of zeroes.
    Eigen::MatrixXd reflectionMatrix;
    bool nodesConverged;
    MinuitBetweenPaths* pathRefiner;
    std::vector< double > minuitInitialSteps;


    // This should move the nodes individually towards whatever the derived
    // class considers to be the optimal tunneling path, by some amount, which
    // is up to the derived class, but it should also note if it is finished
    // moving the nodes around by setting nodesConverged to true.
    virtual void ImproveNodes() = 0;

    // This sets up reflectionMatrix to be the Householder reflection matrix
    // which reflects the axis of field 0 to be parallel to
    // currentParallelComponent.
    void SetUpHouseholderReflection();

    // This takes the numberOfFields-1-dimensional vector and prepends a 0 to
    // make an numberOfFields-dimensional Eigen::VectorXd.
    Eigen::VectorXd UntransformedNode(
                     std::vector< double > const& nodeParameterization ) const;

    // This should set up pathNodes for a new pair of vacua.
    virtual void SetNodesForInitialPath( PotentialMinimum const& falseVacuum,
                                      PotentialMinimum const& trueVacuum ) = 0;

    // This sets the components of currentParallelComponent to be the elements
    // of endNode minus the values those elements have in startNode.
    void SetParallelVector( std::vector< double > const& startNode,
                            std::vector< double > const& endNode );

    // This just sets the tolerance to be minuitToleranceFraction times the
    // difference in potential between the vacua, but could be over-ridden if
    // necessary.
    virtual void
    SetCurrentMinuitTolerance( PotentialMinimum const& falseVacuum,
                               PotentialMinimum const& trueVacuum )
    { currentMinuitTolerance = ( minuitToleranceFraction
          * ( falseVacuum.PotentialValue() - trueVacuum.PotentialValue() ) ); }

    // This just sets the elements of minuitInitialSteps to be 0.5 times the
    // Euclidean length of currentParallelComponent.
    virtual void SetCurrentMinuitSteps();
  };




  // This sets the vacua to be those given, and resets the nodes to describe a
  // straight path between the new vacua, as well as setting pathTemperature
  // and currentMinuitTolerance appropriately.
  inline TunnelPath const* MinimizingPotentialOnHypersurfaces::SetInitialPath(
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                                 double const pathTemperature )
  {
    this->pathTemperature = pathTemperature;
    SetNodesForInitialPath( falseVacuum,
                            trueVacuum );
    SetCurrentMinuitTolerance( falseVacuum,
                               trueVacuum );
    SetCurrentMinuitTolerance( falseVacuum,
                               trueVacuum );
    std::vector< std::vector< double > > straightPath( 2,
                                            falseVacuum.FieldConfiguration() );
    straightPath.back() = trueVacuum.FieldConfiguration();
    return new LinearSplineThroughNodes( straightPath,
                                         std::vector< double >( 0 ),
                                         pathTemperature );
  }

  // This takes the numberOfFields-1-dimensional vector and prepends a 0 to
  // make an numberOfFields-dimensional Eigen::VectorXd.
  inline Eigen::VectorXd MinimizingPotentialOnHypersurfaces::UntransformedNode(
                      std::vector< double > const& nodeParameterization ) const
  {
    Eigen::VectorXd untransformedNode( numberOfFields );
    untransformedNode( 0 ) = 0.0;
    for( size_t fieldIndex( 1 );
        fieldIndex < numberOfFields;
        ++fieldIndex )
    {
      untransformedNode( fieldIndex ) = nodeParameterization[ fieldIndex - 1 ];
    }
    return untransformedNode;
  }

  // This sets the components of currentParallelComponent to be the elements
  // of endNode minus the values those elements have in startNode.
  inline void MinimizingPotentialOnHypersurfaces::SetParallelVector(
                                        std::vector< double > const& startNode,
                                         std::vector< double > const& endNode )
  {
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      currentParallelComponent[ fieldIndex ] = ( endNode[ fieldIndex ]
                                                 - startNode[ fieldIndex ] );
    }
  }

  // This just sets the elements of minuitInitialSteps to be 0.5 times the
  // Euclidean length of currentParallelComponent.
  inline void MinimizingPotentialOnHypersurfaces::SetCurrentMinuitSteps()
  {
    double lengthSquared( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      lengthSquared += ( currentParallelComponent[ fieldIndex ]
                         * currentParallelComponent[ fieldIndex ] );
    }
    minuitInitialSteps.assign( ( numberOfFields - 1 ),
                               sqrt( lengthSquared ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINIMIZINGPOTENTIALONHYPERSURFACES_HPP_ */
