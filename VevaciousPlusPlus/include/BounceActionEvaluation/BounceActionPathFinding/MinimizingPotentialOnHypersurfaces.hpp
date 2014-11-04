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

namespace VevaciousPlusPlus
{

  class MinimizingPotentialOnHypersurfaces : public MinuitPathFinder
  {
  public:
    MinimizingPotentialOnHypersurfaces(
                                    PotentialFunction const& potentialFunction,
                                        size_t const numberOfPathSegments,
                                        unsigned int const minuitStrategy = 1,
                                  double const minuitToleranceFraction = 0.5 );
    virtual ~MinimizingPotentialOnHypersurfaces();


    // It may seem unwise to have this object call Minuit on itself, but really
    // it's just a handy way of keeping the minimization function within the
    // class that ends up finding its minimum. In this case,
    // nodeParameterization is just the parameterization of the node. The
    // default takes nodeParameterization as a vector in the hyperplane
    // perpendicular to field 0, applies reflectionMatrix to it, adds the
    // result to currentParallelComponent, and returns the potential function
    // for that field configuration.
    virtual double
    operator()( std::vector< double > const& nodeParameterization ) const;

    // This sets the vacua to be those given, and resets the nodes to describe
    // a straight path between the new vacua, as well as setting
    // pathTemperature and currentMinuitTolerance appropriately.
    virtual void SetVacuaAndTemperature( PotentialMinimum const& falseVacuum,
                                         PotentialMinimum const& trueVacuum,
                                         double const pathTemperature = 0.0 );


  protected:
    PotentialFunction const& potentialFunction;
    size_t const numberOfFields;
    size_t const numberOfVaryingNodes;
    std::vector< std::vector< double > > pathNodes;
    std::vector< double > currentParallelComponent;
    // This is the vector between the current 2 reference nodes, scaled to be
    // appropriate for the step from the false-side node under a
    // parameterization that is just a set of zeroes.
    std::vector< double > currentHyperplaneOrigin;
    // This is the vector to which the transformation of the vector in the
    // parameterization hyperplane [the hyperplane perpendicular to field 0]
    // gets added to yield the field configuration for the Minuit
    // parameterization. It has to be set by TryToImprovePath before operator()
    // is called (unless operator() has been over-ridden appropriately).
    Eigen::MatrixXd reflectionMatrix;
    std::vector< double > minuitInitialSteps;
    std::vector< double > const nodeZeroParameterization;
    Eigen::VectorXd minuitResultAsUntransformedVector;


    // This sets up reflectionMatrix to be the Householder reflection matrix
    // which reflects the axis of field 0 to be parallel to
    // currentParallelComponent.
    void SetUpHouseholderReflection();

    // This takes the numberOfFields-1-dimensional vector and prepends a 0 to
    // make an numberOfFields-dimensional Eigen::VectorXd.
    Eigen::VectorXd UntransformedNode(
                     std::vector< double > const& nodeParameterization ) const;

    // This should set up pathNodes for a new pair of vacua. By default, it
    // just sets pathNodes.front() to be falseVacuum.FieldConfiguration() and
    // pathNodes.back() to be trueVacuum.FieldConfiguration(), but this might
    // need to be over-ridden.
    virtual void SetNodesForInitialPath( PotentialMinimum const& falseVacuum,
                                         PotentialMinimum const& trueVacuum );

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

    // This just sets the elements of minuitInitialSteps to be fractionOfLength
    // times the Euclidean length of currentParallelComponent.
    virtual void SetCurrentMinuitSteps( double const fractionOfLength );

    // This creates and runs a Minuit2 MnMigrad object and converts the result
    // into a node vector (transformed by reflectionMatrix and added to
    // currentHyperplaneOrigin) and puts that into resultVector.
    virtual void
    RunMigradAndPutTransformedResultIn( std::vector< double >& resultVector );
  };





  // It may seem unwise to have this object call Minuit on itself, but really
  // it's just a handy way of keeping the minimization function within the
  // class that ends up finding its minimum. In this case,
  // nodeParameterization is just the parameterization of the node. The
  // default takes nodeParameterization as a vector in the hyperplane
  // perpendicular to field 0, applies reflectionMatrix to it, adds the
  // result to currentParallelComponent, and returns the potential function
  // for that field configuration.
  inline double MinimizingPotentialOnHypersurfaces::operator()(
                      std::vector< double > const& nodeParameterization ) const
  {
    Eigen::VectorXd const transformedNode( reflectionMatrix
                                 * UntransformedNode( nodeParameterization ) );
    std::vector< double > fieldConfiguration( currentHyperplaneOrigin );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      fieldConfiguration[ fieldIndex ] += transformedNode( fieldIndex );
    }
    return potentialFunction( fieldConfiguration,
                              pathTemperature );
  }

  // This sets the vacua to be those given, and resets the nodes to describe a
  // straight path between the new vacua, as well as setting pathTemperature
  // and currentMinuitTolerance appropriately.
  inline void MinimizingPotentialOnHypersurfaces::SetVacuaAndTemperature(
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                                 double const pathTemperature )
  {
    this->pathTemperature = pathTemperature;
    SetNodesForInitialPath( falseVacuum,
                            trueVacuum );
    SetCurrentMinuitTolerance( falseVacuum,
                               trueVacuum );
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

  // This should set up pathNodes for a new pair of vacua. By default, it
  // just sets pathNodes.front() to be falseVacuum.FieldConfiguration() and
  // pathNodes.back() to be trueVacuum.FieldConfiguration(), but this might
  // need to be over-ridden.
  inline void MinimizingPotentialOnHypersurfaces::SetNodesForInitialPath(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    pathNodes.front() = falseVacuum.FieldConfiguration();
    pathNodes.back() = trueVacuum.FieldConfiguration();
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

  // This just sets the elements of minuitInitialSteps to be fractionOfLength
  // times the Euclidean length of currentParallelComponent.
  inline void MinimizingPotentialOnHypersurfaces::SetCurrentMinuitSteps(
                                                double const fractionOfLength )
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
                               ( fractionOfLength * sqrt( lengthSquared ) ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINIMIZINGPOTENTIALONHYPERSURFACES_HPP_ */
