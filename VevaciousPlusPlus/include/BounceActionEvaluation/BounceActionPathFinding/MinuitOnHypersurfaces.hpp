/*
 * MinuitOnHypersurfaces.hpp
 *
 *  Created on: Sep 2, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITONHYPERSURFACES_HPP_
#define MINUITONHYPERSURFACES_HPP_

#include "CommonIncludes.hpp"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "Eigen/Dense"
#include "MinuitPathFinder.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include "BounceActionEvaluation/PathParameterization/LinearSplineThroughNodes.hpp"
#include "MinuitWrappersAndHelpers/MinuitHypersphereBoundAlternative.hpp"

namespace VevaciousPlusPlus
{

  class MinuitOnHypersurfaces : public MinuitPathFinder
  {
  public:
    // This puts the elements of eigenVector into stlVector.
    static void ConvertEigenToStl( Eigen::VectorXd const& eigenVector,
                                   std::vector< double >& stlVector );

    // This puts the elements of stlVector into eigenVector.
    static void ConvertStlToEigen( std::vector< double > const& stlVector,
                                   Eigen::VectorXd& eigenVector );


    MinuitOnHypersurfaces( size_t const numberOfPathSegments,
                           unsigned int const minuitStrategy = 1,
                           double const minuitToleranceFraction = 0.5 );
    virtual ~MinuitOnHypersurfaces();


    // It may seem unwise to have this object call Minuit on itself, but really
    // it's just a handy way of keeping the minimization function within the
    // class that ends up finding its minimum. In this case,
    // nodeParameterization is just the parameterization of the node. The
    // default takes nodeParameterization as a vector in the hyperplane
    // perpendicular to field 0, applies reflectionMatrix to it, adds the
    // result to currentParallelComponent, and returns the potential function
    // for that field configuration, plus the fourth power of the Euclidean
    // length of nodeParameterization as a penalty to stop the path straying
    // too far from the zero parameterization. The penalty is an attempt to
    // stop jumping between paths to degenerate vacua.
    virtual double
    operator()( std::vector< double > const& nodeParameterization ) const;

    // This sets the potential function and vacua to be those given (also
    // noting the maximum allowed scale as the maximum Euclidean length for
    // field configurations, and also the number of fields which vary for the
    // potential), and resets the nodes to describe a straight path between the
    // new vacua, as well as setting pathTemperature and currentMinuitTolerance
    // appropriately.
    virtual void SetPotentialAndVacuaAndTemperature(
                                    PotentialFunction const& potentialFunction,
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                          double const pathTemperature = 0.0 );


  protected:
    PotentialFunction const* potentialFunction;
    double potentialAtOrigin;
    double maximumFieldVectorLengthSquared;
    size_t numberOfFields;
    size_t const numberOfVaryingNodes;
    double segmentAuxiliaryLength;
    std::vector< std::vector< double > > returnPathNodes;
    Eigen::VectorXd currentParallelComponent;
    // This is the vector between the current 2 reference nodes, scaled to be
    // appropriate for the step from the false-side node under a
    // parameterization that is just a set of zeroes.
    Eigen::VectorXd currentHyperplaneOrigin;
    // This is the vector to which the transformation of the vector in the
    // parameterization hyperplane [the hyperplane perpendicular to field 0]
    // gets added to yield the field configuration for the Minuit
    // parameterization. It has to be set by TryToImprovePath before operator()
    // is called (unless operator() has been over-ridden appropriately).
    Eigen::MatrixXd reflectionMatrix;
    std::vector< double > minuitInitialSteps;
    std::vector< double > const nodeZeroParameterization;
    Eigen::VectorXd minuitResultAsUntransformedVector;


    // This caps the field configuration so that the square of its Euclidean
    // length is not longer than maximumFieldVectorLengthSquared, then
    // evaluates the difference between potentialFunction->operator() on the
    // capped field configuration at pathTemperature and potentialAtOrigin, and
    // returns that difference, plus an extra amount if the field configuration
    // had to be capped, where the extra amount is the square of the amount that
    // the configuration's length squared exceeded
    // maximumFieldVectorLengthSquared.
    double PotentialValue( std::vector< double > fieldConfiguration ) const;

    // This sets up reflectionMatrix to be the Householder reflection matrix
    // which reflects the axis of field 0 to be parallel to
    // currentParallelComponent.
    void SetUpHouseholderReflection();

    // This takes the numberOfFields-1-dimensional vector and prepends a 0 to
    // make an numberOfFields-dimensional Eigen::VectorXd.
    Eigen::VectorXd UntransformedNode(
                     std::vector< double > const& nodeParameterization ) const;

    // This should set up returnPathNodes for a new pair of vacua. By default,
    // it just sets returnPathNodes.front() to be
    // falseVacuum.FieldConfiguration() and returnPathNodes.back() to be
    // trueVacuum.FieldConfiguration(), but this might need to be over-ridden.
    virtual void SetNodesForInitialPath( PotentialMinimum const& falseVacuum,
                                         PotentialMinimum const& trueVacuum );

    // This sets the components of currentParallelComponent to be the elements
    // of endNode minus the values those elements have in startNode.
    void SetParallelVector( std::vector< double > const& startNode,
                            std::vector< double > const& endNode );

    // This sets the components of currentParallelComponent to be the elements
    // of endNode minus the values those elements have in startNode.
    void SetParallelVector( Eigen::VectorXd const& startNode,
                            Eigen::VectorXd const& endNode );

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
    // into a node vector transformed by reflectionMatrix and puts that into
    // displacementVector.
    virtual Eigen::VectorXd RunMigradAndReturnDisplacement();

    // This creates and runs a Minuit2 MnMigrad object and converts the result
    // into a node vector (transformed by reflectionMatrix and added to
    // currentHyperplaneOrigin) and puts that into resultVector.
    virtual void
    RunMigradAndPutTransformedResultIn( std::vector< double >& resultVector );
  };





  // This puts the elements of eigenVector into stlVector.
  inline void
  MinuitOnHypersurfaces::ConvertEigenToStl( Eigen::VectorXd const& eigenVector,
                                            std::vector< double >& stlVector )
  {
    size_t const numberOfFields( eigenVector.rows() );
    stlVector.resize( numberOfFields );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      stlVector[ fieldIndex ] = eigenVector( fieldIndex );
    }
  }

  // This puts the elements of stlVector into eigenVector.
  inline void MinuitOnHypersurfaces::ConvertStlToEigen(
                                        std::vector< double > const& stlVector,
                                                 Eigen::VectorXd& eigenVector )
  {
    size_t const numberOfFields( stlVector.size() );
    eigenVector.resize( numberOfFields );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      eigenVector( fieldIndex ) = stlVector[ fieldIndex ];
    }
  }

  // It may seem unwise to have this object call Minuit on itself, but really
  // it's just a handy way of keeping the minimization function within the
  // class that ends up finding its minimum. In this case,
  // nodeParameterization is just the parameterization of the node. The
  // default takes nodeParameterization as a vector in the hyperplane
  // perpendicular to field 0, applies reflectionMatrix to it, adds the
  // result to currentParallelComponent, and returns the potential function
  // for that field configuration, plus the fourth power of the Euclidean
  // length of nodeParameterization as a penalty to stop the path straying
  // too far from the zero parameterization. The penalty is an attempt to
  // stop jumping between paths to degenerate vacua.
  inline double MinuitOnHypersurfaces::operator()(
                      std::vector< double > const& nodeParameterization ) const
  {
    Eigen::VectorXd const transformedNode( reflectionMatrix
                                 * UntransformedNode( nodeParameterization ) );
    std::vector< double > fieldConfiguration( numberOfFields );
    double parameterizationLengthSquared( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      fieldConfiguration[ fieldIndex ] = ( transformedNode( fieldIndex )
                                     + currentHyperplaneOrigin( fieldIndex ) );
      parameterizationLengthSquared
      += ( transformedNode( fieldIndex ) * transformedNode( fieldIndex ) );
    }
    return ( ( parameterizationLengthSquared * parameterizationLengthSquared )
             + PotentialValue( fieldConfiguration ) );
  }

  // This sets the potential function and vacua to be those given (also
  // noting the maximum allowed scale as the maximum Euclidean length for
  // field configurations, and also the number of fields which vary for the
  // potential), and resets the nodes to describe a straight path between the
  // new vacua, as well as setting pathTemperature and currentMinuitTolerance
  // appropriately.
  inline void MinuitOnHypersurfaces::SetPotentialAndVacuaAndTemperature(
                                    PotentialFunction const& potentialFunction,
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                                 double const pathTemperature )
  {
    this->potentialFunction = &potentialFunction;
    numberOfFields = potentialFunction.NumberOfFieldVariables();
    double const
    maximumFieldVectorLength( potentialFunction.GetLagrangianParameterManager(
                                                  ).MaximumEvaluationScale() );
    maximumFieldVectorLengthSquared
    = ( maximumFieldVectorLength * maximumFieldVectorLength );
    this->pathTemperature = pathTemperature;
    potentialAtOrigin
    = potentialFunction( potentialFunction.FieldValuesOrigin(),
                         pathTemperature );
    SetNodesForInitialPath( falseVacuum,
                            trueVacuum );
    SetCurrentMinuitTolerance( falseVacuum,
                               trueVacuum );
  }

  // This caps the field configuration so that the square of its Euclidean
  // length is not longer than maximumFieldVectorLengthSquared, then
  // evaluates the difference between potentialFunction->operator() on the
  // capped field configuration at pathTemperature and potentialAtOrigin, and
  // returns that difference, plus an extra amount if the field configuration
  // had to be capped, where the extra amount is the square of the amount that
  // the configuration's length squared exceeded
  // maximumFieldVectorLengthSquared.
  inline double MinuitOnHypersurfaces::PotentialValue(
                               std::vector< double > fieldConfiguration ) const
  {
    double const squaredLengthBeyondCap(
                          MinuitHypersphereBoundAlternative::CapVariableVector(
                                                            fieldConfiguration,
                                           maximumFieldVectorLengthSquared ) );
    return ( ( squaredLengthBeyondCap * squaredLengthBeyondCap )
             + (*potentialFunction)( fieldConfiguration,
                                     pathTemperature )
             - potentialAtOrigin );
  }

  // This takes the numberOfFields-1-dimensional vector and prepends a 0 to
  // make an numberOfFields-dimensional Eigen::VectorXd.
  inline Eigen::VectorXd MinuitOnHypersurfaces::UntransformedNode(
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
  // just sets returnPathNodes.front() to be falseVacuum.FieldConfiguration()
  // and returnPathNodes.back() to be trueVacuum.FieldConfiguration(), but this
  // might need to be over-ridden.
  inline void MinuitOnHypersurfaces::SetNodesForInitialPath(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    returnPathNodes.front() = falseVacuum.FieldConfiguration();
    returnPathNodes.back() = trueVacuum.FieldConfiguration();
  }

  // This sets the components of currentParallelComponent to be the elements
  // of endNode minus the values those elements have in startNode.
  inline void MinuitOnHypersurfaces::SetParallelVector(
                                        std::vector< double > const& startNode,
                                         std::vector< double > const& endNode )
  {
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      currentParallelComponent( fieldIndex ) = ( endNode[ fieldIndex ]
                                                 - startNode[ fieldIndex ] );
    }
  }

  // This sets the components of currentParallelComponent to be the elements
  // of endNode minus the values those elements have in startNode.
  inline void MinuitOnHypersurfaces::SetParallelVector(
                                              Eigen::VectorXd const& startNode,
                                               Eigen::VectorXd const& endNode )
  {
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      currentParallelComponent( fieldIndex ) = ( endNode( fieldIndex )
                                                 - startNode( fieldIndex ) );
    }
  }

  // This just sets the elements of minuitInitialSteps to be fractionOfLength
  // times the Euclidean length of currentParallelComponent.
  inline void
  MinuitOnHypersurfaces::SetCurrentMinuitSteps( double const fractionOfLength )
  {
    double lengthSquared( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      lengthSquared += ( currentParallelComponent( fieldIndex )
                         * currentParallelComponent( fieldIndex ) );
    }
    minuitInitialSteps.assign( ( numberOfFields - 1 ),
                               ( fractionOfLength * sqrt( lengthSquared ) ) );
  }

  // This creates and runs a Minuit2 MnMigrad object and converts the result
  // into a node vector (transformed by reflectionMatrix and added to
  // currentHyperplaneOrigin) and puts that into resultVector.
  inline void MinuitOnHypersurfaces::RunMigradAndPutTransformedResultIn(
                                          std::vector< double >& resultVector )
  {
    Eigen::VectorXd const
    displacementVector( RunMigradAndReturnDisplacement() );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      resultVector[ fieldIndex ] = ( currentHyperplaneOrigin( fieldIndex )
                                     + displacementVector( fieldIndex ) );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINUITONHYPERSURFACES_HPP_ */
