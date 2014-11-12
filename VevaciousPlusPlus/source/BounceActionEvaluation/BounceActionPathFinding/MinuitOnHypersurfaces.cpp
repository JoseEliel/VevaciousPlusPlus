/*
 * MinuitOnHypersurfaces.cpp
 *
 *  Created on: Sep 2, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinuitOnHypersurfaces.hpp"

namespace VevaciousPlusPlus
{
  MinuitOnHypersurfaces::MinuitOnHypersurfaces(
                                    PotentialFunction const& potentialFunction,
                                             size_t const numberOfPathSegments,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
    MinuitPathFinder( minuitStrategy,
                      minuitToleranceFraction ),
    potentialFunction( potentialFunction ),
    potentialAtOrigin( -1.0 ),
    numberOfFields( potentialFunction.NumberOfFieldVariables() ),
    numberOfVaryingNodes( numberOfPathSegments - 1 ),
    segmentAuxiliaryLength( 1.0
                            / static_cast< double > ( numberOfPathSegments ) ),
    returnPathNodes( ( numberOfVaryingNodes + 2 ),
                     std::vector< double >( numberOfFields ) ),
    currentParallelComponent( Eigen::VectorXd::Zero( numberOfFields ) ),
    currentHyperplaneOrigin( currentParallelComponent ),
    reflectionMatrix( numberOfFields,
                      numberOfFields ),
    nodeZeroParameterization( ( numberOfFields - 1 ),
                              0.0 ),
    minuitResultAsUntransformedVector( currentParallelComponent )
  {
    // This constructor is just an initialization list.
  }

  MinuitOnHypersurfaces::~MinuitOnHypersurfaces()
  {
    // This does nothing.
  }


  // This sets up reflectionMatrix to be the Householder reflection matrix
  // which reflects the axis of field 0 to be parallel to
  // currentParallelComponent.
  void MinuitOnHypersurfaces::SetUpHouseholderReflection()
  {
    // First we check that targetVector doesn't already lie on the axis of the
    // field with index 0.
    bool alreadyParallel( true );
    for( size_t fieldIndex( 1 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( currentParallelComponent( fieldIndex ) != 0.0 )
      {
        alreadyParallel = false;
        break;
      }
    }
    if( alreadyParallel )
    {
      reflectionMatrix = Eigen::MatrixXd::Identity( numberOfFields,
                                                    numberOfFields );
      reflectionMatrix( 0,
                        0 ) = -1.0;
    }
    else
    {
      double targetLengthSquared( 0.0 );
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        targetLengthSquared += ( currentParallelComponent( fieldIndex )
                                 * currentParallelComponent( fieldIndex ) );
      }
      double const targetNormalization( 1.0 / sqrt( targetLengthSquared ) );
      double const minusInverseOfOneMinusDotProduct( 1.0 /
           ( ( currentParallelComponent( 0 ) * targetNormalization ) - 1.0 ) );
      reflectionMatrix = Eigen::MatrixXd::Zero( numberOfFields,
                                                numberOfFields );
      for( size_t rowIndex( 0 );
           rowIndex < numberOfFields;
           ++rowIndex )
      {
        double rowIndexPart( currentParallelComponent( rowIndex )
                             * targetNormalization );
        if( rowIndex == 0 )
        {
          rowIndexPart -= 1.0;
        }
        for( size_t columnIndex( rowIndex + 1 );
             columnIndex < numberOfFields;
             ++columnIndex )
        {
          double columnIndexPart( currentParallelComponent( columnIndex )
                                  * targetNormalization );
          if( columnIndex == 0 )
          {
            columnIndexPart -= 1.0;
          }
          reflectionMatrix( rowIndex,
                            columnIndex ) = ( rowIndexPart * columnIndexPart
                                          * minusInverseOfOneMinusDotProduct );
          reflectionMatrix( columnIndex,
                            rowIndex ) = ( rowIndexPart * columnIndexPart
                                          * minusInverseOfOneMinusDotProduct );
          // The reflection matrix is symmetric, and we set the diagonal
          // elements outside this loop.
        }
        reflectionMatrix( rowIndex,
                          rowIndex ) = ( 1.0 + ( rowIndexPart * rowIndexPart
                                        * minusInverseOfOneMinusDotProduct ) );
      }
    }
  }

  // This creates and runs a Minuit2 MnMigrad object and converts the result
  // into a node vector transformed by reflectionMatrix and puts that into
  // displacementVector.
  Eigen::VectorXd
  MinuitOnHypersurfaces::RunMigradAndReturnDisplacement()
  {
    ROOT::Minuit2::MnMigrad mnMigrad( *this,
                                      nodeZeroParameterization,
                                      minuitInitialSteps,
                                      minuitStrategy );
    ROOT::Minuit2::FunctionMinimum const minuitResult( mnMigrad( 0,
                                                    currentMinuitTolerance ) );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "minuitResult.Fval() = " << minuitResult.Fval()
    << ", (*this)( nodeZeroParameterization ) = "
    << (*this)( nodeZeroParameterization );
    std::cout << std::endl;/**/

    // We return a zero displacement if Minuit2 failed to minimize operator()
    // better than that.
    if( minuitResult.Fval() > (*this)( nodeZeroParameterization ) )
    {
      return Eigen::VectorXd::Zero( numberOfFields );
    }
    ROOT::Minuit2::MnUserParameters const&
    userParameters( minuitResult.UserParameters() );
    // We assume that minuitResultAsUntransformedVector( 0 ) was set to 0.0
    // in the constructor and never changes.
    for( size_t variableIndex( 1 );
         variableIndex < numberOfFields;
         ++variableIndex )
    {
      minuitResultAsUntransformedVector( variableIndex )
      = userParameters.Value( variableIndex - 1 );
    }
    return ( reflectionMatrix * minuitResultAsUntransformedVector );
  }

} /* namespace VevaciousPlusPlus */
