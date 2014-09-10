/*
 * MinimizingPotentialOnHemispheres.cpp
 *
 *  Created on: Sep 2, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinimizingPotentialOnHemispheres.hpp"

namespace VevaciousPlusPlus
{
  double const
  MinimizingPotentialOnHemispheres::pathForceAlignmentForZeroForce(
                                     -(std::numeric_limits< double >::max()) );

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
    currentHemisphereNode( NULL ),
    currentAngleScaling( NAN ),
    currentDistanceToTrueVacuum( NAN ),
    headingToSaddlePoint( true ),
    lastPotential( NAN )
  {
    // This constructor is just an initialization list.
  }

  MinimizingPotentialOnHemispheres::~MinimizingPotentialOnHemispheres()
  {
    // This does nothing.
  }


  // This sets up the nodes by minimizing the potential on a hypersector of a
  // hypersphere of radius equal to ( 0.5 * stepSize ) from the last set node
  // beginning with the false vacuum, facing the true vacuum, and then makes
  // a step of length stepSize through this point on the hemisphere from the
  // node to create a new node, which is then used as the center of the
  // "hemisphere" for the next node. This continues until the direct distance
  // to the true vacuum is less than stepSize, or until Minuit2 has exceeded
  // movesPerImprovement function calls in total with this call of
  // ImproveNodes(), and nodesConverged is set appropriately.
  void MinimizingPotentialOnHemispheres::ImproveNodes()
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "MinimizingPotentialOnHemispheres::ImproveNodes() called."
    << " pathNodes = { { ";
    for( std::vector< std::vector< double > >::iterator
         pathNode( pathNodes.begin() );
         pathNode < pathNodes.end();
         ++pathNode )
    {
      if( pathNode > pathNodes.begin() )
      {
        std::cout << " }," << std::endl << "{ ";
      }
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        if( fieldIndex > 0 )
        {
          std::cout << ", ";
        }
        std::cout << (*pathNode)[ fieldIndex ];
      }
    }
    std::cout << " } }, lastPotential = " << lastPotential
    << ", headingToSaddlePoint = " << headingToSaddlePoint;
    std::cout << std::endl;/**/

    for( size_t numberOfMovesSoFar( 0 );
         numberOfMovesSoFar < movesPerImprovement;
         ++numberOfMovesSoFar )
    {
      pathNodes.push_back( pathNodes.back() );
      // Now the true vacuum has been copied into
      // pathNodes[ pathNodes.size() - 1 ], that means that we can overwrite
      // pathNodes[ pathNodes.size() - 2 ] (which currently is a copy of the
      // true vacuum) as the node just before the true vacuum, based on the
      // node just before it, which is pathNodes[ pathNodes.size() - 3 ].
      currentHemisphereNode = &(pathNodes[ pathNodes.size() - 3 ]);
      SetParallelVector( *currentHemisphereNode,
                         pathNodes.back() );
      SetUpHouseholderReflection();
      SetCurrentMinuitSteps();
      SetUpCurrentAngleScaling();

      ROOT::Minuit2::MnMigrad mnMigrad( *this,
                                 std::vector< double >( ( numberOfFields - 1 ),
                                                        0.0 ),
                                        minuitInitialSteps,
                                        minuitStrategy );
      MinuitMinimum minuitResult( ( numberOfFields - 1 ),
                                  mnMigrad( 0,
                                            currentMinuitTolerance ) );
      // We use 0 as the "maximum number of function calls" so that Minuit2
      // uses the default number.
      Eigen::VectorXd const
      stepVector( StepVector( minuitResult.VariableValues() ) );
      std::vector< double >& currentNode( pathNodes[ pathNodes.size() - 2 ] );
      UpdateNode( currentNode,
                  *currentHemisphereNode,
                  StepVector( minuitResult.VariableValues() ) );

      // Now we check to see if we have gone over the saddle point, and note if
      // we have so that we can just minimize the potential rather than the
      // alignment for generating the rest of the path.
      double const currentPotential( potentialFunction( currentNode ) );
      if( currentPotential < lastPotential )
      {
        headingToSaddlePoint = false;
      }
      lastPotential = currentPotential;

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "currentParallelComponent = { ";
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        if( fieldIndex > 0 )
        {
          std::cout << ", ";
        }
        std::cout << currentParallelComponent[ fieldIndex ];
      }
      std::cout << " }" << std::endl << "minuitResult = "
      << minuitResult.AsDebuggingString() << std::endl << "reflectionMatrix ="
      << std::endl << reflectionMatrix << std::endl
      << "pathNodes[ pathNodes.size() - 2 ] = { ";
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        if( fieldIndex > 0 )
        {
          std::cout << ", ";
        }
        std::cout << pathNodes[ pathNodes.size() - 2 ][ fieldIndex ];
      }
      std::cout << " }, currentDistanceToTrueVacuum = "
      << currentDistanceToTrueVacuum << std::endl;
      std:: cout << "*currentHemisphereNode = { ";
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        if( fieldIndex > 0 )
        {
          std::cout << ", ";
        }
        std::cout << (*currentHemisphereNode)[ fieldIndex ];
      }
      std::cout << " }" << std::endl;
      std::cout << "stepVector = { ";
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        if( fieldIndex > 0 )
        {
          std::cout << ", ";
        }
        std::cout << stepVector( fieldIndex );
      }
      std::cout << " }, headingToSaddlePoint = " << headingToSaddlePoint
      << ", potential of current node = " << lastPotential
      << ", potential of last node = "
      << potentialFunction( *currentHemisphereNode );
      std::cout << std::endl;
      std::vector< double > gradientPoint( *currentHemisphereNode );
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        gradientPoint += ( 0.5 * stepVector( fieldIndex ) );
      }
      std::vector< double > forceVector;
      potentialFunction.SetAsGradientAt( forceVector,
                                         gradientPoint );
      std::cout << "gradientPoint = { ";
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        if( fieldIndex > 0 )
        {
          std::cout << ", ";
        }
        std::cout << gradientPoint[ fieldIndex ];
      }
      std::cout << "}, forceVector = { ";
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        if( fieldIndex > 0 )
        {
          std::cout << ", ";
        }
        std::cout << forceVector[ fieldIndex ];
      }
      std::cout << " }";/**/

      if( currentDistanceToTrueVacuum <= stepSize )
      {
        nodesConverged = true;
        return;
        // If currentNode is within stepSize of the true vacuum, we note that
        // this stage of path improvement is over by setting nodesConverged to
        // true, and then leaving the function immediately (rather than just
        // breaking out of the loop, as nothing else needs to happen after the
        // loop).
      }
    }
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
    Eigen::VectorXd untransformedNode( numberOfFields );
    double tangentSquared( 0.0 );
    for( std::vector< double >::const_iterator
         nodeParameter( nodeParameterization.begin() );
         nodeParameter < nodeParameterization.end();
         ++nodeParameter )
    {
      tangentSquared += ( (*nodeParameter) * (*nodeParameter) );
    }
    if( tangentSquared > 0.0 )
    {
      double const parameterizationRadius( sqrt( tangentSquared ) );
      double const
      scaledAngle( currentAngleScaling * atan( parameterizationRadius ) );
      double const radialScaling( ( stepSize * sin( scaledAngle ) )
                                  / parameterizationRadius );
      untransformedNode( 0 ) = ( stepSize * cos( scaledAngle ) );
      for( size_t fieldIndex( 1 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        untransformedNode( fieldIndex )
        = ( radialScaling * nodeParameterization[ fieldIndex - 1 ] );
      }
    }
    else
    {
      untransformedNode( 0 ) = stepSize;
      for( size_t fieldIndex( 1 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        untransformedNode( fieldIndex ) = 0.0;
      }
    }

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "MinimizingPotentialOnHemispheres::StepVector("
    << " nodeParameterization = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < ( numberOfFields - 1 );
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << nodeParameterization[ fieldIndex ];
    }
    std::cout << " } ) finishing. tangentSquared = " << tangentSquared
    << ", currentAngleScaling = " << currentAngleScaling
    << ", stepSize = " << stepSize
    << ", untransformedNode =" << std::endl << untransformedNode
    << ", about to return ( reflectionMatrix * untransformedNode ) ="
    << std::endl << ( reflectionMatrix * untransformedNode );
    std::cout << std::endl;*/

    return ( reflectionMatrix * untransformedNode );
  }

  // This returns -(F.S)^2/(F.F), where F is the gradient of the
  // potential at the pointOnHemisphere and S is stepVector. Essentially this
  // is returning -1/(S.S) times the square of the cosine of the angle
  // between the force and the path, and since stepVector is a fixed length,
  // minimizing -(F.S)^2/(F.F) is equivalent to maximizing the alignment of
  // the force with the path.
  double MinimizingPotentialOnHemispheres::NegativeOfForcePathAlignment(
                                std::vector< double > const& pointOnHemisphere,
                                      Eigen::VectorXd const& stepVector ) const
  {
    std::vector< double > forceVector;
    potentialFunction.SetAsGradientAt( forceVector,
                                       pointOnHemisphere );
    double pathDotForce( 0.0 );
    double minusForceDotForce( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      pathDotForce += ( stepVector( fieldIndex ) * forceVector[ fieldIndex ] );
      minusForceDotForce
      -= ( forceVector[ fieldIndex ] * forceVector[ fieldIndex ] );
    }
    if( minusForceDotForce < 0.0 )
    {
      return ( ( pathDotForce * pathDotForce ) / minusForceDotForce );
    }
    else
    {
      return pathForceAlignmentForZeroForce;
    }
  }

} /* namespace VevaciousPlusPlus */
