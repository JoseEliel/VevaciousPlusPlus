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
    currentHemisphereNode( NULL ),
    currentAngleScaling( NAN ),
    currentDistanceToTrueVacuum( NAN )
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
    std::cout << " } }";
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
      UpdateNode( pathNodes[ pathNodes.size() - 2 ],
                  *currentHemisphereNode,
                  StepVector( minuitResult.VariableValues() ) );

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
      std::cout << " }";
      std::cout << std::endl;/**/

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

} /* namespace VevaciousPlusPlus */
