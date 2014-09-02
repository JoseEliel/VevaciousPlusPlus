/*
 * MinimizingPotentialOnHypersurfaces.cpp
 *
 *  Created on: Sep 2, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinimizingPotentialOnHypersurfaces.hpp"

namespace VevaciousPlusPlus
{

  MinimizingPotentialOnHypersurfaces::MinimizingPotentialOnHypersurfaces(
                                    PotentialFunction const& potentialFunction,
                                              size_t const movesPerImprovement,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
    MinuitPathFinder( movesPerImprovement,
                      minuitStrategy,
                      minuitToleranceFraction ),
    potentialFunction( potentialFunction ),
    numberOfFields( potentialFunction.NumberOfFieldVariables() ),
    pathNodes(),
    currentParallelComponent(),
    reflectionMatrix( numberOfFields,
                      numberOfFields ),
    nodesConverged( false )
  {
    // This constructor is just an initialization list.
  }

  MinimizingPotentialOnHypersurfaces::~MinimizingPotentialOnHypersurfaces()
  {
    // This does nothing.
  }



  // This allows Minuit2 to adjust the full path a set number of times to try
  // to minimize the sum of potentials at a set of nodes or bounce action
  // along the adjusted path, and then sets the path.
  TunnelPath const* MinimizingPotentialOnHypersurfaces::ImprovePath()
  {
    if( !nodesConverged )
    {
      ImproveNodes();
      // If this call of ImproveNodes() was the last one for these vacua, we
      // set up the straight path with as many nodes as pathNodes has, so that
      // the final variation of the path between pathNodes and straightNodes
      // can proceed.
      if( nodesConverged )
      {
        // placeholder:
        /**/std::cout << std::endl
        << "Placeholder: "
        << "This should go into pathAverager.UpdateNodes( ... )!";
        std::cout << std::endl;/**/
        straightPath.resize( pathNodes.size(),
                             pathNodes.front() );
        std::vector< double > straightNode( pathNodes.back() );
        for( size_t fieldIndex( 0 );
             fieldIndex < numberOfFields;
             ++fieldIndex )
        {
          straightNode[ fieldIndex ] -= pathNodes.front()[ fieldIndex ];
        }
        double const straightFraction( 1.0
                          / static_cast< double >( straightPath.size() - 1 ) );
        for( size_t nodeIndex( 1 );
             nodeIndex < ( straightPath.size() - 1 );
             ++nodeIndex )
        {
          for( size_t fieldIndex( 0 );
               fieldIndex < numberOfFields;
               ++fieldIndex )
          {
            straightPath[ nodeIndex ][ fieldIndex ]
            += ( nodeIndex * straightFraction * straightNode[ fieldIndex ] );
          }
        }
        straightPath.back() = pathNodes.back();
      }
      return new LinearSplineThroughNodes( pathNodes,
                                           std::vector< double >( 0 ),
                                           pathTemperature );
    }
    else
    {
      return pathAverager.ImprovePath();
    }
  }

  // This sets up reflectionMatrix to be the Householder reflection matrix
  // which reflects the axis of field 0 to be parallel to
  // currentParallelComponent.
  void MinimizingPotentialOnHypersurfaces::SetUpHouseholderReflection()
  {
    // First we check that targetVector doesn't already lie on the axis of the
    // field with index 0.
    bool alreadyParallel( true );
    for( size_t fieldIndex( 1 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( currentParallelComponent[ fieldIndex ] != 0.0 )
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
        targetLengthSquared += ( currentParallelComponent[ fieldIndex ]
                                 * currentParallelComponent[ fieldIndex ] );
      }
      double const targetNormalization( 1.0 / sqrt( targetLengthSquared ) );
      double const minusInverseOfOneMinusDotProduct( 1.0 /
              ( ( currentParallelComponent[ 0 ] * targetNormalization ) - 1.0 ) );
      reflectionMatrix = Eigen::MatrixXd::Zero( numberOfFields,
                                                numberOfFields );
      for( size_t rowIndex( 0 );
           rowIndex < numberOfFields;
           ++rowIndex )
      {
        double rowIndexPart( currentParallelComponent[ rowIndex ]
                             * targetNormalization );
        if( rowIndex == 0 )
        {
          rowIndexPart -= 1.0;
        }
        for( size_t columnIndex( rowIndex + 1 );
             columnIndex < numberOfFields;
             ++columnIndex )
        {
          double columnIndexPart( currentParallelComponent[ columnIndex ]
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

} /* namespace VevaciousPlusPlus */
