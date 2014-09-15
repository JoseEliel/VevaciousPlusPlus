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
                                               MinuitBetweenPaths* pathRefiner,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
    MinuitPathFinder( minuitStrategy,
                      minuitToleranceFraction ),
    potentialFunction( potentialFunction ),
    numberOfFields( potentialFunction.NumberOfFieldVariables() ),
    pathNodes(),
    currentParallelComponent( numberOfFields ),
    reflectionMatrix( numberOfFields,
                      numberOfFields ),
    pathRefiner( pathRefiner ),
    nodeZeroParameterization( ( numberOfFields - 1 ),
                              0.0 )
  {
    // This constructor is just an initialization list.
  }

  MinimizingPotentialOnHypersurfaces::~MinimizingPotentialOnHypersurfaces()
  {
    delete pathRefiner;
  }



  // This allows Minuit2 to adjust the full path a set number of times to try
  // to minimize the sum of potentials at a set of nodes or bounce action
  // along the adjusted path, and then sets the path. Whether or not the last
  // call of TryToImprovePath( ... ) lowered the bounce action is given by
  // lastImprovementWorked, which can be used by a derived class to decide
  // internally whether to change strategy. If no improvement can be made, NULL
  // is returned.
  TunnelPath const* MinimizingPotentialOnHypersurfaces::TryToImprovePath(
                                             bool const lastImprovementWorked )
  {
    if( NodesCanStillBeImproved( lastImprovementWorked ) )
    {
      TryToImproveNodes();
      return new LinearSplineThroughNodes( pathNodes,
                                           nodeZeroParameterization,
                                           pathTemperature );
    }
    else
    {
      if( pathRefiner != NULL )
      {
        return pathRefiner->ImprovePath();
      }
      else
      {
        return NULL;
      }
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
