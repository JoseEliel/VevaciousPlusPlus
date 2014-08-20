/*
 * NodesFromParameterization.cpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/NodesFromParameterization.hpp"

namespace VevaciousPlusPlus
{

  NodesFromParameterization::NodesFromParameterization(
                                                   size_t const numberOfFields,
                                     size_t const numberOfIntermediateNodes ) :
    numberOfFields( numberOfFields ),
    numberOfIntermediateNodes( numberOfIntermediateNodes ),
    pathNodes( ( numberOfIntermediateNodes + 2 ),
               std::vector< double >( numberOfFields ) ),
    zeroFullParameterization( ( ( numberOfFields - 1 )
                                * numberOfIntermediateNodes ),
                              0.0 ),
    zeroNodeParameterization( ( numberOfFields - 1 ),
                              0.0 ),
    initialStepSizes( zeroFullParameterization )
  {
    // This constructor is just an initialization list.
  }

  NodesFromParameterization::~NodesFromParameterization()
  {
    // This does nothing.
  }


  // This sets reflectionMatrix to be the Householder reflection which
  // reflects the axis of field referenceAxis to lie along targetVector.
  void
  NodesFromParameterization::SetAsHouseholderReflectionFromAxisToVector(
                                             Eigen::MatrixXd& reflectionMatrix,
                                                    size_t const referenceAxis,
                                    std::vector< double > const& targetVector )
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "NodesFromParameterization::SetAsHouseholderReflectionFromAxisToVector("
    << " reflectionMatrix = {";
    std::cout << std::endl;
    std::cout << reflectionMatrix;
    std::cout << std::endl;
    std::cout << "}, referenceAxis = " << referenceAxis
    << ", targetVector = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < targetVector.size();
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        std::cout << ", ";
      }
      std::cout << targetVector[ fieldIndex ];
    }
    std::cout << " } ) called.";
    std::cout << std::endl;/**/

    // First we check that targetVector doesn't already lie on referenceAxis.
    bool alreadyParallel( true );
    for( size_t fieldIndex( 0 );
         fieldIndex < targetVector.size();
         ++fieldIndex )
    {
      if( ( fieldIndex != referenceAxis )
          &&
          ( targetVector[ fieldIndex ] != 0.0 ) )
      {
        alreadyParallel = false;
        break;
      }
    }

    if( alreadyParallel )
    {
      reflectionMatrix = Eigen::MatrixXd::Identity( numberOfFields,
                                                    numberOfFields );
      reflectionMatrix( referenceAxis, referenceAxis ) = -1.0;
    }
    else
    {
      double targetLengthSquared( 0.0 );
      for( size_t fieldIndex( 0 );
           fieldIndex < targetVector.size();
           ++fieldIndex )
      {
        targetLengthSquared += ( targetVector[ fieldIndex ]
                                 * targetVector[ fieldIndex ] );
      }
      double const targetNormalization( 1.0 / sqrt( targetLengthSquared ) );
      double const
      minusInverseOfOneMinusDotProduct( 1.0
                                        / ( ( targetVector[ referenceAxis ]
                                              * targetNormalization )
                                            - 1.0 ) );
      reflectionMatrix = Eigen::MatrixXd::Zero( numberOfFields,
                                                numberOfFields );
      for( size_t rowIndex( 0 );
           rowIndex < numberOfFields;
           ++rowIndex )
      {
        double
        rowIndexPart( targetVector[ rowIndex ] * targetNormalization );
        if( rowIndex == referenceAxis )
        {
          rowIndexPart -= 1.0;
        }
        for( size_t columnIndex( rowIndex + 1 );
             columnIndex < numberOfFields;
             ++columnIndex )
        {
          double columnIndexPart( targetVector[ columnIndex ]
                                  * targetNormalization );
          if( columnIndex == referenceAxis )
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
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "reflectionMatrix set to" << std::endl;
      std::cout << reflectionMatrix;
      Eigen::VectorXd targetEigen( numberOfFields );
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        targetEigen( fieldIndex ) = ( targetVector[ fieldIndex ]
                                      * targetNormalization );
      }
      std::cout << std::endl;
      std::cout << "targetEigen =" << std::endl;
      std::cout << targetEigen;
      std::cout << std::endl;
      std::cout << "reflectionMatrix * targetEigen =" << std::endl;
      std::cout << ( reflectionMatrix * targetEigen );
      std::cout << std::endl;/**/
    }
  }

} /* namespace VevaciousPlusPlus */
