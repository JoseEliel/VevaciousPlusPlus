/*
 * PathFromNodes.cpp
 *
 *  Created on: Jun 2, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/TunnelingCalculation.hpp"

namespace VevaciousPlusPlus
{

  PathFromNodes::PathFromNodes( size_t const numberOfFields,
                                size_t const referenceFieldIndex,
                                size_t const numberOfVaryingPathNodes ) :
     numberOfFields( numberOfFields ),
     referenceFieldIndex( referenceFieldIndex ),
     numberOfParameterizationFields( numberOfFields - 1 ),
     numberOfVaryingPathNodes( numberOfVaryingPathNodes ),
     pathStepSize( 1.0 / (double)( numberOfVaryingPathNodes + 1 ) ),
     pathStepInversion( CreatePathStepPowersInverse() )
  {
    // This constructor is just an initialization list.
  }

  PathFromNodes::~PathFromNodes()
  {
    // This does nothing.
  }


  // This turns a flattened matrix of numbers parameterizing the path from
  // the false vacuum to the true vacuum through field space. It assumes that
  // numberOfVaryingPathNodes nodes of numberOfParameterizationFields field
  // values (in the plane where the reference field is 0) are given by
  // pathParameterization, and sets fieldsAsPolynomials appropriately.
  void
  PathFromNodes::operator()( std::vector< double > const& pathParameterization,
                             std::vector< double > const& straightPath,
                             double const straightPathInverseLengthSquared,
                         std::vector< double > const& falseVacuumConfiguration,
                   std::vector< SimplePolynomial >& fieldsAsPolynomials ) const
  {
    // The nodes are taken as being in the plane with reference field = 0 and
    // as being relative to the false vacuum configuration. Now we project them
    // onto planes perpendicular to the straight path from false vacuum to true
    // vacuum. Calling the vector in field space from false vacuum to true
    // vacuum v, for each ( numberOfFields - 1 )-dimensional node vector n, we
    // project it to a final numberOfFields-dimensional field configuration p
    // by p = n + v [ i(n) - ( n.v / v^2 ) ], where i(n) is a fraction given
    // by j / ( numberOfVaryingPathNodes + 2 ) for the jth node given by
    // pathParameterization. Then we will have a set of vectors in field space
    // which are the relative displacements of the projected nodes from the
    // false vacuum.
    Eigen::MatrixXd pathNodes( ( numberOfVaryingPathNodes + 1 ),
                               numberOfFields );
    // We set the last node to be the true vacuum (relative to the false
    // vacuum):
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      pathNodes( numberOfVaryingPathNodes,
                 fieldIndex ) = straightPath[ fieldIndex ];
    }
    size_t actualFieldIndex( 0 );
    double dotProductWithStraightPath( NAN );
    double pathStepsToAdd( NAN );
    std::vector< double > nodeVector( numberOfFields );
    for( size_t nodeIndex( 0 );
         nodeIndex < numberOfVaryingPathNodes;
         nodeIndex += numberOfVaryingPathNodes )
    {
      dotProductWithStraightPath = 0.0;
      nodeVector[ referenceFieldIndex ] = 0.0;
      for( size_t parameterizationFieldIndex( 0 );
           parameterizationFieldIndex < numberOfParameterizationFields;
           ++parameterizationFieldIndex )
      {
        actualFieldIndex = parameterizationFieldIndex;
        if( parameterizationFieldIndex >= referenceFieldIndex )
        {
          ++actualFieldIndex;
        }
        nodeVector[ actualFieldIndex ]
        = pathParameterization[ nodeIndex + parameterizationFieldIndex ];
        dotProductWithStraightPath += ( nodeVector[ actualFieldIndex ]
                                          * straightPath[ actualFieldIndex ] );
        // p = n + v [ i(n) - ( n.v / v^2 ) ]
        // p[ actualFieldIndex ] is going into pathNodes( nodeIndex,
        //                                                actualFieldIndex )
        // n[ actualFieldIndex ] is nodeVector[ actualFieldIndex ]
        // v[ actualFieldIndex ] is straightPath[ actualFieldIndex ]
        // i(n) is nodeIndex * pathStepSize
      }
      pathStepsToAdd = ( ( (double)nodeIndex * pathStepSize )
                         - ( dotProductWithStraightPath
                             * straightPathInverseLengthSquared ) );
      for( unsigned int fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        pathNodes( nodeIndex,
                   fieldIndex )
        = ( nodeVector[ fieldIndex ]
            + ( straightPath[ fieldIndex ] * pathStepsToAdd ) );
      }
    }
    // Once we have a set of nodes, we might want to project them from planes
    // perpendicular to v onto maybe hyperbolae? I dunno, something to maybe
    // allow the path to leave the endpoints in directions that initially
    // increase the distance to the other endpoint... This can be put in later,
    // anyway.

    // Getting the polynomial coefficients from the nodes is as easy as a
    // simple matrix multiplication:
    Eigen::MatrixXd pathCoefficients( pathStepInversion * pathNodes );
    fieldsAsPolynomials.resize( numberOfFields,
                            SimplePolynomial( numberOfVaryingPathNodes + 2 ) );
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      std::vector< double >&
      coefficientVector( fieldsAsPolynomials[ fieldIndex ] );
      coefficientVector[ 0 ] = falseVacuumConfiguration[ fieldIndex ];
      for( unsigned int coefficientIndex( 0 );
           coefficientIndex <= numberOfVaryingPathNodes;
           ++coefficientIndex )
      {
        coefficientVector[ coefficientIndex + 1 ]
        = pathCoefficients( coefficientIndex,
                            fieldIndex );
      }
    }
  }

} /* namespace VevaciousPlusPlus */
