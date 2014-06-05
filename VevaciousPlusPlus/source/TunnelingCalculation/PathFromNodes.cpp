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
     pathStepInversion( CreatePathStepPowersInverse( numberOfVaryingPathNodes,
                                                     pathStepSize ) )
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
  // pathParameterization, and returns a PathFieldsAndPotential constructed
  // from the polynomial fits of the fields.
  PathFieldsAndPotential
  PathFromNodes::operator()( std::vector< double > const& pathParameterization,
                             std::vector< double > const& straightPath,
                             double const straightPathInverseLengthSquared,
                         std::vector< double > const& falseVacuumConfiguration,
                             double const falseVacuumDepth,
                             double const trueVacuumDepth,
                             double const givenTemperature ) const
  {

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "PathFromNodes::operator( pathParameterization = { ";
    for( std::vector< double >::const_iterator
         pathParameter( pathParameterization.begin() );
         pathParameter < pathParameterization.end();
         ++pathParameter )
    {
      if( pathParameter != pathParameterization.begin() )
      {
        std::cout << ", ";
      }
      std::cout << *pathParameter;
    }
    std::cout << " }," << std::endl << "  straightPath = {";
    for( std::vector< double >::const_iterator
         pathParameter( straightPath.begin() );
         pathParameter < straightPath.end();
         ++pathParameter )
    {
      if( pathParameter != straightPath.begin() )
      {
        std::cout << ", ";
      }
      std::cout << *pathParameter;
    }
    std::cout << " }," << std::endl << "  straightPathInverseLengthSquared = "
    << straightPathInverseLengthSquared << "," << std::endl
    << "  falseVacuumConfiguration = {";
    for( std::vector< double >::const_iterator
         fieldValue( falseVacuumConfiguration.begin() );
         fieldValue < falseVacuumConfiguration.end();
         ++fieldValue )
    {
      if( fieldValue != falseVacuumConfiguration.begin() )
      {
        std::cout << ", ";
      }
      std::cout << *fieldValue;
    }
    std::cout << " } ) called.";
    std::cout << std::endl;/**/

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
         ++nodeIndex )
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
        = pathParameterization[ ( nodeIndex * numberOfParameterizationFields )
                                + parameterizationFieldIndex ];
        dotProductWithStraightPath += ( nodeVector[ actualFieldIndex ]
                                        * straightPath[ actualFieldIndex ] );
        // p = n + v [ i(n) - ( n.v / v^2 ) ]
        // p[ actualFieldIndex ] is going into pathNodes( nodeIndex,
        //                                                actualFieldIndex )
        // n[ actualFieldIndex ] is nodeVector[ actualFieldIndex ]
        // v[ actualFieldIndex ] is straightPath[ actualFieldIndex ]
        // i(n) is ( ( nodeIndex + 1 ) * pathStepSize )
      }
      pathStepsToAdd = ( ( (double)( nodeIndex + 1 ) * pathStepSize )
                         - ( dotProductWithStraightPath
                             * straightPathInverseLengthSquared ) );
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "nodeIndex = " << nodeIndex << ", pathStepsToAdd = "
      << pathStepsToAdd << ", nodeVector = { ";
      for( unsigned int fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        if( fieldIndex > 0 )
        {
          std::cout << ", ";
        }
        std::cout << nodeVector[ fieldIndex ];
      }
      std::cout << " }";
      std::cout << std::endl;/**/
      for( unsigned int fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        pathNodes( nodeIndex,
                   fieldIndex )
        = ( nodeVector[ fieldIndex ]
            + ( straightPath[ fieldIndex ] * pathStepsToAdd ) );
      }
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "pathNodes( nodeIndex, . ) = { ";
      for( unsigned int fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        if( fieldIndex > 0 )
        {
          std::cout << ", ";
        }
        std::cout << pathNodes( nodeIndex,
                                fieldIndex );
      }
      std::cout << " }";
    }
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "pathNodes =" << std::endl << pathNodes;
    std::cout << std::endl;/**/
    // Once we have a set of nodes, we might want to project them from planes
    // perpendicular to v onto maybe hyperbolae? I dunno, something to maybe
    // allow the path to leave the endpoints in directions that initially
    // increase the distance to the other endpoint... This can be put in later,
    // anyway.

    // Getting the polynomial coefficients from the nodes is as easy as a
    // simple matrix multiplication:
    return PathFieldsAndPotential( ( pathStepInversion * pathNodes ),
                                   falseVacuumDepth,
                                   trueVacuumDepth,
                                   givenTemperature );
  }

} /* namespace VevaciousPlusPlus */
