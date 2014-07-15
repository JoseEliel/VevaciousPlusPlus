/*
 * PolynomialThroughNodes.cpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/PathParameterization/PolynomialThroughNodes.hpp"

namespace VevaciousPlusPlus
{

  PolynomialThroughNodes::PolynomialThroughNodes(
                         std::vector< std::vector< double > > const& pathNodes,
                                                  double const pathTemperature,
                                    Eigen::MatrixXd const& pathStepsInverse ) :
    TunnelPath( pathNodes.front().size(),
                pathTemperature ),
    fieldPolynomials( numberOfFields ),
    firstDerivatives( numberOfFields ),
    secondDerivatives( numberOfFields )
  {
    Eigen::MatrixXd fieldValuesMatrix( pathNodes.size(),
                                       numberOfFields );
    for( size_t rowIndex( 0 );
         rowIndex < pathNodes.size();
         ++rowIndex )
    {
      for( size_t columnIndex( 0 );
           columnIndex < pathNodes[ rowIndex ].size();
           ++columnIndex )
      {
        fieldValuesMatrix( rowIndex,
                           columnIndex )
        = pathNodes[ rowIndex ][ columnIndex ];
      }
    }
    Eigen::MatrixXd coefficientMatrix( pathStepsInverse * fieldValuesMatrix );
    for( size_t columnIndex( 0 );
         columnIndex < numberOfFields;
         ++columnIndex )
    {
      std::vector< double >&
      coefficientVector( fieldPolynomials[ columnIndex ].CoefficientVector() );
      for( size_t rowIndex( 0 );
           rowIndex < pathNodes.size();
           ++rowIndex )
      {
        coefficientVector[ rowIndex ] = coefficientMatrix( rowIndex,
                                                           columnIndex );
      }
      firstDerivatives[ columnIndex ].BecomeFirstDerivativeOf(
                                             fieldPolynomials[ columnIndex ] );
      secondDerivatives[ columnIndex ].BecomeFirstDerivativeOf(
                                             firstDerivatives[ columnIndex ] );
    }
  }

  PolynomialThroughNodes::~PolynomialThroughNodes()
  {
    // This does nothing.
  }


  // This is for debugging.
  std::string PolynomialThroughNodes::AsDebuggingString() const
  {
    std::stringstream returnStream;
    returnStream
    << "PolynomialPathThroughNodes: fieldPolynomials = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        returnStream << "," << std::endl;
      }
      returnStream << fieldPolynomials[ fieldIndex ].AsDebuggingString();
    }
    returnStream
    << "}," << std::endl
    << "firstDerivatives = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        returnStream << "," << std::endl;
      }
      returnStream << firstDerivatives[ fieldIndex ].AsDebuggingString();
    }
    returnStream
    << "}," << std::endl
    << "secondDerivatives = { ";
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        returnStream << "," << std::endl;
      }
      returnStream << secondDerivatives[ fieldIndex ].AsDebuggingString();
    }
    returnStream << "}";
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
