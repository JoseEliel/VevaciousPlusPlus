/*
 * SlhaInterpolatedParameterFunctionoid.cpp
 *
 *  Created on: Oct 23, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/SlhaInterpolatedParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  SlhaInterpolatedParameterFunctionoid::SlhaInterpolatedParameterFunctionoid(
                                              size_t const indexInValuesVector,
                 LHPC::SLHA::SparseManyIndexedBlock< double > const& slhaBlock,
                                             std::string const& indexString ) :
    slhaBlock( slhaBlock ),
    indexVector( BOL::StringParser::stringToIntVector( indexString ) ),
    scaleLogarithmPowerCoefficients( std::vector< double > ( 1,
                                                             0.0 ) ),
    indexInValuesVector( indexInValuesVector )
  {
    // This constructor is just an initialization list.
  }

  SlhaInterpolatedParameterFunctionoid::SlhaInterpolatedParameterFunctionoid(
                     SlhaInterpolatedParameterFunctionoid const& copySource ) :
    slhaBlock( copySource.slhaBlock ),
    indexVector( copySource.indexVector ),
    scaleLogarithmPowerCoefficients(
                                  copySource.scaleLogarithmPowerCoefficients ),
    indexInValuesVector( copySource.indexInValuesVector )
  {
    // This constructor is just an initialization list.
  }

  SlhaInterpolatedParameterFunctionoid::~SlhaInterpolatedParameterFunctionoid()
  {
    // This does nothing.
  }


  // This re-calculates the coefficients of the polynomial of the logarithm
  // of the scale used in evaluating the functionoid.
  void SlhaInterpolatedParameterFunctionoid::UpdateForNewSlhaParameters()
  {
    // We set up a matrix equation for the coefficients of the polynomial in
    // the logarithm of the scale based on how many explicit values of the
    // parameter at different scales we have.
    size_t const numberOfScales( static_cast< size_t >(
                           slhaBlock.getNumberOfCopiesWithDifferentScale() ) );
    std::vector< double >&
    scaleCoefficients( scaleLogarithmPowerCoefficients.CoefficientVector() );
    scaleCoefficients.resize( numberOfScales );

    if( numberOfScales > 1 )
    {
      double logarithmOfScale;
      Eigen::MatrixXd scaleDependenceMatrix( numberOfScales,
                                             numberOfScales );
      Eigen::VectorXd scaleDependenceVector( numberOfScales );
      for( size_t whichScale( 0 );
           whichScale < numberOfScales;
           ++whichScale )
      {
        scaleDependenceVector( whichScale )
        = slhaBlock[ whichScale + 1 ]( indexVector );
        logarithmOfScale = log( slhaBlock[ whichScale + 1 ].getScale() );
        scaleDependenceMatrix( whichScale,
                               0 ) = 1.0;
        scaleDependenceMatrix( whichScale,
                               1 ) = logarithmOfScale;
        for( size_t whichPower( 2 );
             whichPower < numberOfScales;
             ++whichPower )
        {
          scaleDependenceMatrix( whichScale,
                                 whichPower ) = pow( logarithmOfScale,
                                                     whichPower );
        }
      }

      // Now we solve for the coefficients:
      scaleLogarithmPowerCoefficients.CopyFromEigen(
                             scaleDependenceMatrix.colPivHouseholderQr().solve(
                                                     scaleDependenceVector ) );
    }
    else
    {
      scaleLogarithmPowerCoefficients.CoefficientVector().assign( 1,
                                               slhaBlock[ 1 ]( indexVector ) );
    }
  }

} /* namespace VevaciousPlusPlus */
