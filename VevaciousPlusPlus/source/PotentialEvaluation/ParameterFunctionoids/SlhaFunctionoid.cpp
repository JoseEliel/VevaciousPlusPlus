/*
 * SlhaFunctionoid.cpp
 *
 *  Created on: Mar 3, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/ParameterFunctionoids/SlhaFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  SlhaFunctionoid::SlhaFunctionoid( std::string const& indexString,
                                    std::string const& creationString,
                                    std::string const& pythonParameterName ) :
    ParameterFunctionoid( creationString,
                          pythonParameterName ),
    slhaBlock( NULL ),
    indexVector( BOL::StringParser::stringToIntVector( indexString ) ),
    scaleLogarithmPowerCoefficients( std::vector< double > ( 1,
                                                             0.0 ) )
  {
    // This constructor is just an initialization list.
  }

  SlhaFunctionoid::~SlhaFunctionoid()
  {
    // This does nothing.
  }


  // This re-calculates the coefficients of the polynomial of the logarithm
  // of the scale used in evaluating the functionoid.
  void SlhaFunctionoid::UpdateForNewSlhaParameters()
  {
    // We set up a matrix equation for the coefficients of the polynomial in
    // the logarithm of the scale based on how many explicit values of the
    // parameter at different scales we have.
    unsigned int const
    numberOfScales( slhaBlock->getNumberOfCopiesWithDifferentScale() );
    std::vector< double >&
    scaleCoefficients( scaleLogarithmPowerCoefficients.CoefficientVector() );
    scaleCoefficients.resize( numberOfScales );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "[" << this->AsString()
    << "].UpdateForNewSlhaParameters() called. numberOfScales = "
    << numberOfScales;
    std::cout << std::endl;*/

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
        = (*slhaBlock)[ whichScale + 1 ]( indexVector );
        logarithmOfScale = log( (*slhaBlock)[ whichScale + 1 ].getScale() );
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

      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "scaleDependenceMatrix = " << std::endl
      << scaleDependenceMatrix << std::endl
      << "scaleDependenceVector = " << std::endl
      << scaleDependenceVector << std::endl;
      std::cout << std::endl;*/

      // Now we solve for the coefficients:
      scaleLogarithmPowerCoefficients.CopyFromEigen(
                             scaleDependenceMatrix.colPivHouseholderQr().solve(
                                                     scaleDependenceVector ) );
      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "scaleLogarithmPowerCoefficients = {";
      for( std::vector< double >::const_iterator
           coefficientValue(
                 scaleLogarithmPowerCoefficients.CoefficientVector().begin() );
           coefficientValue
           < scaleLogarithmPowerCoefficients.CoefficientVector().end();
           ++coefficientValue )
      {
        std::cout << " " << *coefficientValue;
      }
      std::cout << " }";
      std::cout << std::endl;*/
    }
    else
    {
      scaleLogarithmPowerCoefficients.CoefficientVector().assign( 1,
                                            (*slhaBlock)[ 1 ]( indexVector ) );
    }
  }

} /* namespace VevaciousPlusPlus */
