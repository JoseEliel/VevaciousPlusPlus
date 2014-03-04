/*
 * SlhaFunctionoid.cpp
 *
 *  Created on: Mar 3, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  SlhaFunctionoid::SlhaFunctionoid( std::string const& indexString ) :
    ParameterFunctionoid(),
    slhaBlock( NULL ),
    indexVector( BOL::StringParser::stringToIntVector( indexString ) ),
    scaleLogarithmPowerCoefficients( 1,
                                     NAN )
  {
    // This constructor is just an initialization list.

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "SlhaFunctionoid::SlhaFunctionoid( \"" << indexString << "\" )"
    << std::endl
    << "slhaBlock = " << slhaBlock << ", indexVector =";
    for( std::vector< int >::const_iterator
         whichIndex( indexVector.begin() );
         whichIndex < indexVector.end();
         ++whichIndex )
    {
      std::cout << " " << *whichIndex;
    }
    std::cout << std::endl;/**/
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
    scaleLogarithmPowerCoefficients.resize( numberOfScales );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "[" << this->AsString()
    << "].UpdateForNewSlhaParameters() called. numberOfScales = "
    << numberOfScales;
    std::cout << std::endl;/**/

    if( numberOfScales > 1 )
    {
      double logarithmOfScale;
      Eigen::MatrixXd scaleDependenceMatrix( numberOfScales,
                                             numberOfScales );
      Eigen::VectorXd scaleDependenceVector( numberOfScales );
      for( unsigned int whichScale( 0 );
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
        for( unsigned int whichPower( 2 );
             whichPower < numberOfScales;
             ++whichPower )
        {
          scaleDependenceMatrix( whichScale,
                                 whichPower ) = pow( logarithmOfScale,
                                                     whichPower );
        }
      }

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "scaleDependenceMatrix = " << std::endl
      << scaleDependenceMatrix << std::endl
      << "scaleDependenceVector = " << std::endl
      << scaleDependenceVector << std::endl;
      std::cout << std::endl;/**/

      // Now we solve for the coefficients:
      Eigen::VectorXd
      solutionVector( scaleDependenceMatrix.colPivHouseholderQr().solve(
                                                     scaleDependenceVector ) );
      for( unsigned int whichScale( 0 );
           whichScale < numberOfScales;
           ++whichScale )
      {
        scaleLogarithmPowerCoefficients[ whichScale ]
        = solutionVector( whichScale );

        // debugging:
        /**/std::cout << std::endl << "debugging:"
        << std::endl
        << "scaleLogarithmPowerCoefficients[ " << whichScale << " ] = "
        << scaleLogarithmPowerCoefficients[ whichScale ];
        std::cout << std::endl;/**/
      }
    }
    else
    {
      scaleLogarithmPowerCoefficients.assign( 1,
                                            (*slhaBlock)[ 1 ]( indexVector ) );
    }
  }

} /* namespace VevaciousPlusPlus */
