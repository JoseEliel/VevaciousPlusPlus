/*
 * SlhaPolynomialFitBlockEntry.cpp
 *
 *  Created on: Oct 28, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "SlhaPolynomialFitBlockEntry.hpp"

namespace VevaciousPlusPlus
{

  SlhaPolynomialFitBlockEntry::SlhaPolynomialFitBlockEntry(
                                              size_t const indexInValuesVector,
                              LHPC::SlhaSimplisticInterpreter const& lhaParser,
                                          std::string const& parameterName  ) :
    SlhaInterpolatedParameterFunctionoid( indexInValuesVector,
                                          lhaParser,
                                          parameterName ),
    scaleLogarithmPowerCoefficients( std::vector< double > ( 1,
                                                             0.0 ) )
  {
    // This constructor is just an initialization list.
  }

  SlhaPolynomialFitBlockEntry::SlhaPolynomialFitBlockEntry(
                              SlhaPolynomialFitBlockEntry const& copySource ) :
    SlhaInterpolatedParameterFunctionoid( copySource ),
    scaleLogarithmPowerCoefficients(
                                   copySource.scaleLogarithmPowerCoefficients )
  {
    // This constructor is just an initialization list.
  }

  SlhaPolynomialFitBlockEntry::~SlhaPolynomialFitBlockEntry()
  {
    // This does nothing.
  }


  // This re-calculates the coefficients of the polynomial of the logarithm
  // of the scale used in evaluating the functionoid.
  void SlhaPolynomialFitBlockEntry::UpdateForNewSlhaParameters()
  {
    // We set up a matrix equation for the coefficients of the polynomial in
    // the logarithm of the scale based on how many explicit values of the
    // parameter at different scales we have.
    std::list< std::pair< double, std::string > >
    scalesWithStrings( lhaParser.getScalesPairedWithValues( parameterName ) );
    size_t const numberOfScales( scalesWithStrings.size() );

    // First we guard against no block found (in which case the value is set
    // to a flat zero for all scales) or only 1 block found (in which case
    // the value is set to be the read value flat across all scales).
    if( !( numberOfScales > 1 ) )
    {
      scaleLogarithmPowerCoefficients.CoefficientVector().assign( 1,
                                                     (( numberOfScales == 0 ) ?
                                                          0.0 :
      BOL::StringParser::stringToDouble( scalesWithStrings.front().second )) );
    }
    else
    {
      // The blocks are ordered as they were read from the SLHA file, which
      // may not necessarily be in ascending order with respect to the scale.
      // However, the values do not need to be ordered for the matrix inversion
      // method to work.
      Eigen::MatrixXd scaleDependenceMatrix( numberOfScales,
                                             numberOfScales );
      Eigen::VectorXd scaleDependenceVector( numberOfScales );
      double logarithmOfScale;
      size_t scaleIndex( 0 );
      std::list< std::pair< double, std::string > >::const_iterator
      listIterator( scalesWithStrings.begin() );
      while( scaleIndex < numberOfScales )
      {
        logarithmOfScale = log( listIterator->first );
        scaleDependenceMatrix( scaleIndex,
                               0 ) = 1.0;
        scaleDependenceMatrix( scaleIndex,
                               1 ) = logarithmOfScale;
        for( size_t powerIndex( 2 );
             powerIndex < numberOfScales;
             ++powerIndex )
        {
          scaleDependenceMatrix( scaleIndex,
                                 powerIndex ) = pow( logarithmOfScale,
                                                     powerIndex );
        }
        scaleDependenceVector( scaleIndex++ )
        = BOL::StringParser::stringToDouble( (listIterator++)->second );
        // Post-increments on last use of each variable. Yey terseness.
      }

      // Now we solve for the coefficients:
      scaleLogarithmPowerCoefficients.CopyFromEigen(
                             scaleDependenceMatrix.colPivHouseholderQr().solve(
                                                     scaleDependenceVector ) );
    }
  }

} /* namespace VevaciousPlusPlus */
