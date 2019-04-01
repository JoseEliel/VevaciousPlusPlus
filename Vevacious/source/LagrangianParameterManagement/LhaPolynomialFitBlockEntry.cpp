/*
 * LhaPolynomialFitBlockEntry.cpp
 *
 *  Created on: Oct 28, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/LhaPolynomialFitBlockEntry.hpp"

namespace VevaciousPlusPlus
{

  LhaPolynomialFitBlockEntry::LhaPolynomialFitBlockEntry(
                                              size_t const indexInValuesVector,
                                        LHPC::SimpleLhaParser const& lhaParser,
                                          std::string const& parameterName  ) :
    LhaInterpolatedParameterFunctionoid( indexInValuesVector,
                                         lhaParser,
                                         parameterName ),
    scaleLogarithmPowerCoefficients( std::vector< double > ( 1,
                                                             0.0 ) )
  {
    // This constructor is just an initialization list.
  }

  LhaPolynomialFitBlockEntry::LhaPolynomialFitBlockEntry(
                              LhaPolynomialFitBlockEntry const& copySource ) :
    LhaInterpolatedParameterFunctionoid( copySource ),
    scaleLogarithmPowerCoefficients(
                                   copySource.scaleLogarithmPowerCoefficients )
  {
    // This constructor is just an initialization list.
  }

  LhaPolynomialFitBlockEntry::~LhaPolynomialFitBlockEntry()
  {
    // This does nothing.
  }


  // This re-calculates the coefficients of the polynomial of the logarithm
  // of the scale used in evaluating the functionoid.
  void LhaPolynomialFitBlockEntry::UpdateForNewLhaParameters()
  {
    // We set up a matrix equation for the coefficients of the polynomial in
    // the logarithm of the scale based on how many explicit values of the
    // parameter at different scales we have.
    std::list< std::pair< std::string, double > > entriesAtScales;
    (*lhaParser)( parameterName,
                  entriesAtScales,
                  true );
    size_t const numberOfScales( entriesAtScales.size() );

    // First we guard against no block found (in which case the value is set
    // to a flat zero for all scales) or only 1 block found (in which case
    // the value is set to be the read value flat across all scales).
    if( !( numberOfScales > 1 ) )
    {
      double constantValue( 0.0 );
      if( !(entriesAtScales.empty()) )
      {
        constantValue = LHPC::ParsingUtilities::StringToDouble(
                                              entriesAtScales.back().first );
      }
      else
      {
        // If there were no entries with explicit scales, we check for entries
        // in blocks without scales, which will be assumed to be constant over
        // all scales.
        (*lhaParser)( parameterName,
                      entriesAtScales,
                      false );

        if( !(entriesAtScales.empty()) )
        {
          constantValue = LHPC::ParsingUtilities::StringToDouble(
                                                entriesAtScales.back().first );
        }
      }

      scaleLogarithmPowerCoefficients.CoefficientVector().assign( 1,
                                                               constantValue );
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
      std::list< std::pair< std::string, double > >::const_iterator
      listIterator( entriesAtScales.begin() );
      while( scaleIndex < numberOfScales )
      {
        logarithmOfScale = log( listIterator->second );
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
        = LHPC::ParsingUtilities::StringToDouble( (listIterator++)->first );
        // Post-increments on last use of each variable. Yey terseness.
      }

      // Now we solve for the coefficients:
      scaleLogarithmPowerCoefficients.CopyFromEigen(
                             scaleDependenceMatrix.colPivHouseholderQr().solve(
                                                     scaleDependenceVector ) );
    }
  }

} /* namespace VevaciousPlusPlus */
