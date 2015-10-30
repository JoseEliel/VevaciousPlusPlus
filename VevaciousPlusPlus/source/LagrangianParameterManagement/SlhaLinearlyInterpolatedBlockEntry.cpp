/*
 * SlhaLinearlyInterpolatedBlockEntry.cpp
 *
 *  Created on: Oct 28, 2015
 *      Author: bol
 */

#include "LagrangianParameterManagement/SlhaLinearlyInterpolatedBlockEntry.hpp"

namespace VevaciousPlusPlus
{

  SlhaLinearlyInterpolatedBlockEntry::SlhaLinearlyInterpolatedBlockEntry(
                                              size_t const indexInValuesVector,
                              LHPC::SlhaSimplisticInterpreter const& lhaParser,
                                           std::string const& parameterName ) :
    SlhaInterpolatedParameterFunctionoid( indexInValuesVector,
                                          lhaParser,
                                          parameterName ),
    logScalesWithValues(),
    lastIndex( 0 )
  {
    // This constructor is just an initialization list.
  }

  SlhaLinearlyInterpolatedBlockEntry::SlhaLinearlyInterpolatedBlockEntry(
                      SlhaLinearlyInterpolatedBlockEntry const& copySource  ) :
    SlhaInterpolatedParameterFunctionoid( copySource ),
    logScalesWithValues( copySource.logScalesWithValues ),
    lastIndex( copySource.lastIndex )
  {
    // This constructor is just an initialization list.
  }

  SlhaLinearlyInterpolatedBlockEntry::~SlhaLinearlyInterpolatedBlockEntry()
  {
    // This does nothing.
  }


  // This re-assigns the vector of values paired with logarithms of the
  // block's scale according to the current status of the block.
  void SlhaLinearlyInterpolatedBlockEntry::UpdateForNewSlhaParameters()
  {
    std::list< std::pair< double, std::string > >
    scalesWithStrings( lhaParser.getScalesPairedWithValues( parameterName ) );
    size_t const numberOfScales( scalesWithStrings.size() );

    // First we guard against no block found (in which case the value is set
    // to a flat zero for all scales) or only 1 block found (in which case
    // the value is set to be the read value flat across all scales).
    if( !( numberOfScales > 1 ) )
    {
      logScalesWithValues.resize( 2 );
      lastIndex = 1;
      logScalesWithValues[ 0 ].first = 0.0;
      logScalesWithValues[ 1 ].first = 1.0;

      if( numberOfScales == 0 )
      {
        logScalesWithValues[ 0 ].second = 0.0;
        logScalesWithValues[ 1 ].second = 0.0;
      }
      else
      {
        logScalesWithValues[ 0 ].second = BOL::StringParser::stringToDouble(
                                            scalesWithStrings.front().second );
        logScalesWithValues[ 1 ].second = logScalesWithValues[ 0 ].second;
      }
    }
    else
    {
      // The blocks are ordered as they were read from the SLHA file, which
      // may not necessarily be in ascending order with respect to the scale.
      scalesWithStrings.sort( &(FirstPairDotFirstIsLower< std::string >) );
      logScalesWithValues.resize( numberOfScales );
      size_t scaleIndex( 0 );
      std::list< std::pair< double, std::string > >::const_iterator
      listIterator( scalesWithStrings.begin() );
      while( scaleIndex < numberOfScales )
      {
        logScalesWithValues[ scaleIndex ].first = log( listIterator->first );
        logScalesWithValues[ scaleIndex++ ].second
        = BOL::StringParser::stringToDouble( (listIterator++)->second );
        // Post-increments on last use of each variable. Yey terseness.
      }
      lastIndex = ( numberOfScales - 1 );
    }
  }

  // This is for creating a Python version of the potential.
  std::string SlhaLinearlyInterpolatedBlockEntry::PythonParameterEvaluation(
                                            int const indentationSpaces ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << std::setprecision( 12 )
    << PythonIndent( indentationSpaces );

    for( size_t whichIndex( 1 );
         whichIndex <= lastIndex;
         ++whichIndex )
    {
      // The first iteration should print "if (...):".
      // The subsequent iterations should have newlines before printing.
      // Iterations after the first but before the last should print
      // "elif(...):" after the newline.
      // The last iterations should print "else:" after the newline.
      if( whichIndex > 1 )
      {
        stringBuilder
        << PythonNewlineThenIndent( indentationSpaces ) << "el";
      }
      if( whichIndex < lastIndex )
      {
        stringBuilder
        << "if ( lnQ < " << logScalesWithValues[ whichIndex ].first << " ):";
      }
      else
      {
        stringBuilder << "se:";
      }

      stringBuilder
      << PythonNewlineThenIndent( indentationSpaces + 2 )
      << "parameterValues[ " << IndexInValuesVector() << " ] = ( "
      << logScalesWithValues[ whichIndex ].second << " + ( ( "
      << ( logScalesWithValues[ whichIndex ].second
            - logScalesWithValues[ whichIndex - 1 ].second )
      << " * ( lnQ - " << logScalesWithValues[ whichIndex ].first
      << " ) ) / "
      << ( logScalesWithValues[ whichIndex ].first
           - logScalesWithValues[ whichIndex - 1 ].first )
      << " ) )";
      // The above should print
      // [el]if ( lnQ < [lnQi] ): [or just "else" for the last iteration]
      //   parameterValues[ n ]
      //   = ( [Vi]
      //       + ( ( [Vi-Vbefore] * ( lnQ - lnQi ) )
      //           / [lnQi-lnQbefore] ) )
    }
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
