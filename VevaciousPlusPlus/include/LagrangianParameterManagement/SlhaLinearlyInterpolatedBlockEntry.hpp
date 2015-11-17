/*
 * SlhaLinearlyInterpolatedBlockEntry.hpp
 *
 *  Created on: Oct 28, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHALINEARLYINTERPOLATEDBLOCKENTRY_HPP_
#define SLHALINEARLYINTERPOLATEDBLOCKENTRY_HPP_

#include "CommonIncludes.hpp"
#include "SlhaInterpolatedParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaLinearlyInterpolatedBlockEntry :
                                    public SlhaInterpolatedParameterFunctionoid
  {
  public:
    template< typename secondType > static bool FirstPairDotFirstIsLower(
                              std::pair< double, secondType > const& firstPair,
                            std::pair< double, secondType > const& secondPair )
    { return ( firstPair.first < secondPair.first ); }

    SlhaLinearlyInterpolatedBlockEntry( size_t const indexInValuesVector,
                              LHPC::SlhaSimplisticInterpreter const& lhaParser,
                                        std::string const& parameterName );
    SlhaLinearlyInterpolatedBlockEntry(
                        SlhaLinearlyInterpolatedBlockEntry const& copySource );
    virtual ~SlhaLinearlyInterpolatedBlockEntry();


    // This returns the value of the functionoid for the given logarithm of the
    // scale.
    virtual double operator()( double const logarithmOfScale ) const
    { return InterpolateOrExtrapolate( logarithmOfScale ); }

    // This returns the value of the functionoid for the given logarithm of the
    // scale. It ignores the values of the other parameters.
    virtual double operator()( double const logarithmOfScale,
                       std::vector< double > const& interpolatedValues ) const
    { return InterpolateOrExtrapolate( logarithmOfScale ); }

    // This re-assigns the vector of values paired with logarithms of the
    // block's scale according to the current status of the block.
    virtual void UpdateForNewSlhaParameters();

    // This is for creating a Python version of the potential.
    virtual std::string
    PythonParameterEvaluation( int const indentationSpaces ) const;

    // This is mainly for debugging.
    virtual std::string AsDebuggingString() const;


  protected:
    std::vector< std::pair< double, double > > logScalesWithValues;
    size_t lastIndex;


    // This finds the index of the smallest logarithm of the scale which is
    // larger than logarithmOfScale and returns the interpolation (or possibly
    // extrapolation) based on the logarithm-value pairs at the found index and
    // the index just before it. It starts with index 1 so that there is always
    // an index just before.
    double InterpolateOrExtrapolate( double const logarithmOfScale ) const;

    // This takes the points in logScalesWithValues at indexOfGreaterLog and
    // one before it, and interpolates (or possibly extrapolates) the value
    // for logarithmOfScale on a straight line through the points.
    double InterpolateOrExtrapolate( size_t const indexOfGreaterLog,
                                     double const logarithmOfScale ) const
    { return ( logScalesWithValues[ indexOfGreaterLog ].second
               + ( ( ( logScalesWithValues[ indexOfGreaterLog ].second
                       - logScalesWithValues[ indexOfGreaterLog - 1 ].second )
                     / ( logScalesWithValues[ indexOfGreaterLog ].first
                       - logScalesWithValues[ indexOfGreaterLog - 1 ].first ) )
                   * ( logarithmOfScale
                      - logScalesWithValues[ indexOfGreaterLog ].first ) ) ); }
  };





  // This is mainly for debugging.
  inline std::string
  SlhaLinearlyInterpolatedBlockEntry::AsDebuggingString() const
  {
    std::stringstream stringBuilder;
    stringBuilder
    << "lastIndex = " << lastIndex << ", logScalesWithValues = { [ ";
    for( std::vector< std::pair< double, double > >::const_iterator
         logScaleWithValue( logScalesWithValues.begin() );
         logScaleWithValue < logScalesWithValues.end();
         ++logScaleWithValue )
    {
      if( logScaleWithValue != logScalesWithValues.begin() )
      {
        stringBuilder << " ], [ ";
      }
      stringBuilder
      << logScaleWithValue->first << " , " << logScaleWithValue->second;
    }
    stringBuilder << " ] }";
    return stringBuilder.str();
  }

  // This finds the index of the smallest logarithm of the scale which is
  // larger than logarithmOfScale and returns the interpolation (or possibly
  // extrapolation) based on the logarithm-value pairs at the found index and
  // the index just before it. It starts with index 1 so that there is always
  // an index just before.
  inline double SlhaLinearlyInterpolatedBlockEntry::InterpolateOrExtrapolate(
                                          double const logarithmOfScale ) const
  {
    for( size_t whichIndex( 1 );
         whichIndex < lastIndex;
         ++whichIndex )
    {
      if( logarithmOfScale < logScalesWithValues[ whichIndex ].first )
      {
        return InterpolateOrExtrapolate( whichIndex,
                                         logarithmOfScale );
      }
    }

    // If the loop ends without returning a value, then we extrapolate from
    // the last 2 points.
    return InterpolateOrExtrapolate( lastIndex,
                                     logarithmOfScale );
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHALINEARLYINTERPOLATEDBLOCKENTRY_HPP_ */
