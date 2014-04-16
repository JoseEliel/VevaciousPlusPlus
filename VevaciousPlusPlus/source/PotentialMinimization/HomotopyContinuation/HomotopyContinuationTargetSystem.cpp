/*
 * HomotopyContinuationTargetSystem.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  HomotopyContinuationTargetSystem::HomotopyContinuationTargetSystem()
  {
    // This constructor is just an initialization list.
  }

  HomotopyContinuationTargetSystem::~HomotopyContinuationTargetSystem()
  {
    // This does nothing.
  }


  // This appends all sign-flip combinations of solutionConfiguration to
  // realSolutions if the combination is both purely real and solves the
  // homotopy continuation target system (checked by considering whether
  // each entry of the gradient changes sign when its variable is changed by
  // +- resolutionSize). It does not append any solutions that are too close
  // to those already in realSolutions, by being within a hypercube of side
  // resolutionSize centered on the existing entry.
  void
  HomotopyContinuationTargetSystem::AppendPureRealSolutionAndValidSignFlips(
            std::vector< std::complex< double > > const& solutionConfiguration,
                           std::vector< std::vector< double > >& realSolutions,
                                                  double const resolutionSize )
  {
    // First we return from the function having done nothing if the solution is
    // not purely real.
    for( std::vector< std::complex< double > >::const_iterator
         fieldValue( solutionConfiguration.begin() );
         fieldValue < solutionConfiguration.end();
         ++fieldValue )
    {
      if( ( fieldValue->imag() > resolutionSize )
          ||
          ( fieldValue->imag() < -resolutionSize ) )
      {
        return;
      }
    }

    // Next we make a vector of vectors with all possible sign flips, since we
    // got to this point because solutionConfiguration is purely real (or real
    // enough).
    unsigned int const numberOfValues( solutionConfiguration.size() );
    std::vector< std::complex< double > > currentValues( numberOfValues,
                                                   std::complex< double >( 0.0,
                                                                       0.0 ) );
    std::vector< std::vector< std::complex< double > > > signFlips;
    signFlips.push_back( solutionConfiguration );
    for( unsigned int flipIndex( 0 );
        flipIndex < numberOfValues;
         ++flipIndex )
    {
      unsigned int const vectorSize( signFlips.size() );
      for( unsigned int vectorIndex( 0 );
           vectorIndex < vectorSize;
           ++vectorIndex )
      {
        if( signFlips[ vectorIndex ][ flipIndex ] != 0.0 )
        {
          currentValues = signFlips[ vectorIndex ];
          currentValues[ flipIndex ] = -(currentValues[ flipIndex ]);
          signFlips.push_back( currentValues );
        }
      }
    }

    // Now we check to see what sign flips solve the homotopy continuation
    // target system (within tolerance).
    std::vector< std::complex< double > > offsetValues( currentValues );
    double positivePartialSlope( 0.0 );
    double negativePartialSlope( 0.0 );
    for( std::vector< std::vector< std::complex< double > > >::iterator
         signFlip( signFlips.begin() );
         signFlip < signFlips.end();
         ++signFlip )
    {
      bool validSolution( true );
      for( unsigned int whichIndex( 0 );
           whichIndex < numberOfValues;
           ++whichIndex )
      {
        offsetValues = *signFlip;
        offsetValues[ whichIndex ].real() += resolutionSize;
        HomotopyContinuationSystemValues( offsetValues,
                                          currentValues );
        positivePartialSlope = currentValues[ whichIndex ].real();
        offsetValues[ whichIndex ].real()
        -= ( resolutionSize + resolutionSize );
        HomotopyContinuationSystemValues( offsetValues,
                                          currentValues );
        negativePartialSlope = currentValues[ whichIndex ].real();
        validSolution = ( ( positivePartialSlope == 0.0 )
                          ||
                          ( ( negativePartialSlope
                              / positivePartialSlope ) <= 0.0 ) );
        // We leave the inner loop noting that we should not append *signFlip
        // to realSolutions if it doesn't solve the target system.
        if( !validSolution )
        {
          break;
        }
      }
      if( validSolution )
      {
        double valueDifference( 0.0 );
        for( std::vector< std::vector< double > >::iterator
             existingSolution( realSolutions.begin() );
             existingSolution < realSolutions.end();
             ++existingSolution )
        {
          validSolution = false;
          for( unsigned int valueIndex( 0 );
               valueIndex < numberOfValues;
               ++valueIndex )
          {
            valueDifference = ( (*existingSolution)[ valueIndex ]
                                - (*signFlip)[ valueIndex ].real() );
            if( valueDifference < 0.0 )
            {
              valueDifference = -valueDifference;
            }
            validSolution = ( valueDifference > resolutionSize );
            // We leave the inner loop noting that *signFlip is sufficiently
            // far away from existingSolution.
            if( validSolution )
            {
              break;
            }
          }
          // At this point, validSolution is true if *signFlip is sufficiently
          // far away from existingSolution. If not, we should leave this loop
          // over existing solutions, leaving validSolution as false.
          if( !validSolution )
          {
            break;
          }
        }
        // If validSolution is true after leaving the loop over existing
        // solutions, then we should append the real values of *signFlip to
        // realSolutions.
        if( validSolution )
        {
          realSolutions.push_back( std::vector< double >( numberOfValues ) );
          for( unsigned int valueIndex( 0 );
               valueIndex < numberOfValues;
               ++valueIndex )
          {
            realSolutions.back()[ valueIndex ]
            = (*signFlip)[ valueIndex ].real();
          }
        }
      }
    }
  }

} /* namespace VevaciousPlusPlus */
