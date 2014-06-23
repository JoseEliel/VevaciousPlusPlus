/*
 * HomotopyContinuationTargetSystem.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  HomotopyContinuationTargetSystem::HomotopyContinuationTargetSystem(
                                   SlhaUpdatePropagator& previousPropagator ) :
    SlhaUpdatePropagator( previousPropagator )
  {
    // This constructor is just an initialization list.
  }

  HomotopyContinuationTargetSystem::HomotopyContinuationTargetSystem(
                                                   SlhaManager& slhaManager ) :
    SlhaUpdatePropagator( slhaManager )
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
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "HomotopyContinuationTargetSystem::"
    << "AppendPureRealSolutionAndValidSignFlips( {";
    for( std::vector< std::complex< double > >::const_iterator
         fieldValue( solutionConfiguration.begin() );
         fieldValue < solutionConfiguration.end();
         ++fieldValue )
    {
      std::cout << " ( " << fieldValue->real() << ", " << fieldValue->imag()
          << " )";
    }
    std::cout << " }, ... ) called.";
    std::cout << std::endl;*/

    // First we return from the function having done nothing if the solution is
    // not purely real.
    size_t const numberOfValues( solutionConfiguration.size() );
    std::vector< double > realSolution( numberOfValues );
    for( size_t fieldIndex( 0 );
         fieldIndex < solutionConfiguration.size();
         ++fieldIndex )
    {
      if( ( solutionConfiguration[ fieldIndex ].imag() > resolutionSize )
          ||
          ( solutionConfiguration[ fieldIndex ].imag() < -resolutionSize ) )
      {
        // debugging:
        /*std::cout << std::endl << "debugging:"
        << std::endl
        << "Not purely real (" << fieldValue->imag() << "i).";
        std::cout << std::endl;*/
        return;
      }
      realSolution[ fieldIndex ] = solutionConfiguration[ fieldIndex ].real();
    }

    // Next we make a vector of vectors with all possible sign flips, since we
    // got to this point because solutionConfiguration is purely real (or real
    // enough).
    std::vector< double > currentValues( numberOfValues,
                                         0.0 );
    std::vector< std::vector< double > > signFlips;
    signFlips.push_back( realSolution );
    for( size_t flipIndex( 0 );
         flipIndex < numberOfValues;
         ++flipIndex )
    {
      size_t const vectorSize( signFlips.size() );
      for( size_t vectorIndex( 0 );
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
    std::vector< double > offsetValues( currentValues );
    double positivePartialSlope( 0.0 );
    double negativePartialSlope( 0.0 );
    for( std::vector< std::vector< double > >::iterator
         signFlip( signFlips.begin() );
         signFlip < signFlips.end();
         ++signFlip )
    {
      bool validSolution( true );
      for( size_t whichIndex( 0 );
           whichIndex < numberOfValues;
           ++whichIndex )
      {
        offsetValues = *signFlip;
        offsetValues[ whichIndex ] += resolutionSize;
        HomotopyContinuationSystemValues( offsetValues,
                                          currentValues );
        positivePartialSlope = currentValues[ whichIndex ];
        offsetValues[ whichIndex ] -= ( resolutionSize + resolutionSize );
        HomotopyContinuationSystemValues( offsetValues,
                                          currentValues );
        negativePartialSlope = currentValues[ whichIndex ];
        validSolution = ( ( positivePartialSlope == 0.0 )
                          ||
                          ( ( negativePartialSlope
                              / positivePartialSlope ) <= 0.0 ) );
        // We leave the inner loop noting that we should not append *signFlip
        // to realSolutions if it doesn't solve the target system.
        if( !validSolution )
        {
          // debugging:
          /*std::cout << std::endl << "debugging:"
          << std::endl
          << "signFlip {";
          for( std::vector< double >::const_iterator
               fieldValue( signFlip->begin() );
               fieldValue < signFlip->end();
               ++fieldValue )
          {
            std::cout << " ( " << fieldValue->real() << ", "
            << fieldValue->imag() << " )";
          }
          std::cout << " } failed to solve target system.";
          std::cout << std::endl;*/

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
                                - (*signFlip)[ valueIndex ] );
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
        // solutions, then we should append *signFlip to realSolutions, unless
        // a derived class has a further veto on this solution (such as
        // requiring it to be a minimum of the tree-level potential rather than
        // just an extremum, for PolynomialGradientTarget for a
        // PotentialFromPolynomialAndMasses instance, for example).
        validSolution = AllowedSolution( *signFlip );
        if( validSolution )
        {
          realSolutions.push_back( *signFlip );
        }
      }
    }
  }

} /* namespace VevaciousPlusPlus */
