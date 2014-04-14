/*
 * HomotopyContinuationAndMinuit.cpp
 *
 *  Created on: Apr 14, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  HomotopyContinuationAndMinuit::HomotopyContinuationAndMinuit(
                                          PotentialFunction& potentialFunction,
                        HomotopyContinuationSolver& homotopyContinuationSolver,
                           double const extremumSeparationThresholdFraction ) :
    HomotopyContinuationAndGradient( potentialFunction,
                                     homotopyContinuationSolver ),
    potentialForMinuit( potentialFunction ),
    minuitManager( potentialForMinuit ),
    extremumSeparationThresholdFraction( extremumSeparationThresholdFraction )
  {
    // This constructor is just an initialization list.
  }

  HomotopyContinuationAndMinuit::~HomotopyContinuationAndMinuit()
  {
    // This does nothing.
  }


  // This uses Minuit2 to minimize potentialForMinuit starting from the
  // values in purelyRealSolutionSets.
  void HomotopyContinuationAndMinuit::RollAndSortExtrema()
  {
    double const
    thresholdSepartionSquared( ( extremumSeparationThresholdFraction
                                 * extremumSeparationThresholdFraction
                                 * dsbVacuum.LengthSquared() ) + 1.0 );
    double const thresholdSepartion( sqrt( thresholdSepartionSquared ) );
    PotentialMinimum foundMinimum;
    for( std::vector< std::vector< double > >::iterator
         realSolution( purelyRealSolutionSets.begin() );
         realSolution < purelyRealSolutionSets.end();
         ++realSolution )
    {
      foundMinimum = minuitManager( *realSolution );
      foundMinima.push_back( foundMinimum );

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "found [" << foundMinimum.AsDebuggingString() << "]"
      << std::endl;
      std::cout << std::endl;/**/

      if( ( ( foundMinimum.FunctionValue() + foundMinimum.FunctionError() )
            < dsbVacuum.FunctionValue() )
          &&
          ( foundMinimum.SquareDistanceTo( dsbVacuum )
            > thresholdSepartionSquared )
          &&
          ( IsNotPhaseRotationOfDsbVacuum( foundMinimum,
                                           thresholdSepartion ) ) )
      {
        if( panicVacua.empty()
            ||
            ( foundMinimum.SquareDistanceTo( dsbVacuum )
              < panicVacuum.SquareDistanceTo( dsbVacuum ) ) )
        {
          panicVacuum = foundMinimum;
        }
        panicVacua.push_back( foundMinimum );
      }
      else
      {
        // placeholder:
        /**/std::cout << std::endl
        << "Placeholder: "
        << "Should nudge this extremum if it's a saddle point.";
        std::cout << std::endl;/**/
      }
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "dsbVacuum = " << dsbVacuum.AsDebuggingString()
    << std::endl
    << "panicVacuum = " << panicVacuum.AsDebuggingString()
    << std::endl
    << "panicVacua.size() = " << panicVacua.size();
    for( unsigned int panicIndex( 0 );
         panicIndex < panicVacua.size();
         ++panicIndex )
    {
      std::cout << std::endl << "panicVacua[ " << panicIndex << " ] = "
      << panicVacua[ panicIndex ].AsDebuggingString();
    }
    std::cout << std::endl;
    std::cout << std::endl;/**/
  }
} /* namespace VevaciousPlusPlus */
