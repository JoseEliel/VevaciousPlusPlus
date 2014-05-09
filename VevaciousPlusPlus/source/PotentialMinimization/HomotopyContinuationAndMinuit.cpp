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
                                    PotentialFunction const& potentialFunction,
                        HomotopyContinuationSolver& homotopyContinuationSolver,
                              double const extremumSeparationThresholdFraction,
                               double const nonDsbRollingToDsbScalingFactor ) :
    HomotopyContinuationAndGradient( potentialFunction,
                                     homotopyContinuationSolver ),
    potentialForMinuit( potentialFunction ),
    minuitManager( potentialForMinuit ),
    extremumSeparationThresholdFraction( extremumSeparationThresholdFraction ),
    nonDsbRollingToDsbScalingFactor( nonDsbRollingToDsbScalingFactor )
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
    for( std::vector< std::vector< double > >::const_iterator
         realSolution( purelyRealSolutionSets.begin() );
         realSolution < purelyRealSolutionSets.end();
         ++realSolution )
    {
      std::cout
      << std::endl
      << "Minuit2 starting point: "
      << potentialFunction.FieldConfigurationAsMathematica( *realSolution );
      std::cout << std::endl;

      foundMinimum = minuitManager( *realSolution );

      std::cout
      << "Rolled to "
      << foundMinimum.AsMathematica( potentialFunction.FieldNames() );
      std::cout << std::endl;

      bool rolledToDsbOrSignFlip( ( foundMinimum.SquareDistanceTo( dsbVacuum )
                                    < thresholdSepartionSquared )
                                  ||
                                !( IsNotPhaseRotationOfDsbVacuum( foundMinimum,
                                                      thresholdSepartion ) ) );

      // We check to see if a homotopy continuation solution that was not the
      // DSB minimum rolled to the DSB minimum: if so, we scale the homotopy
      // continuation solution's fields by a factor and pass the scaled set of
      // field values to Minuit to roll, and then carry on based on this new
      // Minuit minimum. (We discovered in explorations with Vevacious 1 that
      // it could happen that the basin of attraction of the DSB minimum at
      // 1-loop level could grow so large that it would encompass tree-level
      // minima that belong in some sense to other 1-loop minima, which moved
      // very far away due to loop corrections, so even though their basins of
      // attraction also grew very large in the same way that of the DSB
      // minimum did, they moved enough that their tree-level minima were left
      // out.)
      if( rolledToDsbOrSignFlip
          &&
          ( dsbVacuum.SquareDistanceTo( *realSolution )
            > thresholdSepartionSquared ) )
      {
        // We don't want to bother re-rolling the field origin, so we keep note
        // of how far away from the field origin realSolution is.
        double lengthSquared( 0.0 );
        std::vector< double > scaledPoint( *realSolution );
        for( std::vector< double >::iterator
             scaledField( scaledPoint.begin() );
             scaledField < scaledPoint.end();
             ++scaledField )
        {
          lengthSquared += ( (*scaledField) * (*scaledField) );
          *scaledField *= nonDsbRollingToDsbScalingFactor;
        }
        if( lengthSquared > thresholdSepartionSquared )
        {
          std::cout
          << "Starting point from homotopy continuation rolled to the DSB"
          << " minimum, or a phase rotation, using the full potential. Trying"
          << " a scaled starting point: "
          << potentialFunction.FieldConfigurationAsMathematica( scaledPoint );
          std::cout << std::endl;

          foundMinimum = minuitManager( scaledPoint );
          rolledToDsbOrSignFlip = ( foundMinimum.SquareDistanceTo( dsbVacuum )
                                    < thresholdSepartionSquared )
                                    ||
                                !( IsNotPhaseRotationOfDsbVacuum( foundMinimum,
                                                        thresholdSepartion ) );
          std::cout
          << "Rolled to "
          << foundMinimum.AsMathematica( potentialFunction.FieldNames() );
          std::cout << std::endl;
        }
      }

      foundMinima.push_back( foundMinimum );
      if( ( ( foundMinimum.FunctionValue() + foundMinimum.FunctionError() )
            < dsbVacuum.FunctionValue() )
          &&
          !rolledToDsbOrSignFlip )
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
    }

    std::cout
    << std::endl
    << "DSB vacuum = "
    << dsbVacuum.AsMathematica( potentialFunction.FieldNames() ) << std::endl;
    if( panicVacua.empty() )
    {
      std::cout
      << "DSB vacuum is stable as far as the model file allows." << std::endl;
    }
    else
    {
      std::cout << "Panic vacuum = "
      << panicVacuum.AsMathematica( potentialFunction.FieldNames() )
      << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;


    // debugging:
    /*std::cout << std::endl << "debugging:"
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
    std::cout << std::endl;*/
  }

} /* namespace VevaciousPlusPlus */
