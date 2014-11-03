/*
 * BounceAlongPathWithThreshold.cpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/BounceActionTunneling/BounceAlongPathWithThreshold.hpp"

namespace VevaciousPlusPlus
{

  BounceAlongPathWithThreshold::BounceAlongPathWithThreshold(
                                          PotentialFunction& potentialFunction,
                           std::vector< BouncePathFinder* > const& pathFinders,
                                BounceActionCalculator* const actionCalculator,
                TunnelingCalculator::TunnelingStrategy const tunnelingStrategy,
                                     double const survivalProbabilityThreshold,
                                     size_t const thermalIntegrationResolution,
                                           size_t const temperatureAccuracy ) :
    BounceActionTunneler( potentialFunction,
                          tunnelingStrategy,
                          survivalProbabilityThreshold,
                          temperatureAccuracy ),
    pathFinders( pathFinders ),
    actionCalculator( actionCalculator ),
    thermalIntegrationResolution( thermalIntegrationResolution )
  {
    // This constructor is just an initialization list.
  }

  BounceAlongPathWithThreshold::~BounceAlongPathWithThreshold()
  {
    delete actionCalculator;

    for( size_t deletionIndex( 0 );
         deletionIndex < pathFinders.size();
         ++deletionIndex )
    {
      delete pathFinders[ deletionIndex ];
    }
  }


  // This sets thermalSurvivalProbability by numerically integrating from the
  // critical temperature for tunneling to be possible down to T = 0 unless
  // the integral already passes a threshold, and sets
  // dominantTemperatureInGigaElectronVolts to be the temperature with the
  // lowest survival probability.
  void BounceAlongPathWithThreshold::ContinueThermalTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum,
                              double const potentialAtOriginAtZeroTemperature )
  {
    // Now we start at almost the critical temperature, and sum up for
    // decreasing temperatures:
    double partialDecayWidth( 0.0 );
    // The partial decay width scaled by the volume of the observable Universe
    // is recorded in partialDecayWidth so that the bounce action threshold for
    // each temperature can be calculated taking into account the contributions
    // from higher temperatures.
    double const temperatureStep( rangeOfMaxTemperatureForOriginToTrue.first
                 / static_cast< double >( thermalIntegrationResolution + 1 ) );
    double currentTemperature( rangeOfMaxTemperatureForOriginToTrue.first
                               - temperatureStep );
    thermalPotentialMinimizer.SetTemperature( currentTemperature );
    PotentialMinimum thermalFalseVacuum( potentialFunction.FieldValuesOrigin(),
                      potentialFunction( potentialFunction.FieldValuesOrigin(),
                                                           currentTemperature )
                                - thermalPotentialMinimizer.FunctionOffset() );
    if( currentTemperature < rangeOfMaxTemperatureForOriginToFalse.second )
    {
      thermalFalseVacuum = thermalPotentialMinimizer(
                                            falseVacuum.FieldConfiguration() );
    }
    double const thresholdDecayWidth( -log( survivalProbabilityThreshold )
                 / ( temperatureStep * exp( lnOfThermalIntegrationFactor ) ) );
    double actionThreshold( -currentTemperature
      * log( currentTemperature * currentTemperature * thresholdDecayWidth ) );
    double bounceOverTemperature( BoundedBounceAction( thermalFalseVacuum,
                  thermalPotentialMinimizer( trueVacuum.FieldConfiguration() ),
                                                       currentTemperature,
                                                       actionThreshold )
                                  / currentTemperature );
    if( bounceOverTemperature < maximumPowerOfNaturalExponent )
    {
      partialDecayWidth += ( exp( -bounceOverTemperature )
                             / ( currentTemperature * currentTemperature ) );
    }

    double smallestExponent( bounceOverTemperature );
    dominantTemperatureInGigaElectronVolts = 0.0;

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "currentTemperature = " << currentTemperature
    << ", bounceOverTemperature = " << bounceOverTemperature
    << ", partialDecayWidth = " << partialDecayWidth
    << ", thresholdPartialWidth = " << thresholdDecayWidth;
    std::cout << std::endl;/**/

    for( size_t whichStep( 1 );
         whichStep < thermalIntegrationResolution;
         ++whichStep )
    {
      if( partialDecayWidth > thresholdDecayWidth )
      {
        // We don't bother calculating the rest of the contributions to the
        // integral of the decay width if it is already large enough that the
        // survival probability is below the threshold.
        break;
      }
      currentTemperature -= temperatureStep;
      thermalPotentialMinimizer.SetTemperature( currentTemperature );
      if( currentTemperature >= rangeOfMaxTemperatureForOriginToFalse.second )
      {
        thermalFalseVacuum
        = PotentialMinimum( potentialFunction.FieldValuesOrigin(),
                      potentialFunction( potentialFunction.FieldValuesOrigin(),
                                         currentTemperature )
                                - thermalPotentialMinimizer.FunctionOffset() );
      }
      else
      {
        thermalFalseVacuum = thermalPotentialMinimizer(
                                            falseVacuum.FieldConfiguration() );
      }
      // If the bounce action from the next step is sufficiently small, then
      // the contribution to survivalExponent could make survivalExponent >
      // thresholdDecayWidth, which would mean that the survival probability is
      // definitely lower than survivalProbabilityThreshold.
      actionThreshold = ( -currentTemperature
                          * log( currentTemperature * currentTemperature
                             * ( thresholdDecayWidth - partialDecayWidth ) ) );
      bounceOverTemperature = ( BoundedBounceAction( thermalFalseVacuum,
                  thermalPotentialMinimizer( trueVacuum.FieldConfiguration() ),
                                                     currentTemperature,
                                                     actionThreshold )
                                / currentTemperature );

      if( bounceOverTemperature < maximumPowerOfNaturalExponent )
      {
        partialDecayWidth += ( exp( -bounceOverTemperature )
                               / ( currentTemperature * currentTemperature ) );
      }
      if( bounceOverTemperature < smallestExponent )
      {
        smallestExponent = bounceOverTemperature;
        dominantTemperatureInGigaElectronVolts = currentTemperature;
      }
    }
    if( partialDecayWidth > 0.0 )
    {
      logOfMinusLogOfThermalProbability = ( lnOfThermalIntegrationFactor
                                + log( partialDecayWidth * temperatureStep ) );
    }
    else
    {
      logOfMinusLogOfThermalProbability
      = -exp( maximumPowerOfNaturalExponent );
      std::cout
      << std::endl
      << "Warning! The calculated integrated thermal decay width was so close"
      << " to 0 that taking its logarithm would be problematic, so setting the"
      << " logarithm of the negative of the logarithm of the thermal survival"
      << " probability to " << logOfMinusLogOfThermalProbability << ".";
      std::cout << std::endl;
    }
    SetThermalSurvivalProbability();
  }

  // This returns either the dimensionless bounce action integrated over four
  // dimensions (for zero temperature) or the dimensionful bounce action
  // integrated over three dimensions (for non-zero temperature) for tunneling
  // from falseVacuum to trueVacuum at temperature tunnelingTemperature, or an
  // upper bound if the upper bound drops below actionThreshold during the
  // course of the calculation. The vacua are assumed to already be the minima
  // at tunnelingTemperature.
  double BounceAlongPathWithThreshold::BoundedBounceAction(
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                             double const tunnelingTemperature,
                                           double const actionThreshold ) const
  {
    std::vector< std::vector< double > > straightPath( 2,
                                            falseVacuum.FieldConfiguration() );
    straightPath.back() = trueVacuum.FieldConfiguration();
    TunnelPath const* bestPath( new LinearSplineThroughNodes( straightPath,
                                                    std::vector< double >( 0 ),
                                                      tunnelingTemperature ) );

    actionCalculator->ResetVacua( falseVacuum,
                                  trueVacuum,
                                  tunnelingTemperature );
    BubbleProfile const* bestBubble( (*actionCalculator)( *bestPath ) );

    std::cout << std::endl
    << "Initial path bounce action = " << bestBubble->BounceAction();
    if( bestPath->NonZeroTemperature() )
    {
      std::cout << " GeV";
    }
    std::cout << ", threshold is " << actionThreshold;
    if( bestPath->NonZeroTemperature() )
    {
      std::cout << " GeV";
    }
    std::cout << ".";
    std::cout << std::endl;

    // debugging:
    /*std::string straightPathPicture( "StraightBubbleProfile.eps" );
    std::cout << std::endl << "debugging:"
    << std::endl
    << "Initial straight path being plotted in " << straightPathPicture << ".";
    std::cout << std::endl;
    std::vector< std::string > fieldColors;
    fieldColors.push_back( "red" );
    fieldColors.push_back( "brown" );
    fieldColors.push_back( "blue" );
    fieldColors.push_back( "purple" );
    fieldColors.push_back( "green" );
    fieldColors.push_back( "cyan" );
    actionCalculator->PlotBounceConfiguration( *bestPath,
                                               fieldColors,
                                               straightPathPicture );*/

    for( std::vector< BouncePathFinder* >::iterator
         pathFinder( pathFinders.begin() );
         pathFinder < pathFinders.end();
         ++pathFinder )
    {
      (*pathFinder)->SetVacuaAndTemperature( falseVacuum,
                                             trueVacuum,
                                             tunnelingTemperature );
      TunnelPath const* currentPath( bestPath );
      BubbleProfile const* currentBubble( bestBubble );

      // Keeping track of a best path and bubble separately from the last-used
      // path and bubble without copying any instances requires a bit of
      // book-keeping. These two pointers get set to NULL if the current path
      // and bubble happen to also be the best ones so far, so they can be
      // deleted without problem, otherwise they delete the (not best) current
      // path and bubble once they have been used to create the next path and
      // bubble.
      TunnelPath const* currentPathDeleter( NULL );
      BubbleProfile const* currentBubbleDeleter( NULL );

      while( ( bestBubble->BounceAction() > actionThreshold )
             &&
             (*pathFinder)->PathCanBeImproved( *currentBubble ) )
      {
        // The nextPath and nextBubble pointers are not strictly necessary,
        // but they make the logic of the code clearer and will probably be
        // optimized away by the compiler anyway.
        TunnelPath const*
        nextPath( (*pathFinder)->TryToImprovePath( *currentPath,
                                                   *currentBubble ) );
        BubbleProfile const* nextBubble( (*actionCalculator)( *nextPath ) );
        delete currentBubbleDeleter;
        delete currentPathDeleter;

        if( nextBubble->BounceAction() < bestBubble->BounceAction() )
        {
          delete bestBubble;
          bestBubble = nextBubble;
          currentBubble = nextBubble;
          currentBubbleDeleter = NULL;
          delete bestPath;
          bestPath = nextPath;
          currentPath = nextPath;
          currentPathDeleter = NULL;
        }
        else
        {
          currentBubble = nextBubble;
          currentPath = nextPath;

          // If nextBubble wasn't an improvement on bestBubble, it and nextPath
          // will be deleted after being used to generate the nextPath and
          // nextBubble of the next iteration of the loop (or after the loop if
          // this ends up being the last iteration) through these pointers.
          currentBubbleDeleter = nextBubble;
          currentPathDeleter = nextPath;
        }

        std::cout << std::endl
        << "Improved path bounce action = " << currentBubble->BounceAction();
        if( currentPath->NonZeroTemperature() )
        {
          std::cout << " GeV";
        }
        std::cout << ", lowest bounce action so far = "
        << bestBubble->BounceAction();
        if( currentPath->NonZeroTemperature() )
        {
          std::cout << " GeV";
        }
        std::cout << ", threshold is " << actionThreshold;
        if( currentPath->NonZeroTemperature() )
        {
          std::cout << " GeV";
        }
        std::cout << ".";
        std::cout << std::endl;
      }
      // At the end of the loop, these either point at the last tried path and
      // bubble which did not end up as the best ones, so they should be
      // deleted, or they are NULL (because the last tried path and bubble
      // happened to be the best ones), so deleting them is no problem.
      delete currentBubbleDeleter;
      delete currentPathDeleter;
    }
    // debugging:
    /*std::string finalPathPicture( "FinalBubbleProfile.eps" );
    std::cout << std::endl << "debugging:"
    << std::endl
    << "Final deformed path being plotted in " << finalPathPicture << ".";
    std::cout << std::endl;
    actionCalculator->PlotBounceConfiguration( *bestPath,
                                               fieldColors,
                                               finalPathPicture );*/

    std::cout << std::endl
    << "Lowest path bounce action at " << tunnelingTemperature << " GeV was "
    << bestBubble->BounceAction();
    if( bestPath->NonZeroTemperature() )
    {
      std::cout << " GeV";
    }
    std::cout << ", threshold is " << actionThreshold;
    if( bestPath->NonZeroTemperature() )
    {
      std::cout << " GeV";
    }
    std::cout << ".";
    std::cout << std::endl;

    delete bestBubble;
    delete bestPath;
    return bestBubble->BounceAction();
  }

} /* namespace VevaciousPlusPlus */
