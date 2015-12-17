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
                           std::vector< BouncePathFinder* > const& pathFinders,
                                BounceActionCalculator* const actionCalculator,
                TunnelingCalculator::TunnelingStrategy const tunnelingStrategy,
                                     double const survivalProbabilityThreshold,
                               unsigned int const thermalIntegrationResolution,
                                        unsigned int const temperatureAccuracy,
                                    unsigned int const pathPotentialResolution,
                                      double const vacuumSeparationFraction ) :
    BounceActionTunneler( tunnelingStrategy,
                          survivalProbabilityThreshold,
                          temperatureAccuracy,
                          vacuumSeparationFraction ),
    pathFinders( pathFinders ),
    actionCalculator( actionCalculator ),
    thermalIntegrationResolution( thermalIntegrationResolution ),
    pathPotentialResolution( pathPotentialResolution )
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


  // This sets thermalSurvivalProbability by numerically integrating up to the
  // critical temperature for tunneling to be possible from T = 0 unless the
  // integral already passes a threshold, and sets
  // dominantTemperatureInGigaElectronVolts to be the temperature with the
  // lowest survival probability.
  void BounceAlongPathWithThreshold::ContinueThermalTunneling(
                                    PotentialFunction const& potentialFunction,
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                              double const potentialAtOriginAtZeroTemperature )
  {
    // First we set up the (square of the) threshold distance that we demand
    // between the vacua at every temperature to trust the tunneling
    // calculation.
    double const thresholdSeparationSquared( vacuumSeparationFractionSquared
                                * falseVacuum.SquareDistanceTo( trueVacuum ) );

    // We sum up decay widths over increasing temperatures.
    double partialDecayWidth( 0.0 );
    // The partial decay width scaled by the volume of the observable Universe
    // is recorded in partialDecayWidth so that the bounce action threshold for
    // each temperature can be calculated taking into account the contributions
    // from higher temperatures.
    double const temperatureStep( rangeOfMaxTemperatureForOriginToTrue.first
                 / static_cast< double >( thermalIntegrationResolution + 1 ) );
    double currentTemperature( 0.0 );
    MinuitPotentialMinimizer thermalPotentialMinimizer( potentialFunction );
    thermalPotentialMinimizer.SetTemperature( currentTemperature );
    PotentialMinimum thermalFalseVacuum( falseVacuum );
    PotentialMinimum thermalTrueVacuum( trueVacuum );
    double const thresholdDecayWidth( -log( survivalProbabilityThreshold )
                 / ( temperatureStep * exp( lnOfThermalIntegrationFactor ) ) );

    double smallestExponent( maximumPowerOfNaturalExponent );
    dominantTemperatureInGigaElectronVolts = 0.0;

    for( unsigned int whichStep( 0 );
         whichStep < thermalIntegrationResolution;
         ++whichStep )
    {
      currentTemperature += temperatureStep;
      thermalPotentialMinimizer.SetTemperature( currentTemperature );
      // We update the positions of the thermal vacua based on their positions
      // at the last temperature step.
      thermalFalseVacuum
      = thermalPotentialMinimizer( thermalFalseVacuum.FieldConfiguration() );
      // We have to keep checking to see if the field origin should be the
      // thermal false vacuum. The result of thermalPotentialMinimizer already
      // has the value of the potential at the field origin subtracted, so we
      // just compare with zero.
      if( thermalFalseVacuum.PotentialValue() > 0.0 )
      {
        thermalFalseVacuum
        = PotentialMinimum( potentialFunction.FieldValuesOrigin(),
                            0.0 );
      }
      thermalTrueVacuum
      = thermalPotentialMinimizer( thermalTrueVacuum.FieldConfiguration() );

      if( !( thermalTrueVacuum.FunctionValue()
             < thermalFalseVacuum.FunctionValue() ) )
      {
        // If the thermal vacua are the wrong way around in depth order
        // (possibly due to the thermal minimizer failing to converge properly)
        // then we break without adding in any decay width for this
        // temperature.
        std::stringstream warningBuilder;
        warningBuilder << "At temperature " << currentTemperature
        << " GeV, minimizer rolled from panic vacuum from a lower temperature"
        << " to configuration which is not deeper than the configuration to"
        << " which it rolled from the DSB vacuum from that lower temperature."
        << " Skipping the contribution of this temperature.";
        WarningLogger::LogWarning( warningBuilder.str() );
        continue;
      }
      else if( thermalTrueVacuum.SquareDistanceTo( thermalFalseVacuum )
          < thresholdSeparationSquared )
      {
        // If the thermal vacua have gotten so close that a tunneling
        // calculation is suspect, we break without adding in any decay width
        // for this temperature.
        std::stringstream warningBuilder;
        warningBuilder << "At temperature " << currentTemperature
        << " GeV, minimizer found DSB vacuum and panic vacuum to be so close"
        << " that a tunneling calculation is not trustworthy. Skipping the"
        << " contribution of this temperature and higher temperatures.";
        WarningLogger::LogWarning( warningBuilder.str() );
        break;
      }

      double const actionThreshold( -currentTemperature
                                    * log( currentTemperature
                                           * currentTemperature
                                           * ( thresholdDecayWidth
                                               - partialDecayWidth ) ) );
      double const
      bounceOverTemperature( BoundedBounceAction( potentialFunction,
                                                  thermalFalseVacuum,
                                                  thermalTrueVacuum,
                                                  currentTemperature,
                                                  actionThreshold,
                                                  thresholdSeparationSquared )
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

      if( partialDecayWidth > thresholdDecayWidth )
      {
        // We don't bother calculating the rest of the contributions to the
        // integral of the decay width if it is already large enough that the
        // survival probability is below the threshold.
        break;
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
      std::stringstream warningBuilder;
      warningBuilder
      << "The calculated integrated thermal decay width was so close to zero"
      << " that taking its logarithm would be problematic, so setting the"
      << " logarithm of the negative of the logarithm of the thermal survival"
      << " probability to " << logOfMinusLogOfThermalProbability << ".";
      WarningLogger::LogWarning( warningBuilder.str() );
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
                                    PotentialFunction const& potentialFunction,
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                             double const tunnelingTemperature,
                                                  double const actionThreshold,
                                 double const requiredVacuumSeparationSquared )
  {
    std::vector< std::vector< double > > straightPath( 2,
                                            falseVacuum.FieldConfiguration() );
    straightPath.back() = trueVacuum.FieldConfiguration();
    TunnelPath const* bestPath( new LinearSplineThroughNodes( straightPath,
                                                    std::vector< double >( 0 ),
                                                      tunnelingTemperature ) );

    actionCalculator->ResetVacua( potentialFunction,
                                  falseVacuum,
                                  trueVacuum,
                                  tunnelingTemperature );

    SplinePotential pathPotential( potentialFunction,
                                   *bestPath,
                                   pathPotentialResolution,
                                   requiredVacuumSeparationSquared );

    if( !(pathPotential.EnergyBarrierWasResolved()) )
    {
      std::stringstream warningBuilder;
      warningBuilder << "Unable to resolve an energy barrier between false"
      << " vacuum and true vacuum: returning bounce action of zero (which"
      << " should be sufficient to exclude the parameter point).";
      WarningLogger::LogWarning( warningBuilder.str() );
      return 0.0;
    }

    BubbleProfile const* bestBubble( (*actionCalculator)( *bestPath,
                                                          pathPotential ) );

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

    for( std::vector< BouncePathFinder* >::iterator
         pathFinder( pathFinders.begin() );
         pathFinder < pathFinders.end();
         ++pathFinder )
    {
      std::cout << std::endl
      << "Passing best path so far to next path finder.";
      std::cout << std::endl;

      (*pathFinder)->SetPotentialAndVacuaAndTemperature( potentialFunction,
                                                         falseVacuum,
                                                         trueVacuum,
                                                        tunnelingTemperature );
      TunnelPath const* currentPath( bestPath );
      BubbleProfile const* currentBubble( bestBubble );

      // The paths produced in sequence by pathFinder are kept separate from
      // bestPath to give more freedom to pathFinder internally (though I
      // cannot right now think of any way in which it would actually be
      // useful, apart from maybe some crazy MCMC which might try to get out of
      // a local minimum, which wouldn't work if it was sent back to the local
      // minimum at each step).

      // Keeping track of a best path and bubble separately from the last-used
      // path and bubble without copying any instances requires a bit of
      // book-keeping. Each iteration of the loop below will produce new
      // instances of a path and a bubble, and either the new path and bubble
      // need to be deleted or the previous best need to be deleted. These two
      // pointers keep track of which is to be deleted. They start as NULL,
      // allowing the first iteration to delete them before any path has been
      // marked for deletion.
      TunnelPath const* pathDeleter( NULL );
      BubbleProfile const* bubbleDeleter( NULL );

      // This loop will get a path from pathFinder and then repeat if
      // pathFinder decides that the path can be improved once the bubble
      // profile is obtained, as long as the bounce action has not dropped
      // below the threshold.
      do
      {
        // The nextPath and nextBubble pointers are not strictly necessary,
        // but they make the logic of the code clearer and will probably be
        // optimized away by the compiler anyway.
        TunnelPath const*
        nextPath( (*pathFinder)->TryToImprovePath( *currentPath,
                                                   *currentBubble ) );

        SplinePotential potentialApproximation( potentialFunction,
                                                *nextPath,
                                                pathPotentialResolution,
                                             requiredVacuumSeparationSquared );

        BubbleProfile const* nextBubble( (*actionCalculator)( *nextPath,
                                                    potentialApproximation ) );

        // On the first iteration of the loop, these pointers are NULL, so
        // it's no problem to delete them. On subsequent iterations, they point
        // at the higher-action path and bubble from the last iterations
        // comparison between its nextPath and bestPath.
        delete bubbleDeleter;
        delete pathDeleter;

        if( nextBubble->BounceAction() < bestBubble->BounceAction() )
        {
          // If nextBubble was an improvement on bestBubble, what bestPath
          // currently points at gets marked for deletion either on the next
          // iteration of this loop or just after the loop, and then bestPath
          // is set to point at nextPath, with the corresponding operations for
          // the bubble pointers.
          bubbleDeleter = bestBubble;
          bestBubble = nextBubble;
          pathDeleter = bestPath;
          bestPath = nextPath;
        }
        else
        {
          // If nextBubble wasn't an improvement on bestBubble, it and nextPath
          // will be deleted after being used to generate the nextPath and
          // nextBubble of the next iteration of the loop (or after the loop if
          // this ends up being the last iteration) through these pointers.
          bubbleDeleter = nextBubble;
          pathDeleter = nextPath;
        }
        currentBubble = nextBubble;
        currentPath = nextPath;

        std::cout << std::endl
        << "bounce action for new path = " << currentBubble->BounceAction();
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
      } while( ( bestBubble->BounceAction() > actionThreshold )
               &&
               (*pathFinder)->PathCanBeImproved( *currentBubble ) );
      // At the end of the loop, these point at the last tried path and bubble
      // which did not end up as the best ones, so deleting them is no problem.
      delete bubbleDeleter;
      delete pathDeleter;

      // We don't bother with the rest of the path finders if the action has
      // already dropped below the threshold.
      if( bestBubble->BounceAction() < actionThreshold )
      {
        std::cout
        << std::endl
        << "Bounce action dropped below threshold, breaking off from looking"
        << " for further path improvements.";
        std::cout << std::endl;

        break;
      }
    }

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

    double const bounceAction( bestBubble->BounceAction() );
    delete bestBubble;
    delete bestPath;
    return bounceAction;
  }

} /* namespace VevaciousPlusPlus */
