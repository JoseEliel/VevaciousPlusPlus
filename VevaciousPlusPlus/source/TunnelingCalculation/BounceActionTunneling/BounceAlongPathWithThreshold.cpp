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
                                            BouncePathFinder* const pathFinder,
                                BounceActionCalculator* const actionCalculator,
                TunnelingCalculator::TunnelingStrategy const tunnelingStrategy,
                                     double const survivalProbabilityThreshold,
                                              size_t const temperatureAccuracy,
                                            size_t const evaporationResolution,
                                  size_t const thermalIntegrationResolution ) :
    BounceActionTunneler( potentialFunction,
                          tunnelingStrategy,
                          survivalProbabilityThreshold,
                          temperatureAccuracy,
                          evaporationResolution ),
    pathFinder( pathFinder ),
    actionCalculator( actionCalculator ),
    thermalIntegrationResolution( thermalIntegrationResolution )
  {
    // This constructor is just an initialization list.
  }

  BounceAlongPathWithThreshold::~BounceAlongPathWithThreshold()
  {
    delete actionCalculator;
    delete pathFinder;
  }


  // This sets thermalSurvivalProbability by numerically integrating from the
  // critical temperature for tunneling to be possible down to T = 0 unless
  // the integral already passes a threshold, and sets
  // dominantTemperatureInGigaElectronVolts to be the temperature with the
  // lowest survival probability.
  void BounceAlongPathWithThreshold::CalculateThermalTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    // First we find the temperature at which the DSB vacuum evaporates, and
    // possibly exclude the parameter point based on DSB being less deep than
    // origin.
    std::vector< double > const&
    fieldOrigin( potentialFunction.FieldValuesOrigin() );
    double const potentialAtOrigin( potentialFunction( fieldOrigin ) );
    if( potentialFunction( falseVacuum.FieldConfiguration() )
        > potentialAtOrigin )
    {
      std::cout
      << std::endl
      << "DSB vacuum has higher energy density than vacuum with no non-zero"
      << " VEVs! Assuming that it is implausible that the Universe cooled into"
      << " this false vacuum from the symmetric phase, and so setting survival"
      << " probability to 0.";
      std::cout << std::endl;
      dominantTemperatureInGigaElectronVolts = 0.0;
      thermalSurvivalProbability = 0.0;
      return;
    }

    // We note the temperatures where the DSB vacuum evaporates and where
    // tunneling from the true vacuum to the field origin becomes impossible :
    criticalRatherThanEvaporation = false;
    evaporationMinimum = falseVacuum;
    double const falseEvaporationTemperature( CriticalOrEvaporationTemperature(
                                                         potentialAtOrigin ) );

    criticalRatherThanEvaporation = true;
    criticalMinimum = trueVacuum;
    double const criticalTunnelingTemperature(
                       CriticalOrEvaporationTemperature( potentialAtOrigin ) );


    // Now we start at almost the critical temperature, and sum up for
    // decreasing temperatures:
    double survivalExponent( 0.0 );
    double const temperatureStep( criticalTunnelingTemperature
                 / static_cast< double >( thermalIntegrationResolution + 1 ) );
    double
    currentTemperature( criticalTunnelingTemperature - temperatureStep );
    thermalPotentialMinimizer.SetTemperature( currentTemperature );
    PotentialMinimum thermalFalseVacuum( potentialFunction.FieldValuesOrigin(),
                                         potentialAtOrigin );
    if( currentTemperature < falseEvaporationTemperature )
    {
      thermalFalseVacuum = thermalPotentialMinimizer(
                                            falseVacuum.FieldConfiguration() );
    }
    double bounceOverTemperature( BoundedBounceAction( thermalFalseVacuum,
                  thermalPotentialMinimizer( trueVacuum.FieldConfiguration() ),
                                                       currentTemperature,
                                                 lnOfThermalIntegrationFactor )
                                  / currentTemperature );
    if( bounceOverTemperature < maximumPowerOfNaturalExponent )
    {
      survivalExponent += ( exp( -bounceOverTemperature )
                             / ( currentTemperature * currentTemperature ) );
    }
    double smallestExponent( bounceOverTemperature );
    dominantTemperatureInGigaElectronVolts = 0.0;
    double const thresholdExponent( -log( survivalProbabilityThreshold ) );
    for( size_t whichStep( 1 );
         whichStep < thermalIntegrationResolution;
         ++whichStep )
    {
      if( survivalExponent > thresholdExponent )
      {
        // We don't bother calculating the rest of the contributions to the
        // integral of the decay width if it is already large enough that the
        // survival probability is below the threshold.
        break;
      }
      currentTemperature -= temperatureStep;
      thermalPotentialMinimizer.SetTemperature( currentTemperature );
      if( currentTemperature >= falseEvaporationTemperature )
      {
        thermalFalseVacuum = PotentialMinimum( fieldOrigin,
                                               potentialFunction( fieldOrigin,
                                                        currentTemperature ) );
      }
      else
      {
        thermalFalseVacuum = thermalPotentialMinimizer(
                                            falseVacuum.FieldConfiguration() );
      }
      // If the bounce action from the next step is sufficiently small, then
      // the contribution to survivalExponent could make survivalExponent >
      // thresholdExponent, which would mean that the survival probability is
      // definitely lower than survivalProbabilityThreshold.
      double const actionThreshold( currentTemperature
            * log( temperatureStep / ( ( thresholdExponent - survivalExponent )
                                * currentTemperature * currentTemperature) ) );
      bounceOverTemperature = ( BoundedBounceAction( thermalFalseVacuum,
                  thermalPotentialMinimizer( trueVacuum.FieldConfiguration() ),
                                                     currentTemperature,
                                                     actionThreshold )
                                / currentTemperature );

      if( bounceOverTemperature < maximumPowerOfNaturalExponent )
      {
        survivalExponent
        += ( ( temperatureStep * exp( -bounceOverTemperature ) )
             / ( currentTemperature * currentTemperature ) );
      }
      if( bounceOverTemperature < smallestExponent )
      {
        smallestExponent = bounceOverTemperature;
        dominantTemperatureInGigaElectronVolts = currentTemperature;
      }
    }
    if( survivalExponent > 0.0 )
    {
      logOfMinusLogOfThermalProbability
      = ( lnOfThermalIntegrationFactor - log( survivalExponent ) );
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
    actionCalculator->ResetVacua( falseVacuum,
                                  trueVacuum,
                                  tunnelingTemperature );
    TunnelPath const* bestPath( pathFinder->SetInitialPath( falseVacuum,
                                                            trueVacuum,
                                                      tunnelingTemperature ) );
    double bestBounceAction( (*actionCalculator)( *bestPath ) );
    double currentBounceAction( bestBounceAction );
    double lastBounceAction( 2.0 * currentBounceAction );

    // debugging:
    /**/std::string straightPathPicture( "StraightBubbleProfile.eps" );
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
                                               straightPathPicture );/**/

    while( ( bestBounceAction > actionThreshold )
           &&
           pathFinder->PathCanBeImproved() )
    {
      TunnelPath const* currentPath( pathFinder->TryToImprovePath(
                                    currentBounceAction < lastBounceAction ) );
      if( currentPath == NULL )
      {
        break;
      }
      lastBounceAction = currentBounceAction;
      currentBounceAction = (*actionCalculator)( *currentPath );

      if( currentBounceAction < bestBounceAction )
      {
        bestBounceAction = currentBounceAction;
        delete bestPath;
        bestPath = currentPath;
      }
      else
      {
        delete currentPath;
      }

      std::cout << std::endl
      << "Improved path bounce action = " << currentBounceAction;
      if( currentPath->NonZeroTemperature() )
      {
        std::cout << " GeV";
      }
      std::cout << ", lowest bounce action so far = " << bestBounceAction;
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

    // debugging:
    /**/std::string finalPathPicture( "FinalBubbleProfile.eps" );
    std::cout << std::endl << "debugging:"
    << std::endl
    << "Final deformed path being plotted in " << finalPathPicture << ".";
    std::cout << std::endl;
    actionCalculator->PlotBounceConfiguration( *bestPath,
                                               fieldColors,
                                               finalPathPicture );/**/

    std::cout << std::endl
    << "Lowest path bounce action at " << tunnelingTemperature << " GeV was "
    << bestBounceAction;
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

    delete bestPath;
    return bestBounceAction;
  }

} /* namespace VevaciousPlusPlus */
