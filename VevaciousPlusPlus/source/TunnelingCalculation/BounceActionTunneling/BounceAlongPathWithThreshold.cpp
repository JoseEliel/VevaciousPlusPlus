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
                                    PotentialFunction const& potentialFunction,
                                            std::string const& xmlArguments ) :
    BounceActionTunneler( potentialFunction,
                          xmlArguments ),
    pathFinder( NULL ),
    actionCalculator( NULL ),
    actionThreshold( NAN ),
    thermalIntegrationResolution( 5 )
  {
    BOL::AsciiXmlParser outerArgumentParser;
    BOL::AsciiXmlParser innerArgumentParser;
    outerArgumentParser.loadString( xmlArguments );
    std::string actionCalculatorType( "BubbleShootingOnSpline" );
    std::string actionCalculatorArguments( "" );
    std::string pathFinderType( "MinuitNodePotentialMinimizer" );
    std::string pathFinderArguments( "" );
    while( outerArgumentParser.readNextElement() )
    {
      if( outerArgumentParser.currentElementNameMatches(
                                                       "BouncePotentialFit" ) )
      {
        innerArgumentParser.loadString(
                       outerArgumentParser.getTrimmedCurrentElementContent() );
        innerArgumentParser.readNextElement();
        actionCalculatorType.assign(
                       innerArgumentParser.getTrimmedCurrentElementContent() );
        innerArgumentParser.readNextElement();
        actionCalculatorArguments.assign(
                       innerArgumentParser.getTrimmedCurrentElementContent() );
      }
      else if( outerArgumentParser.currentElementNameMatches(
                                                         "TunnelPathFinder" ) )
      {
        innerArgumentParser.loadString(
                       outerArgumentParser.getTrimmedCurrentElementContent() );
        innerArgumentParser.readNextElement();
        pathFinderType.assign(
                       innerArgumentParser.getTrimmedCurrentElementContent() );
        innerArgumentParser.readNextElement();
        pathFinderArguments.assign(
                       innerArgumentParser.getTrimmedCurrentElementContent() );
      }
    }
    if( pathFinderType.compare( "MinuitPathBounceMinimizer" ) == 0 )
    {
      pathFinder = new MinuitPathBounceMinimizer( potentialFunction,
                                                  pathFinderArguments );
    }
    else if( pathFinderType.compare( "MinuitPathPotentialMinimizer" ) == 0 )
    {
      pathFinder = new MinuitPathPotentialMinimizer( potentialFunction,
                                                     pathFinderArguments );
    }
    else if( pathFinderType.compare( "MinuitNodePotentialMinimizer" ) == 0 )
    {
      pathFinder = new MinuitNodePotentialMinimizer( potentialFunction,
                                                     pathFinderArguments );
    }
    if( actionCalculatorType.compare( "BubbleShootingOnSpline" ) == 0 )
    {
      actionCalculator = new BubbleShootingOnSpline( potentialFunction,
                                                   actionCalculatorArguments );
    }
  }

  BounceAlongPathWithThreshold::~BounceAlongPathWithThreshold()
  {
    delete actionCalculator;
    delete pathFinder;
  }

  // This should set thermalSurvivalProbability and
  // dominantTemperatureInGigaElectronVolts appropriately.
  void BounceAlongPathWithThreshold::CalculateThermalTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    double survivalExponent( 0.0 );

    // First we find the temperature at which the DSB vacuum evaporates, and
    // possibly exclude the parameter point based on DSB being less deep than
    // origin.
    double const potentialAtOrigin( potentialFunction(
                                     potentialFunction.FieldValuesOrigin() ) );
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
    criticalRatherThanEvaporation = false;
    double const
    thresholdSeparationSquared( 0.05 * 0.05 * falseVacuum.LengthSquared() );
    evaporationMinimum = falseVacuum;
    double const falseEvaporationTemperature( CriticalOrEvaporationTemperature(
                                                         potentialAtOrigin ) );
    criticalRatherThanEvaporation = true;
    criticalMinimum = trueVacuum;
    // We need the temperature where tunneling from the true vacuum to the
    // field origin becomes impossible.
    double const criticalTunnelingTemperature(
                       CriticalOrEvaporationTemperature( potentialAtOrigin ) );
    double const temperatureStep( criticalTunnelingTemperature
                              / (double)( thermalIntegrationResolution + 1 ) );
    double
    currentTemperature( criticalTunnelingTemperature - temperatureStep );
    double bounceOverTemperature( BounceAction(
                 thermalPotentialMinimizer( falseVacuum.FieldConfiguration() ),
                  thermalPotentialMinimizer( trueVacuum.FieldConfiguration() ),
                                                currentTemperature )
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
      // If the bounce action from the next step is sufficiently small, then
      // the contribution to survivalExponent could make survivalExponent >
      // thresholdExponent, which would mean that the survival probability is
      // definitely lower than survivalProbabilityThreshold.
      actionThreshold = ( currentTemperature
          * log( temperatureStep / ( ( thresholdExponent - survivalExponent )
                               * currentTemperature * currentTemperature) ) );
      bounceOverTemperature = ( BounceAction(
                 thermalPotentialMinimizer( falseVacuum.FieldConfiguration() ),
                  thermalPotentialMinimizer( trueVacuum.FieldConfiguration() ),
                                              currentTemperature )
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

} /* namespace VevaciousPlusPlus */
