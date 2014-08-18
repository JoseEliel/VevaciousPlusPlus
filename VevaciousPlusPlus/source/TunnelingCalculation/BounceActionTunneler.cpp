/*
 * BounceActionTunneler.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/BounceActionTunneler.hpp"

namespace VevaciousPlusPlus
{
  double const BounceActionTunneler::maximumPowerOfNaturalExponent(
                           log( 0.5 * std::numeric_limits< double >::max() ) );
  double const
  BounceActionTunneler::hBarInGigaElectronVoltSeconds( 6.58211928E-25 );
  double const BounceActionTunneler::ageOfKnownUniverseInSeconds( 4.3E+17 );
  double const
  BounceActionTunneler::ageOfKnownUniverseInInverseGigaElectronVolts(
                              BounceActionTunneler::ageOfKnownUniverseInSeconds
                       / BounceActionTunneler::hBarInGigaElectronVoltSeconds );
  double const BounceActionTunneler::fourVolumeOfKnownUniverseOverGevFourth(
             BounceActionTunneler::ageOfKnownUniverseInInverseGigaElectronVolts
           * BounceActionTunneler::ageOfKnownUniverseInInverseGigaElectronVolts
           * BounceActionTunneler::ageOfKnownUniverseInInverseGigaElectronVolts
        * BounceActionTunneler::ageOfKnownUniverseInInverseGigaElectronVolts );
  double const BounceActionTunneler::lnOfThermalIntegrationFactor( 244.53 );
  // Based on correspondence with Alexander Kusenko and discussion with Bjoern
  // Garbrecht:
  // Taking [decay width per horizon]
  // = [horizon volume] * [solitonic coefficient] * exp(-[thermal action]/T)
  // at temperature T, where [horizon volume] = ( M_Plank / T^2 )^3, and
  // taking [solitonic coefficient] to be T^4,
  // the survival probability per horizon =
  // exp( -integral of [time at T] with respect to [decay time] )
  // = exp( -integral of [decay width per horizon] dT * [factor] ) )
  // which exponents for N horizons to
  // exp( -N * integral * [factor] ) )
  // and [decay width per horizon] = M_Plank^3 T^(-2) exp(-S_3(T)/T)
  // where [thermal action at temperature T] = S_3(T)
  // exp( -N * integral * [factor] ) ) can be written, from entropy
  // conservation and so on, as
  // exp( -N * integral of C T^(-2) exp(-S_3(T)/T) dT ) )
  // where C = [reduced Planck mass] * [solitonic coefficient/T^4]
  // * sqrt[45/(4 pi^3 g_star(T))] * [g_star^now/g_star(T)] * (T_now/H_now)^3
  // and we take g_star(T) to be 105.75 (what it is for the SM above
  // temperatures of m_top) and conservatively take it as constant from T = 0
  // to T_dom. Hence we have
  // exp( -1.581E+106 GeV * integral of T^(-2) exp(-S_3(T)/T) dT ) )
  // integrated from T = 0 to T_dom (as the contribution from higher
  // temperatures drops off very quickly).
  // 1.581E+106 is exp( 244.53 = lnOfThermalIntegrationFactor ) which is
  // in agreement with the value of 240 quoted in the CosmoTransitions manual
  // for an estimate of the threshold S_3(T)/T for T= 100 GeV.
  // Kusenko (and others in the literature, including Wainwright implicitly in
  // the CosmoTransitions manual as just mentioned) took the integral of
  // exp(-S_3(T)/T) T^(-2) to be exp( S_3(T_dom)/T_dom) T_dom^(-1) where T_dom
  // is the optimal tunneling temperature which dominates the integral.
  // This might be a bit aggressive, and taking S_3(T) to be approximated by
  // S_3(0) + T S' leads to the integral being
  // exp( -S_3(T_dom)/T_dom ) / S_3(0)
  // Assuming that S_3(0) < S_3(T_dom) (which should hold for all cases of
  // interest), the full integral should be between
  // exp( -S_3(T_dom)/T_dom ) / S_3(T_dom)
  // and exp( -S_3(T_dom)/T_dom ) / T_dom
  // For a threshold survival probability P,
  // 1.581E+106 GeV * integral of T^(-2) exp(-S_3(T)/T) dT )
  // should be larger than ln(1/P). Hence we compare
  // (S_3(T_dom)/T_dom) + ln( x / GeV ) to
  // lnOfThermalIntegrationFactor - ln( ln(1/P) ) where x is either
  // S_3(T_dom) or T_dom.


  BounceActionTunneler::BounceActionTunneler(
                                          PotentialFunction& potentialFunction,
                TunnelingCalculator::TunnelingStrategy const tunnelingStrategy,
                                     double const survivalProbabilityThreshold,
                                              size_t const temperatureAccuracy,
                                         size_t const evaporationResolution ) :
    TunnelingCalculator( potentialFunction,
                         tunnelingStrategy,
                         survivalProbabilityThreshold ),
    potentialFunction( potentialFunction ),
    temperatureAccuracy( temperatureAccuracy ),
    evaporationResolution( evaporationResolution ),
    thermalPotentialMinimizer( potentialFunction ),
    evaporationMinimum(),
    criticalMinimum(),
    criticalRatherThanEvaporation( true )
  {
    // This constructor is just an initialization list.
  }

  BounceActionTunneler::~BounceActionTunneler()
  {
    // This does nothing.
  }


  // This decides what virtual tunneling calculation functions to call based
  // on tunnelingStrategy.
  void BounceActionTunneler::CalculateTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    // First we set all variables to their "not calculated" values.
    quantumSurvivalProbability = -1.0;
    quantumLifetimeInSeconds = -1.0;
    thermalSurvivalProbability = -1.0;
    dominantTemperatureInGigaElectronVolts = -1.0;

    if( tunnelingStrategy == NoTunneling )
    {
      std::cout
      << std::endl
      << "Not tunneling as tunneling strategy is \"NoTunneling\"";
      std::cout << std::endl;

      return;
    }
    PrepareCommonExtras();
    if( tunnelingStrategy == JustQuantum )
    {
      CalculateQuantumTunneling( falseVacuum,
                                 trueVacuum );
    }
    else if( tunnelingStrategy == JustThermal )
    {
      CalculateThermalTunneling( falseVacuum,
                                 trueVacuum );
    }
    else if( tunnelingStrategy == QuantumThenThermal )
    {
      CalculateQuantumTunneling( falseVacuum,
                                 trueVacuum );
      if( quantumSurvivalProbability > survivalProbabilityThreshold )
      {
        CalculateThermalTunneling( falseVacuum,
                                   trueVacuum );
      }
    }
    else if( tunnelingStrategy == ThermalThenQuantum )
    {
      CalculateThermalTunneling( falseVacuum,
                                 trueVacuum );
      if( thermalSurvivalProbability > survivalProbabilityThreshold )
      {
        CalculateQuantumTunneling( falseVacuum,
                                   trueVacuum );
      }
    }
    else
    {
      std::cout
      << std::endl
      << "No valid tunneling strategy was set, so treating it as"
      << " \"NoTunneling\"!";
      std::cout << std::endl;
    }
  }

  // This sets quantumSurvivalProbability and quantumLifetimeInSeconds
  // appropriately.
  void BounceActionTunneler::CalculateQuantumTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    double quantumAction( BounceAction( falseVacuum,
                                        trueVacuum,
                                        0.0 ) );
    double const fourthRootOfSolitonicFactor( sqrt(
                potentialFunction.ScaleSquaredRelevantToTunneling( falseVacuum,
                                                              trueVacuum ) ) );
    logOfMinusLogOfQuantumProbability
    = ( log( ( ageOfKnownUniverseInSeconds * hBarInGigaElectronVoltSeconds )
             / fourthRootOfSolitonicFactor ) - ( 0.25 * quantumAction ) );
    if( quantumAction >= maximumPowerOfNaturalExponent )
    {
      quantumLifetimeInSeconds = 1.0;
      quantumSurvivalProbability = 0.0;
      std::cout
      << std::endl
      << "Warning! The calculated bounce action was so large and positive that"
      << " exponentiating it would result in an overflow error, so capping the"
      << " lifetime at " << quantumLifetimeInSeconds
      << " and the survival probability at " << quantumSurvivalProbability;
      std::cout << std::endl;
      return;
    }
    else if( quantumAction <= -maximumPowerOfNaturalExponent )
    {
      quantumLifetimeInSeconds = 1.0E+100;
      quantumSurvivalProbability = 1.0;
      std::cout
      << std::endl
      << "Warning! The calculated bounce action was so large and negative that"
      << " exponentiating it would result in an overflow error, so capping the"
      << " lifetime at " << quantumLifetimeInSeconds
      << " and the survival probability at " << quantumSurvivalProbability;
      std::cout << std::endl;
      return;
    }
    quantumLifetimeInSeconds
    = ( ( exp( 0.25 * quantumAction )
          * hBarInGigaElectronVoltSeconds )
        / fourthRootOfSolitonicFactor );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "quantumLifetimeInSeconds = " << quantumLifetimeInSeconds;
    std::cout << std::endl;/**/

    double
    survivalExponent( ageOfKnownUniverseInSeconds / quantumLifetimeInSeconds );
    if( survivalExponent >= maximumPowerOfNaturalExponent )
    {
      std::cout
      << std::endl
      << "Warning! The calculated decay width was so large that"
      << " exponentiating it would result in an overflow error, so setting the"
      << " survival probability to zero.";
      quantumSurvivalProbability = 0.0;
    }
    else
    {
      quantumSurvivalProbability = exp( -survivalExponent );
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "quantumSurvivalProbability = " << quantumSurvivalProbability;
    std::cout << std::endl;/**/
  }

  // This calculates the temperature at which either tunneling from
  // givenVacuum to the field origin becomes impossible if
  // criticalRatherThanEvaporation is true or the temperature at which
  // givenVacuum evaporates if false.
  double BounceActionTunneler::CriticalOrEvaporationTemperature(
                                               double const potentialAtOrigin )
  {
    PotentialMinimum* minimumPointer( NULL );
    std::cout
    << std::endl
    << "Looking for temperature ";
    if( criticalRatherThanEvaporation )
    {
      std::cout <<  "at which tunneling from the field origin to "
      << criticalMinimum.AsMathematica( potentialFunction.FieldNames() )
      << " becomes impossible.";
      minimumPointer = &criticalMinimum;
    }
    else
    {
      std::cout <<  "at which "
      << evaporationMinimum.AsMathematica( potentialFunction.FieldNames() )
      << " evaporates.";
      minimumPointer = &evaporationMinimum;
    }
    std::cout << std::endl;

    // The corrections are ( T^4 / ( 2 pi^2 ) ) * sum of J functions, and the
    // values of the J functions are about 2 for massless bosonic & fermionic
    // degrees of freedom, & there are ~100 degrees of freedom in the SM. Hence
    // we take the coefficient of T^4 to be 100 / ( 2 pi^2 ) ~= 5.
    double temperatureGuess( pow( ( 0.2 * ( potentialAtOrigin
               - potentialFunction( minimumPointer->FieldConfiguration() ) ) ),
                                              0.25 ) );
    // We aim to have a pair of temperatures, one above the sought temperature,
    // the other below. If the initial guess was below the sought temperature,
    // we start doubling the temperature, recording the previous temperature
    // each time. If it was above, we start halving the temperature, recording
    // the previous temperature each time.
    std::cout << "Trying " << temperatureGuess << " GeV.";
    std::cout << std::endl;

    while( BelowCriticalOrEvaporation( temperatureGuess ) )
    {
      temperatureGuess += temperatureGuess;
      std::cout << "... too low. Trying " << temperatureGuess << " GeV.";
      std::cout << std::endl;
    }
    // Now temperatureGuess is definitely about the sought temperature, so we
    // halve temperatureGuess and see if it is still too high, & if so, keep
    // halving.
    temperatureGuess = ( 0.5 * temperatureGuess );
    while( !(BelowCriticalOrEvaporation( temperatureGuess )) )
    {
      temperatureGuess = ( 0.5 * temperatureGuess );
      std::cout << "... too high. Trying " << temperatureGuess << " GeV.";
      std::cout << std::endl;
    }
    // At this point, temperatureGuess should be between 0.5 and 1.0 times the
    // value which we should return.
    double justBelowTemperature( temperatureGuess );
    double justAboveTemperature( temperatureGuess + temperatureGuess );
    // We aim to be within a factor of 2^( -temperatureAccuracy ) of the
    // critical temperature,
    // hence temperatureAccuracy iterations of this loop.
    for( size_t narrowingStep( 0 );
         narrowingStep < temperatureAccuracy;
         ++narrowingStep )
    {
      temperatureGuess = sqrt( justBelowTemperature * justAboveTemperature );
      std::cout << "Trying " << temperatureGuess << " GeV.";
      std::cout << std::endl;
      if( BelowCriticalOrEvaporation( temperatureGuess ) )
      {
        justBelowTemperature = temperatureGuess;
      }
      else
      {
        justAboveTemperature = temperatureGuess;
      }
    }
    std::cout << "Temperature lies between " << justBelowTemperature
    << " GeV and " << justAboveTemperature << " GeV.";
    std::cout << std::endl;
    return sqrt( justBelowTemperature * justAboveTemperature );
  }

  // This returns true if the temperature is below that at which
  // evaporationMinimum is no longer separated from the field origin by an
  // energy barrier.
  bool BounceActionTunneler::BelowEvaporationTemperature(
                                                double const temperatureGuess )
  {
    double const stepFraction( 1.0
                            / static_cast< double >( evaporationResolution ) );
    std::vector< double >
    fieldConfiguration( potentialFunction.FieldValuesOrigin() );
    double lastPotentialValue( potentialFunction( fieldConfiguration,
                                                  temperatureGuess ) );
    double currentPotentialValue( lastPotentialValue );
    // We compare the potential at a series of points from the field origin
    // to evaporationMinimum. If the potential ever goes down on the way, there
    // is an energy barrier, so we note that we are below the evaporation
    // temperature.
    for( size_t whichStep( 1 );
         whichStep <= evaporationResolution;
         ++whichStep )
    {
      for( size_t fieldIndex( 0 );
           fieldIndex < fieldConfiguration.size();
           ++fieldIndex )
      {
        fieldConfiguration[ fieldIndex ]
        = ( whichStep * stepFraction
            * evaporationMinimum.FieldConfiguration()[ fieldIndex ] );
      }
      currentPotentialValue = potentialFunction( fieldConfiguration,
                                                 temperatureGuess );

      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "step " << whichStep << ", fieldConfiguration = "
      << potentialFunction.FieldConfigurationAsMathematica(
                                                           fieldConfiguration )
      << ", lastPotentialValue = " << lastPotentialValue
      << ", currentPotentialValue = " << currentPotentialValue;
      std::cout << std::endl;*/

      if( currentPotentialValue < lastPotentialValue )
      {
        return true;
      }
      lastPotentialValue = currentPotentialValue;
    }
    // If the loop ends without returning, there was no energy barrier.
    return false;
  }

  // This ensures that thermalSurvivalProbability is set correctly from
  // logOfMinusLogOfThermalProbability.
  void BounceActionTunneler::SetThermalSurvivalProbability()
  {
    if( logOfMinusLogOfThermalProbability >= maximumPowerOfNaturalExponent )
    {
      std::cout
      << std::endl
      << "Warning! The calculated bounce action was so large and positive that"
      << " exponentiating it would result in an overflow error, so setting the"
      << " survival probability to 0.";
      std::cout << std::endl;
      thermalSurvivalProbability = 0.0;
    }
    else if( logOfMinusLogOfThermalProbability
             <= -maximumPowerOfNaturalExponent )
    {
      std::cout
      << std::endl
      << "Warning! The calculated bounce action was so large and negative that"
      << " exponentiating it would result in an overflow error, so setting the"
      << " survival probability to 1.";
      std::cout << std::endl;
      thermalSurvivalProbability = 1.0;
    }
    else if( exp( logOfMinusLogOfThermalProbability )
             >= maximumPowerOfNaturalExponent )
    {
      std::cout
      << std::endl
      << "Warning! The calculated integrated decay width was so large and"
      << " positive that exponentiating it would result in an overflow error,"
      << " so setting the survival probability to 0.";
      std::cout << std::endl;
      thermalSurvivalProbability = 0.0;
    }
    else
    {
      thermalSurvivalProbability
      = exp( -exp( logOfMinusLogOfThermalProbability ) );
    }
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "thermalSurvivalProbability = " << thermalSurvivalProbability;
    std::cout << std::endl;/**/
  }

} /* namespace VevaciousPlusPlus */
