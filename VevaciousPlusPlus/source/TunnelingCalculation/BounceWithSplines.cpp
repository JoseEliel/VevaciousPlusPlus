/*
 * BounceWithSplines.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{
  double const BounceWithSplines::maximumPowerOfNaturalExponent(
                           log( 0.5 * std::numeric_limits< double >::max() ) );
  double const
  BounceWithSplines::hBarInGigaElectronVoltSeconds( 6.58211928E-25 );
  double const BounceWithSplines::ageOfKnownUniverseInSeconds( 4.3E+17 );
  double const BounceWithSplines::lnOfThermalIntegrationFactor( 238.553 );
// Based on correspondence with Alexander Kusenko and discussion with Bjoern
// Garbrecht:
// Taking [decay width per horizon]
// = [horizon volume] * [solitonic coefficient] * exp(-[thermal action]/T)
// at temperature T, where [horizon volume] = ( M_Plank / T^2 )^3, & taking
// [solitonic coefficient] to be T^4,
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
// where C = M_Planck * [solitonic coefficient/T^4]
// * sqrt[45/(4 pi^2 g_star(T))] * [g_star^now/g_star(T)] * (T_now/H_now)^3
// and we take g_star(T) to be 105.75 (what it is for the SM above
// temperatures of m_top) and conservatively take it as constant from T = 0
// to T_opt. Hence we have
// exp( -4.00163E+103 GeV * integral of T^(-2) exp(-S_3(T)/T) dT ) )
// integrated from T = 0 to T_opt (as the contribution from higher
// temperatures drops off very quickly).
// 4.00163E+103 is exp( 238.553 = lnOfThermalIntegrationFactor ) which is in
// agreement with the value of 240 quoted in the CosmoTransitions manual for
// an estimate of the threshold S_3(T)/T for T= 100 GeV.
// Kusenko (and others in the literature, including Wainwright implicitly in
// the CosmoTransitions manual as just mentioned) took the integral of
// exp(-S_3(T)/T) T^(-2) to be exp( S_3(T_opt)/T_opt) T_opt^(-1) where T_opt
// is the optimal tunneling temperature which dominates the integral.
// This might be a bit aggressive, and taking S_3(T) to be approximated by
// S_3(0) + T S' leads to the integral being
// exp( -S_3(T_opt)/T_opt ) / S_3(0)
// Assuming that S_3(0) < S_3(T_opt) (which should hold for all cases of
// interest), the full integral should be between
// exp( -S_3(T_opt)/T_opt ) / S_3(T_opt)
// and exp( -S_3(T_opt)/T_opt ) / T_opt
// For a threshold survival probability P,
// 4.00163E+103 GeV * integral of T^(-2) exp(-S_3(T)/T) dT )
// should be larger than ln(1/P). Hence we compare
// (S_3(T_opt)/T_opt) + ln( x / GeV ) to
// lnOfThermalIntegrationFactor - ln( ln(1/P) ) where x is either
// S_3(T_opt) or T_opt.


  BounceWithSplines::BounceWithSplines( PotentialFunction& potentialFunction,
                                     TunnelingStrategy const tunnelingStrategy,
                                  double const survivalProbabilityThreshold ) :
    TunnelingCalculator( potentialFunction,
                         tunnelingStrategy,
                         survivalProbabilityThreshold ),
    potentialFunction( potentialFunction )
  {
    // This constructor is just an initialization list.
  }

  BounceWithSplines::~BounceWithSplines()
  {
    // This does nothing.
  }


  // This does something.
  void BounceWithSplines::CalculateTunneling(
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

} /* namespace VevaciousPlusPlus */
