/*
 * MinuitBounceActionMinimizer.cpp
 *
 *  Created on: Apr 15, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/BounceActionTunneling/MinuitBounceActionMinimizer.hpp"

namespace VevaciousPlusPlus
{

  MinuitBounceActionMinimizer::MinuitBounceActionMinimizer(
                                          PotentialFunction& potentialFunction,
                                     TunnelingStrategy const tunnelingStrategy,
                                     double const survivalProbabilityThreshold,
                                                    size_t const numberOfNodes,
                                                  double const initialStepSize,
                                             unsigned int const minuitStrategy,
                                               double const minuitTolerance ) :
    BounceWithSplines( potentialFunction,
                       tunnelingStrategy,
                       survivalProbabilityThreshold ),
    numberOfNodes( numberOfNodes ),
    initialStepSize( initialStepSize ),
    minuitStrategy( minuitStrategy ),
    minuitTolerance( minuitTolerance )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::MinuitBounceActionMinimizer( ... )";
    std::cout << std::endl;/**/
  }

  MinuitBounceActionMinimizer::~MinuitBounceActionMinimizer()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::~MinuitBounceActionMinimizer()";
    std::cout << std::endl;/**/
  }


  // This doesn't do anything here.
  void MinuitBounceActionMinimizer::UpdateSelfForNewSlha(
                                               SlhaManager const& slhaManager )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::UpdateSelfForNewSlha(...)";
    std::cout << std::endl;/**/
  }

  // This returns either the dimensionless bounce action integrated over
  // four dimensions (for zero temperature) or the dimensionful bounce
  // action integrated over three dimensions (for non-zero temperature) for
  // tunneling from falseVacuum to trueVacuum at temperature
  // tunnelingTemperature.
  double MinuitBounceActionMinimizer::BounceAction(
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                      double const tunnelingTemperature ) const
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "weird path check:";
    std::cout << std::endl;
    std::vector< double > fieldConfiguration( 4 );
    fieldConfiguration[ 0 ] = -67.1847;
    fieldConfiguration[ 1 ] = 216.362;
    fieldConfiguration[ 2 ] = 216.675;
    fieldConfiguration[ 3 ] = -35.2998;
    std::cout << "p = 0.01: field configuration = { ";
    for( std::vector< double >::const_iterator
         fieldValue( fieldConfiguration.begin() );
         fieldValue < fieldConfiguration.end();
         ++fieldValue )
    {
      if( fieldValue > fieldConfiguration.begin() )
      {
        std::cout << ", ";
      }
      std::cout <<  *fieldValue;
    }
    std::cout << " } has potential = "
    << potentialFunction( fieldConfiguration ) << std::endl;
    fieldConfiguration[ 0 ] = -50.0763;
    fieldConfiguration[ 1 ] = 236.54;
    fieldConfiguration[ 2 ] = 245.357;
    fieldConfiguration[ 3 ] = -15.7476;
    std::cout << "p = 0.02: field configuration = { ";
    for( std::vector< double >::const_iterator
         fieldValue( fieldConfiguration.begin() );
         fieldValue < fieldConfiguration.end();
         ++fieldValue )
    {
      if( fieldValue > fieldConfiguration.begin() )
      {
        std::cout << ", ";
      }
      std::cout <<  *fieldValue;
    }
    std::cout << " } has potential = "
    << potentialFunction( fieldConfiguration ) << std::endl;
    fieldConfiguration[ 0 ] = -1.17451;
    fieldConfiguration[ 1 ] = 269.903;
    fieldConfiguration[ 2 ] = 207.585;
    fieldConfiguration[ 3 ] = 23.2112;
    std::cout << "p = 0.03: field configuration = { ";
    for( std::vector< double >::const_iterator
         fieldValue( fieldConfiguration.begin() );
         fieldValue < fieldConfiguration.end();
         ++fieldValue )
    {
      if( fieldValue > fieldConfiguration.begin() )
      {
        std::cout << ", ";
      }
      std::cout <<  *fieldValue;
    }
    std::cout << " } has potential = "
    << potentialFunction( fieldConfiguration ) << std::endl;
    fieldConfiguration[ 0 ] = 51.6088;
    fieldConfiguration[ 1 ] = 304.704;
    fieldConfiguration[ 2 ] = 162.002;
    fieldConfiguration[ 3 ] = 64.4757;
    std::cout << "p = 0.04: field configuration = { ";
    for( std::vector< double >::const_iterator
         fieldValue( fieldConfiguration.begin() );
         fieldValue < fieldConfiguration.end();
         ++fieldValue )
    {
      if( fieldValue > fieldConfiguration.begin() )
      {
        std::cout << ", ";
      }
      std::cout <<  *fieldValue;
    }
    std::cout << " } has potential = "
    << potentialFunction( fieldConfiguration ) << std::endl;
    fieldConfiguration[ 0 ] = 97.1669;
    fieldConfiguration[ 1 ] = 336.344;
    fieldConfiguration[ 2 ] = 131.835;
    fieldConfiguration[ 3 ] = 101.27;
    std::cout << "p = 0.05: field configuration = { ";
    for( std::vector< double >::const_iterator
         fieldValue( fieldConfiguration.begin() );
         fieldValue < fieldConfiguration.end();
         ++fieldValue )
    {
      if( fieldValue > fieldConfiguration.begin() )
      {
        std::cout << ", ";
      }
      std::cout <<  *fieldValue;
    }
    std::cout << " } has potential = "
    << potentialFunction( fieldConfiguration ) << std::endl;
    fieldConfiguration[ 0 ] = 133.094;
    fieldConfiguration[ 1 ] = 363.89;
    fieldConfiguration[ 2 ] = 122.015;
    fieldConfiguration[ 3 ] = 132.152;
    std::cout << "p = 0.06: field configuration = { ";
    for( std::vector< double >::const_iterator
         fieldValue( fieldConfiguration.begin() );
         fieldValue < fieldConfiguration.end();
         ++fieldValue )
    {
      if( fieldValue > fieldConfiguration.begin() )
      {
        std::cout << ", ";
      }
      std::cout <<  *fieldValue;
    }
    std::cout << " } has potential = "
    << potentialFunction( fieldConfiguration ) << std::endl;
    std::cout << std::endl;/**/
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::BounceAction(...)";
    std::cout << std::endl;/**/

    ModifiedBounceForMinuit modifiedBounceForMinuit( potentialFunction,
                                                     numberOfNodes,
                                                     16,
                                                     falseVacuum,
                                                     trueVacuum,
                                                     tunnelingTemperature,
                                                     32 );
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "Remember to write in taking shoot attempts, path resolution, MINUIT"
    << " strategy, MINUIT tolerance, etc.,"
    << " from input XML.";
    std::cout << std::endl;/**/
    std::vector< double > pathParameterization;
    std::vector< double > initialStepSizes;
    modifiedBounceForMinuit.SetUpStraightPathForMinuit( pathParameterization,
                                                        initialStepSizes,
                                                        tunnelingTemperature,
                                                        initialStepSize );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "straight path bounce = "
    << modifiedBounceForMinuit( pathParameterization );
    std::cout << std::endl;/**/
    std::vector< std::string > plotColors;
    plotColors.push_back( "red" );
    plotColors.push_back( "purple" );
    plotColors.push_back( "blue" );
    plotColors.push_back( "green" );
    plotColors.push_back( "orange" );
    plotColors.push_back( "gold" );
    modifiedBounceForMinuit.PlotBubbleProfile( pathParameterization,
                                               plotColors,
                                               "StraightBubbleProfile.eps" );
    // more debugging:
    /**/if( !(initialStepSizes.empty()) )
    {
      for( size_t parameterIndex( 0 );
           parameterIndex < pathParameterization.size();
           ++parameterIndex )
      {
        pathParameterization[ parameterIndex ]
        = ( initialStepSizes[ parameterIndex ]
            * (double)( 1
                        + ( parameterIndex
                    / ( potentialFunction.NumberOfFieldVariables() - 1 ) ) ) );
      }
      std::cout << "test path bounce = "
      << modifiedBounceForMinuit( pathParameterization );
      modifiedBounceForMinuit.PlotBubbleProfile( pathParameterization,
                                                 plotColors,
                                                 "WigglyBubbleProfile.eps" );
      for( size_t parameterIndex( 0 );
           parameterIndex < pathParameterization.size();
           ++parameterIndex )
      {
        pathParameterization[ parameterIndex ] = 0.0;
      }
    }
    std::cout << std::endl;/**/

    // Should now do migrad...

    ROOT::Minuit2::MnMigrad mnMigrad( modifiedBounceForMinuit,
                                      pathParameterization,
                                      initialStepSizes,
                                      minuitStrategy );

    size_t const numberOfMinuitVariables( pathParameterization.size() );
    size_t const minuitStepsPerThresholdCheck( 100 * numberOfMinuitVariables );

    ROOT::Minuit2::FunctionMinimum
    minuitMinimum( mnMigrad( minuitStepsPerThresholdCheck,
                             100.0 ) );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "after 1st migrad( " << minuitStepsPerThresholdCheck << ", "
    << minuitTolerance << "), minuitMinimum =" << std::endl << minuitMinimum;
    std::cout << std::endl;/**/

    double thresholdBounceAction( NAN );
    if( tunnelingTemperature > 0.0 )
    {
      // If we're calculating thermal tunneling, we start with the usual
      // estimate for the threshold (arXiv:1405.7376) assuming that the
      // logarithm vanishes, then adjust until
      // thresholdBounceAction / T + ln( thresholdBounceAction )
      // <= lnOfThermalIntegrationFactor (= 244.53).
      // Actually, we lazily drop the factor of T until the end.
      thresholdBounceAction = lnOfThermalIntegrationFactor;
      double const thresholdComparison( lnOfThermalIntegrationFactor
                                        - log( tunnelingTemperature ) );
      while( ( thresholdBounceAction + log( thresholdBounceAction ) )
             >= thresholdComparison )
      {
        thresholdBounceAction -= 1.0;
      }
      thresholdBounceAction *= tunnelingTemperature;
    }
    else
    {
      double const squareOfSolitonicFactor(
                potentialFunction.ScaleSquaredRelevantToTunneling( falseVacuum,
                                                                trueVacuum ) );
      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "-log( survivalProbabilityThreshold ) = "
      << -log( survivalProbabilityThreshold )
      << ", fourVolumeOfKnownUniverseOverGevFourth = "
      << fourVolumeOfKnownUniverseOverGevFourth
      << ", squareOfSolitonicFactor = " << squareOfSolitonicFactor
      << ", ( fourVolumeOfKnownUniverseOverGevFourth * squareOfSolitonicFactor"
      << " * squareOfSolitonicFactor ) = "
      << ( fourVolumeOfKnownUniverseOverGevFourth * squareOfSolitonicFactor
                                                  * squareOfSolitonicFactor )
      << ", -ln(P) / (A*4-vol) = "
      << ( -log( survivalProbabilityThreshold )
          / ( fourVolumeOfKnownUniverseOverGevFourth
              * squareOfSolitonicFactor
              * squareOfSolitonicFactor ) );
      std::cout << std::endl;*/
      thresholdBounceAction = -log( -log( survivalProbabilityThreshold )
                                    / ( fourVolumeOfKnownUniverseOverGevFourth
                                        * squareOfSolitonicFactor
                                        * squareOfSolitonicFactor ) );
    }

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "thresholdBounceAction = " << thresholdBounceAction;
    std::cout << std::endl;/**/

    // We keep looping while MINUIT reports that it has not yet found a
    // satisfactory minimum within the last minuitStepsPerThresholdCheck steps.
    // If minuitMinimum.IsValid() is false for any other reason than this, we
    // stop looping as well. If we find that the (modified) bounce action has
    // already dropped below the threshold even though MINUIT is not at a
    // minimum, we still break out of the loop.
    while( !(minuitMinimum.IsValid())
           &&
           minuitMinimum.HasReachedCallLimit() )
    {
      if( minuitMinimum.Fval() < thresholdBounceAction )
      {
        std::cout
        << std::endl
        << "Bounce action below threshold before MINUIT finished, ending"
        << " migrad() early.";
        std::cout << std::endl;

        break;
      }
      minuitMinimum = mnMigrad( minuitStepsPerThresholdCheck,
                                minuitTolerance );
      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "after loop migrad( " << minuitStepsPerThresholdCheck << ", "
      << minuitTolerance << "), minuitMinimum =" << std::endl << minuitMinimum;
      std::cout << std::endl;/**/
    }
    MinuitMinimum minuitResult( pathParameterization.size(),
                                minuitMinimum );

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "minimized path bounce = "
    << minuitResult.FunctionValue();
    modifiedBounceForMinuit.PlotBubbleProfile( minuitResult.VariableValues(),
                                               plotColors,
                                               "MinuitBubbleProfile.eps" );
    std::cout << std::endl;/**/
    return minuitResult.FunctionValue();
  }

  // This should set thermalSurvivalProbability and
  // dominantTemperatureInGigaElectronVolts appropriately.
  void MinuitBounceActionMinimizer::CalculateThermalTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinuitBounceActionMinimizer::CalculateThermalTunneling(...)";
    std::cout << std::endl;/**/
  }

} /* namespace VevaciousPlusPlus */
