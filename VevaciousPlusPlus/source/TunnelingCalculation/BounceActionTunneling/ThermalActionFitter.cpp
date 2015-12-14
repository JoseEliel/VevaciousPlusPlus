/*
 * ThermalActionFitter.cpp
 *
 *  Created on: Apr 28, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "TunnelingCalculation/BounceActionTunneling/ThermalActionFitter.hpp"

namespace VevaciousPlusPlus
{

  ThermalActionFitter::ThermalActionFitter(
                                  std::vector< double > const& fitTemperatures,
                                    std::vector< double > const& fittedActions,
                                           double const criticalTemperature ) :
    fitCoefficients( fittedActions.size() ),
    criticalTemperature( criticalTemperature )
  {
    std::vector::size_type const numberOfNodes( fittedActions.size() );
    if( numberOfNodes > 1 )
    {
      Eigen::MatrixXd fitMatrix( numberOfNodes,
                                 numberOfNodes );
      Eigen::VectorXd actionTimesTemperatureDifferenceVector( numberOfNodes );
      double currentTemperature( 0.0 );
      double temperatureDifference( 0.0 );
      for( std::vector::size_type whichNode( 0 );
           whichNode < numberOfNodes;
           ++whichNode )
      {
        currentTemperature = fitTemperatures[ whichNode ];
        temperatureDifference = ( criticalTemperature - currentTemperature );
        actionTimesTemperatureDifferenceVector( whichNode )
         = ( fittedActions[ whichNode ]
             * temperatureDifference * temperatureDifference );
        fitMatrix( whichNode,
                   0 ) = 1.0;
        fitMatrix( whichNode,
                   1 ) = currentTemperature;
        for( unsigned int whichPower( 2 );
             whichPower < numberOfNodes;
             ++whichPower )
        {
          fitMatrix( whichNode,
                     whichPower ) = pow( currentTemperature,
                                         whichPower );
        }
      }
      // Now we solve for the coefficients:
      Eigen::VectorXd
      solutionVector( fitMatrix.colPivHouseholderQr().solve(
                                    actionTimesTemperatureDifferenceVector ) );
      for( std::vector::size_type whichNode( 0 );
           whichNode < numberOfNodes;
           ++whichNode )
      {
        fitCoefficients[ whichNode ] = solutionVector( whichNode );
      }
    }
    else
    {
      double const
      temperatureDifference( criticalTemperature - fitTemperatures[ 0 ] );
      fitCoefficients[ 0 ] = ( fittedActions[ 0 ]
                               * temperatureDifference
                               * temperatureDifference );
    }
  }

  ThermalActionFitter::~ThermalActionFitter()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
