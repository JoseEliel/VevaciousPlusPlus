/*
 * ThermalActionFitter.hpp
 *
 *  Created on: Apr 28, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef THERMALACTIONFITTER_HPP_
#define THERMALACTIONFITTER_HPP_

#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include <vector>
#include <cmath>
#include <cstddef>
#include "Eigen/Dense"

namespace VevaciousPlusPlus
{
  // This class is just for fitting the thermal action from a set of actions
  // with temperatures to a function of the form of a polynomial divided by
  // the square of the difference of the temperature from the critical
  // temperature.
  class ThermalActionFitter : public ROOT::Minuit2::FCNBase
  {
  public:
    ThermalActionFitter( std::vector< double > const& fitTemperatures,
                         std::vector< double > const& fittedActions,
                         double const criticalTemperature );
    virtual ~ThermalActionFitter();


    // This implements operator() for FCNBase, the function that MINUIT will
    // minimize. It actually returns ( ( S / T ) + ln( S ) ) for a fitted
    // action S for the temperature T which should be the only element of
    // temperatureVector.
    virtual double
    operator()( std::vector< double > const& temperatureVector ) const;

    // This implements Up() for FCNBase just to stick to a basic value.
    virtual double Up() const { return 1.0; }


  protected:
    std::vector< double > fitCoefficients;
    double const criticalTemperature;
  };




  // This implements operator() for FCNBase, the function that MINUIT will
  // minimize. It actually returns ( ( S / T ) + ln( S ) ) for a fitted
  // action S for the temperature T which should be the only element of
  // temperatureVector.
  inline double ThermalActionFitter::operator()(
                         std::vector< double > const& temperatureVector ) const
  {
    double const
    temperatureDifference( criticalTemperature - temperatureVector[ 0 ] );
    double fittedAction( fitCoefficients[ 0 ] );
    double termContribution( 0.0 );
    for( size_t whichPower( 1 );
         whichPower < fitCoefficients.size();
         ++whichPower )
    {
      termContribution = fitCoefficients[ whichPower ];
      for( size_t powerCount( 0 );
           powerCount < whichPower;
           ++powerCount )
      {
        termContribution *= temperatureVector[ 0 ];
      }
      fittedAction += termContribution;
    }
    fittedAction
    = ( fittedAction / ( temperatureDifference * temperatureDifference ) );
    return ( ( fittedAction / temperatureVector[ 0 ] ) + log( fittedAction ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* THERMALACTIONFITTER_HPP_ */
