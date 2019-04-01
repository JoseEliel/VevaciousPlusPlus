/*
 * OdeintBubbleDerivatives.hpp
 *
 *  Created on: Jun 5, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef ODEINTBUBBLEDERIVATIVES_HPP_
#define ODEINTBUBBLEDERIVATIVES_HPP_

#include "OneDimensionalPotentialAlongPath.hpp"
#include "PathParameterization/TunnelPath.hpp"
#include <vector>

namespace VevaciousPlusPlus
{

  class OdeintBubbleDerivatives
  {
  public:
    OdeintBubbleDerivatives(
                         OneDimensionalPotentialAlongPath const& pathPotential,
                             TunnelPath const& tunnelPath ) :
      pathPotential( pathPotential ),
      tunnelPath( tunnelPath ),
      dampingFactor( tunnelPath.NonZeroTemperature() ? 2.0 : 3.0 ) {}

    virtual ~OdeintBubbleDerivatives() {}

    // This puts the first and second derivatives based on
    // auxiliaryAndFirstDerivative into firstAndSecondDerivatives, in the form
    // required for the Boost odeint package.
    void operator()( std::vector< double > const& auxiliaryAndFirstDerivative,
                     std::vector< double >& firstAndSecondDerivatives,
                     double const radialValue );


  protected:
    OneDimensionalPotentialAlongPath const& pathPotential;
    TunnelPath const& tunnelPath;
    double dampingFactor;
  };




  // This is in the form required for the Boost odeint package.
  inline void OdeintBubbleDerivatives::operator()(
                      std::vector< double > const& auxiliaryAndFirstDerivative,
                              std::vector< double >& firstAndSecondDerivatives,
                                                   double const radialValue )
  {
    double const auxiliaryValue( auxiliaryAndFirstDerivative[ 0 ] );
    // This cheats if there has already been an overshoot, to try to avoid the
    // integration going to small step sizes to resolve the oscillations of the
    // field going off to infinity.
    if( auxiliaryValue < 0.0 )
    {
      firstAndSecondDerivatives[ 0 ] = 0.0;
      firstAndSecondDerivatives[ 1 ] = 0.0;
      return;
    }
    double const auxiliaryDerivative( auxiliaryAndFirstDerivative[ 1 ] );
    firstAndSecondDerivatives[ 0 ] = auxiliaryDerivative;
    firstAndSecondDerivatives[ 1 ]
     = ( ( ( pathPotential.FirstDerivative( auxiliaryValue )
             - ( tunnelPath.SlopeDotAcceleration( auxiliaryValue )
                 * auxiliaryDerivative * auxiliaryDerivative ) )
           / tunnelPath.SlopeSquared( auxiliaryValue ) )
         - ( ( dampingFactor * auxiliaryDerivative ) / radialValue ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* ODEINTBUBBLEDERIVATIVES_HPP_ */
