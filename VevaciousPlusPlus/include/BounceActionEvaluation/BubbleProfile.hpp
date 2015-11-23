/*
 * BubbleProfile.hpp
 *
 *  Created on: May 19, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BUBBLEPROFILE_HPP_
#define BUBBLEPROFILE_HPP_

#include "CommonIncludes.hpp"
#include "boost/numeric/odeint/integrate/integrate.hpp"
#include "boost/math/special_functions/bessel.hpp"
#include "BasicFunctions/SimplePolynomial.hpp"
#include "SplinePotential.hpp"
#include "PathParameterization/TunnelPath.hpp"
#include "OdeintBubbleDerivatives.hpp"
#include "OdeintBubbleObserver.hpp"

namespace VevaciousPlusPlus
{

  class BubbleProfile
  {
  public:
    BubbleProfile();
    virtual ~BubbleProfile();


    // The bounce action for the bubble profile is kept along with it, to save
    // having to constantly pass 2 arguments when improving the tunneling path.
    double BounceAction() const { return bounceAction; }
    double& BounceAction() { return bounceAction; }


    // This should set up the bubble profile in terms of the auxiliary variable
    // and its slope, at values of the radial variable.
    virtual void CalculateProfile( TunnelPath const& tunnelPath,
                   OneDimensionalPotentialAlongPath const& pathPotential ) = 0;

    // This should return the value that the auxiliary variable should have for
    // the radial value given by radialValue.
    virtual double AuxiliaryAt( double const radialValue ) const = 0;


    // This should return the derivative of the path auxiliary value with
    // respect to the radial variable, evaluated at a value of radialValue for
    // the radial variable.
    virtual double AuxiliarySlopeAt( double const radialValue ) const = 0;

    // This should return the maximum radial value up to which a profile can be
    // reasonably plotted.
    virtual double MaximumPlotRadius() const = 0;


  protected:
    double bounceAction;
  };

} /* namespace VevaciousPlusPlus */
#endif /* BUBBLEPROFILE_HPP_ */
