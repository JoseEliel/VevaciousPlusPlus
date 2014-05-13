/*
 * BubbleRadiusFromAuxiliary.hpp
 *
 *  Created on: May 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BUBBLERADIUSFROMAUXILIARY_HPP_
#define BUBBLERADIUSFROMAUXILIARY_HPP_

namespace VevaciousPlusPlus
{

  class BubbleRadiusFromAuxiliary
  {
  public:
    BubbleRadiusFromAuxiliary( double const appropriateScaleInGev );
    virtual
    ~BubbleRadiusFromAuxiliary();

    // This returns the radius corresponding to the auxiliary variable having
    // the value auxiliaryValue, in units of GeV^(-1).
    virtual double RadiusValue( double const auxiliaryValue ) const
    { return
      ( appropriateScaleInInverseGev * ( ( 1.0 / auxiliaryValue ) - 1.0 ) ); }

    // This returns the derivative of the radius with respect to auxiliary
    // variable, evaluated at the auxiliary variable having the value
    // auxiliaryValue, in units of GeV^(-1).
    virtual double RadiusDerivative( double const auxiliaryValue ) const
    { return ( -appropriateScaleInInverseGev
               / ( auxiliaryValue * auxiliaryValue ) ); }


  protected:
    double const appropriateScaleInInverseGev;
  };

} /* namespace VevaciousPlusPlus */
#endif /* BUBBLERADIUSFROMAUXILIARY_HPP_ */
