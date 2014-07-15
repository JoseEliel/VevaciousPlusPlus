/*
 * BubbleRadialValueDescription.hpp
 *
 *  Created on: May 22, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BUBBLERADIALVALUEDESCRIPTION_HPP_
#define BUBBLERADIALVALUEDESCRIPTION_HPP_

namespace VevaciousPlusPlus
{

  class BubbleRadialValueDescription
  {
  public:
    static bool
    SortByRadialValue( BubbleRadialValueDescription const& firstReference,
                       BubbleRadialValueDescription const& secondReference )
    { return ( firstReference.radialValue < secondReference.radialValue ); }

    BubbleRadialValueDescription();
    BubbleRadialValueDescription( double const radialValue,
                                  double const auxiliaryValue,
                                  double const auxiliarySlope );
    BubbleRadialValueDescription(
                              BubbleRadialValueDescription const& copySource );
    virtual ~BubbleRadialValueDescription();

    double radialValue;
    double auxiliaryValue;
    double auxiliarySlope;
  };

} /* namespace VevaciousPlusPlus */
#endif /* BUBBLERADIALVALUEDESCRIPTION_HPP_ */
