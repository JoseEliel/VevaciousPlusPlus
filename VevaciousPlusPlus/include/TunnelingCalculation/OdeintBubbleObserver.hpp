/*
 * OdeintBubbleObserver.hpp
 *
 *  Created on: May 19, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef ODEINTBUBBLEOBSERVER_HPP_
#define ODEINTBUBBLEOBSERVER_HPP_

#include "../CommonIncludes.hpp"
#include "BubbleRadialValueDescription.hpp"

namespace VevaciousPlusPlus
{

  class OdeintBubbleObserver
  {
  public:
    OdeintBubbleObserver(
              std::vector< BubbleRadialValueDescription >& bubbleDescription );
    virtual
    ~OdeintBubbleObserver();


    // This sets needsOrdering to true if a radius is passed out of order, sets
    // definitelyUndershot to true if the auxiliary slope is positive or
    // definitelyOvershot to true if the auxiliary value is negative (the
    // auxiliary variable should decrease monotonically to zero for a perfect
    // shot), and pushes back the radial value, the auxiliary value, and its
    // slope into bubbleDescription.
    void operator()( std::vector< double > const& auxiliaryAndFirstDerivative,
                     double const radialValue );

    std::vector< BubbleRadialValueDescription > const&
    BubbleDescription() const{ return bubbleDescription; }

    // This returns definitelyUndershot which is set true if the slope of the
    // auxiliary variable goes positive for any radial value.
    bool DefinitelyUndershot() const{ return definitelyUndershot; }

    // This returns definitelyOvershot which is set true if the auxiliary
    // variable goes negative for any radial value.
    bool DefinitelyOvershot() const{ return definitelyOvershot; }

    // This is the lowest index of bubbleDescription where the auxiliary value
    // is negative.
    size_t OvershootIndex() const{ return overshootIndex; }

    // This gets set to true if a radius is inserted out of order.
    bool NeedsOrdering() const{ return needsOrdering; }

    // This resets everything for a new integration, including recording the
    // values before the first Euclidean step performed before the
    // odeint::integrate.
    void ResetValues( double const initialAuxiliary );

    // This sorts the vectors by increasing radial value if needsOrdering is
    // true.
    void SortByRadialValue();

    // This is for debugging.
    std::string AsDebuggingString() const;


  protected:
    std::vector< BubbleRadialValueDescription >& bubbleDescription;
    bool definitelyUndershot;
    bool definitelyOvershot;
    size_t overshootIndex;
    bool needsOrdering;
  };




  // This sets needsOrdering to true if a radius is passed out of order, sets
  // definitelyUndershot to true if the auxiliary slope is positive or
  // definitelyOvershot to true if the auxiliary value is negative (the
  // auxiliary variable should decrease monotonically to zero for a perfect
  // shot), and pushes back the radial value, the auxiliary value, and its
  // slope into bubbleDescription.
  inline void OdeintBubbleObserver::operator()(
                      std::vector< double > const& auxiliaryAndFirstDerivative,
                                                double const radialValue )
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "OdeintBubbleObserver::operator()( { "
    << auxiliaryAndFirstDerivative[ 0 ] << ", "
    << auxiliaryAndFirstDerivative[ 1 ] << " }, " << radialValue
    << " ) called. bubbleDescription.size() = " << bubbleDescription.size();
    std::cout << std::endl;/**/
    if( radialValue < bubbleDescription.back().radialValue )
    {
      needsOrdering = true;
    }
    bubbleDescription.push_back( BubbleRadialValueDescription( radialValue,
                                              auxiliaryAndFirstDerivative[ 0 ],
                                          auxiliaryAndFirstDerivative[ 1 ] ) );
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "bubbleDescription.size() = " << bubbleDescription.size()
    << std::endl
    << "bubbleDescription.back() = { r = "
    << bubbleDescription.back().radialValue << ", p = "
    << bubbleDescription.back().auxiliaryValue << ", dp/dr = "
    << bubbleDescription.back().auxiliarySlope << " }";
    std::cout << std::endl;/**/
    if( !definitelyUndershot )
    {
      definitelyUndershot = ( auxiliaryAndFirstDerivative[ 1 ] > 0.0 );
    }
    if( !definitelyOvershot )
    {
      definitelyOvershot = ( auxiliaryAndFirstDerivative[ 0 ] < 0.0 );
      if( definitelyOvershot )
      {
        overshootIndex = ( bubbleDescription.size() - 1 );
      }
    }
  }

  // This resets everything for a new integration, including recording the
  // values before the first Euclidean step performed before the
  // odeint::integrate.
  inline void
  OdeintBubbleObserver::ResetValues( double const initialAuxiliary )
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "OdeintBubbleObserver::ResetValues( " << initialAuxiliary
    << " ) called.";
    std::cout << std::endl;/**/
    bubbleDescription.resize( 1,
                              BubbleRadialValueDescription( 0.0,
                                                            initialAuxiliary,
                                                            0.0 ) );
    definitelyUndershot = false;
    definitelyOvershot = false;
    overshootIndex = 0;
    needsOrdering = false;
  }

  // This sorts the vectors by increasing radial value if needsOrdering is
  // true.
  inline void OdeintBubbleObserver::SortByRadialValue()
  {
    if( needsOrdering )
    {
      std::list< BubbleRadialValueDescription >
      sortableBubble( bubbleDescription.begin(),
                      bubbleDescription.end() );
      sortableBubble.sort(
                          &(BubbleRadialValueDescription::SortByRadialValue) );
      bubbleDescription.assign( sortableBubble.begin(),
                                sortableBubble.end() );
    }
  }

  // This is for debugging.
  inline std::string OdeintBubbleObserver::AsDebuggingString() const
  {
    std::stringstream returnStream;
    returnStream << "ordered properly: ";
    if( needsOrdering )
    {
      returnStream << " no";
    }
    else
    {
      returnStream << " yes";
    }
    returnStream << ", undershot: ";
    if( definitelyUndershot )
    {
      returnStream << " yes";
    }
    else
    {
      returnStream << " dunno";
    }
    returnStream << ", overshot: ";
    if( definitelyOvershot )
    {
      returnStream << " yes [" << overshootIndex << "]";
    }
    else
    {
      returnStream << " dunno";
    }
    returnStream << std::endl << "profile (" << bubbleDescription.size()
    << " elements): { ";
    for( std::vector< BubbleRadialValueDescription >::const_iterator
         bubbleBit( bubbleDescription.begin() );
         bubbleBit < bubbleDescription.end();
         ++bubbleBit )
    {
      if( bubbleBit != bubbleDescription.begin() )
      {
        returnStream << ", ";
      }
      returnStream << "[ r = " << bubbleBit->radialValue << ", p = "
      << bubbleBit->auxiliaryValue << ", dp/dr = "
      << bubbleBit->auxiliarySlope << " ]";
    }
    returnStream << " }";
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* ODEINTBUBBLEOBSERVER_HPP_ */
