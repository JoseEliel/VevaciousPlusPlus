/*
 * SlhaDerivedFunctionoid.hpp
 *
 *  Created on: Oct 23, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHADERIVEDFUNCTIONOID_HPP_
#define SLHADERIVEDFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"
#include "SlhaSourcedParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaDerivedFunctionoid : public SlhaSourcedParameterFunctionoid
  {
  public:
    SlhaDerivedFunctionoid();
    virtual ~SlhaDerivedFunctionoid();


    // This is false until the functionoid has been assigned an index as an
    // active Lagrangian parameter.
    bool IsActive() const { return isActive; }

    // This sets indexInValuesVector to activeIndex and sets isActive to true.
    void MakeActive( size_t const activeIndex );


  protected:
    bool isActive;
  };





  // This sets indexInValuesVector to activeIndex and sets isActive to true.
  inline void SlhaDerivedFunctionoid::MakeActive( size_t const activeIndex )
  {
    indexInValuesVector = activeIndex;
    isActive = true;
  }

} /* namespace VevaciousPlusPlus */

#endif /* SLHADERIVEDFUNCTIONOID_HPP_ */
