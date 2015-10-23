/*
 * SlhaDerivedFunctionoid.hpp
 *
 *  Created on: Oct 23, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHADERIVEDFUNCTIONOID_HPP_
#define SLHADERIVEDFUNCTIONOID_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class SlhaDerivedFunctionoid
  {
  public:
    SlhaDerivedFunctionoid();
    virtual ~SlhaDerivedFunctionoid();


    size_t IndexInValuesVector() const { return indexInValuesVector; }

    // This returns the value of the functionoid for the given logarithm of the
    // scale.
    virtual double operator()( double const logarithmOfScale ) const = 0;

    // This is for creating a Python version of the potential.
    virtual std::string PythonParameterEvaluation() const = 0;

    // This is false until the functionoid has been assigned an index as an
    // active Lagrangian parameter.
    bool IsActive() const { return isActive; }

    // This sets indexInValuesVector to activeIndex and sets isActive to true.
    void MakeActive( size_t const activeIndex );


  protected:
    size_t indexInValuesVector;
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
