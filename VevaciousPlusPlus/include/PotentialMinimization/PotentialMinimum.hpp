/*
 * PotentialMinimum.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POTENTIALMINIMUM_HPP_
#define POTENTIALMINIMUM_HPP_

#include "../StandardIncludes.hpp"

namespace VevaciousPlusPlus
{
  class PotentialMinimum
  {
  public:
    PotentialMinimum();
    virtual
    ~PotentialMinimum();


    // This returns the sum of the squares of the differences in the field
    // values of this PotentialMinimum with comparisonMinimum.
    double SquareDistanceTo( PotentialMinimum const& comparisonMinimum ) const;
  };

} /* namespace VevaciousPlusPlus */
#endif /* POTENTIALMINIMUM_HPP_ */
