/*
 * StartingPointFinder.hpp
 *
 *  Created on: Jun 30, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef STARTINGPOINTFINDER_HPP_
#define STARTINGPOINTFINDER_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class StartingPointFinder
  {
  public:
    StartingPointFinder() { ; /* This does nothing. */ }
    virtual ~StartingPointFinder() { ; /* This does nothing. */ }


    // This should find all the starting points for the gradient-based
    // minimizer and put them into startingPoints.
    virtual void operator()(
              std::vector< std::vector< double > >& startingPoints ) const = 0;
  };

} /* namespace VevaciousPlusPlus */
#endif /* STARTINGPOINTFINDER_HPP_ */
