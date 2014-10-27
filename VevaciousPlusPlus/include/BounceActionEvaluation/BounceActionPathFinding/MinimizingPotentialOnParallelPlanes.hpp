/*
 * MinimizingPotentialOnParallelPlanes.hpp
 *
 *  Created on: Oct 27, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINIMIZINGPOTENTIALONPARALLELPLANES_HPP_
#define MINIMIZINGPOTENTIALONPARALLELPLANES_HPP_

namespace VevaciousPlusPlus
{

  class MinimizingPotentialOnParallelPlanes
  {
  public:
    MinimizingPotentialOnParallelPlanes();
    virtual
    ~MinimizingPotentialOnParallelPlanes();

    // This should minimize on a series of hyperplanes all perpendicular to
    // the vector difference of the vacua, but maybe to avoid "zig-zag-iness"
    // from Minuit2 coming to rest on a jittery line reasonably far from the
    // starting path, each plane should have its starting point be (previous
    // node + vector difference of previous node from the node before it), with
    // special cases for the first and last varying nodes, of course.
  };

} /* namespace VevaciousPlusPlus */
#endif /* MINIMIZINGPOTENTIALONPARALLELPLANES_HPP_ */
