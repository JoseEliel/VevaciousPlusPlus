/*
 * MinimizingNormalForcesOnBisections.hpp
 *
 *  Created on: Oct 27, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINIMIZINGNORMALFORCESONBISECTIONS_HPP_
#define MINIMIZINGNORMALFORCESONBISECTIONS_HPP_

namespace VevaciousPlusPlus
{

  class MinimizingNormalForcesOnBisections
  {
  public:
    MinimizingNormalForcesOnBisections();
    virtual
    ~MinimizingNormalForcesOnBisections();

    // This has to do something taking into account the bubble profile.
    // Should move path false vacuum based on hyperplane perpendicular to line
    // from falseVacuum to pathFalseVacuum, ignoring bubble profile; also
    // should move bubble center configuration on hyperplane perpendicular to
    // line from trueVacuum to pathEnd, also ignoring bubble profile. No other
    // nodes beyond range of bubble profile, and full resolution of nodes
    // between pathFalseVacuum and pathEnd.
  };

} /* namespace VevaciousPlusPlus */
#endif /* MINIMIZINGNORMALFORCESONBISECTIONS_HPP_ */
