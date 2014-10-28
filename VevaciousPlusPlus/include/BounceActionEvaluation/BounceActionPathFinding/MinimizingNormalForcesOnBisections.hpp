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

    // Something for all bisections: since by construction the nodes from the
    // TunnelPath are equally-spaced along the path in field space, the nodes
    // either side of any given node will be equidistant if the kinks in the
    // path between each pair average out as the same. Hence every intermediate
    // node is distant from its neighbors by an amount that should be
    // proportional to the rate of change of the path, so that's probably even
    // better than taking the proper bisection.

    // Back to stuff specific to MinimizingNormalForcesOnBisections:
    // 1) Remember to get sign of acceleration force correct! (Should go on
    // inner side of ridge, not outer side of valley!)
    // 2) Since we have dp/dr at start node and neighboring nodes, we can
    // construct df/dr from quadratic df/dp through end nodes plus varying
    // node, combined with dp/dr assuming that it is constant enough.
  };

} /* namespace VevaciousPlusPlus */
#endif /* MINIMIZINGNORMALFORCESONBISECTIONS_HPP_ */
