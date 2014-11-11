/*
 * MinuitOnPathNormalInertialPotential.hpp
 *
 *  Created on: Nov 10, 2014
 *      Author: bol
 */

#ifndef MINUITONPATHNORMALINERTIALPOTENTIAL_HPP_
#define MINUITONPATHNORMALINERTIALPOTENTIAL_HPP_

#include "CommonIncludes.hpp"
#include "MinuitOnPotentialPerpendicularToPath.hpp"

namespace VevaciousPlusPlus
{

  class MinuitOnPathNormalInertialPotential
      : public MinuitOnPotentialPerpendicularToPath
  {
  public:
    MinuitOnPathNormalInertialPotential();
    virtual ~MinuitOnPathNormalInertialPotential();

    // This class should add a
    // (function of dp/dr and d^2p/dr^2) * (the radius of the Minuit vector)
    // to the potential (or subtract - make sure that the sign is correct!)
    // such that the path-velocity-and-acceleration-dependent equations of
    // motion perpendicular to the path are accounted for, assuming that they
    // are constant near lastPath and given by bubbleFromLastPath. The only
    // other thing it should do beyond MinuitOnPotentialPerpendicularToPath is
    // to note bubbleFromLastPath when TryToImprovePath is called, for the
    // above. (Since we have dp/dr at start node and neighboring nodes, we can
    // construct df/dr from quadratic df/dp through end nodes plus varying
    // node, combined with dp/dr assuming that it is constant enough.)
  };

} /* namespace VevaciousPlusPlus */

#endif /* MINUITONPATHNORMALINERTIALPOTENTIAL_HPP_ */
