/*
 * PolynomialPathThroughNodes.hpp
 *
 *  Created on: Jul 4, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALPATHTHROUGHNODES_HPP_
#define POLYNOMIALPATHTHROUGHNODES_HPP_

#include "CommonIncludes.hpp"
#include "PathFromNodes.hpp"

namespace VevaciousPlusPlus
{

  class PolynomialPathThroughNodes : public PathFromNodes
  {
  public:
    PolynomialPathThroughNodes();
    virtual
    ~PolynomialPathThroughNodes();
  };

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALPATHTHROUGHNODES_HPP_ */
