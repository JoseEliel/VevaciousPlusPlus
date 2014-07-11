/*
 * MinuitPathPotentialMinimizer.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITPATHPOTENTIALMINIMIZER_HPP_
#define MINUITPATHPOTENTIALMINIMIZER_HPP_

#include "CommonIncludes.hpp"
#include "BouncePathFinder.hpp"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "PotentialMinimization/GradientBasedMinimization/MinuitMinimum.hpp"
#include "PathFromNodesFactory.hpp"
#include "NodesFromParameterization.hpp"
#include "LinearSplineThroughNodesFactory.hpp"
#include "QuadraticSplineThroughNodesFactory.hpp"
#include "PolynomialThroughNodesFactory.hpp"

namespace VevaciousPlusPlus
{

  class MinuitPathPotentialMinimizer
  {
  public:
    MinuitPathPotentialMinimizer();
    virtual
    ~MinuitPathPotentialMinimizer();
  };

} /* namespace VevaciousPlusPlus */
#endif /* MINUITPATHPOTENTIALMINIMIZER_HPP_ */
