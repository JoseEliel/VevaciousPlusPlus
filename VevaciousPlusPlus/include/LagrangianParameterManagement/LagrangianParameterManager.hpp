/*
 * LagrangianParameterManager.hpp
 *
 *  Created on: Oct 9, 2015
 *      Author: bol
 */

#ifndef LAGRANGIANPARAMETERMANAGER_HPP_
#define LAGRANGIANPARAMETERMANAGER_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  class LagrangianParameterManager : public BOL::BasicObserved
  {
  public:
    LagrangianParameterManager();
    virtual ~LagrangianParameterManager();
  };

} /* namespace VevaciousPlusPlus */

#endif /* LAGRANGIANPARAMETERMANAGER_HPP_ */
