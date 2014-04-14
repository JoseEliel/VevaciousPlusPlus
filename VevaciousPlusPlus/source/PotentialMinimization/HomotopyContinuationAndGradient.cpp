/*
 * HomotopyContinuationAndGradient.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  HomotopyContinuationAndGradient::HomotopyContinuationAndGradient(
                                          PotentialFunction& potentialFunction,
                     HomotopyContinuationSolver& homotopyContinuationSolver ) :
    PotentialMinimizer( potentialFunction ),
    homotopyContinuationSolver( homotopyContinuationSolver ),
    purelyRealSolutionSets()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "HomotopyContinuationAndGradient::HomotopyContinuationAndGradient("
    << " ... )";
    std::cout << std::endl;/**/
  }

  HomotopyContinuationAndGradient::~HomotopyContinuationAndGradient()
  {
    // This does nothing.
  }

} /* namespace VevaciousPlusPlus */
