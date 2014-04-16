/*
 * BasicPolynomialHomotopyContinuation.cpp
 *
 *  Created on: Apr 14, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  BasicPolynomialHomotopyContinuation::BasicPolynomialHomotopyContinuation(
                    HomotopyContinuationTargetSystem& polynomialPotential ) :
    HomotopyContinuationSolver( polynomialPotential )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "BasicPolynomialHomotopyContinuation::"
    << "BasicPolynomialHomotopyContinuation( ... )";
    std::cout << std::endl;/**/
  }

  BasicPolynomialHomotopyContinuation::~BasicPolynomialHomotopyContinuation()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "BasicPolynomialHomotopyContinuation::"
    << "~BasicPolynomialHomotopyContinuation()";
    std::cout << std::endl;/**/
  }


  // This should find all the extrema of homotopyContinuationPotential and
  // put the purely real solutions into purelyRealSolutionSets.
  void BasicPolynomialHomotopyContinuation::FindTreeLevelExtrema(
                 std::vector< std::vector< double > >& purelyRealSolutionSets )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "BasicPolynomialHomotopyContinuation::FindTreeLevelExtrema( ... )";
    std::cout << std::endl;/**/
  }
} /* namespace VevaciousPlusPlus */
