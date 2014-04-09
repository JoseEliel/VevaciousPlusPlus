/*
 * FixedScaleOneLoopPotential.cpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  FixedScaleOneLoopPotential::FixedScaleOneLoopPotential(
                                              std::string const& modelFilename,
                           RunningParameterManager& runningParameterManager ) :
    PotentialFromPolynomialAndMasses( modelFilename,
                                      runningParameterManager )
  {
    // This constructor is just an initialization list.
  }

  FixedScaleOneLoopPotential::FixedScaleOneLoopPotential(
         PotentialFromPolynomialAndMasses& potentialFromPolynomialAndMasses ) :
    PotentialFromPolynomialAndMasses( potentialFromPolynomialAndMasses )
  {
    // This constructor is just an initialization list.
  }

  FixedScaleOneLoopPotential::~FixedScaleOneLoopPotential()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "FixedScaleOneLoopPotential::~FixedScaleOneLoopPotential()";
    std::cout << std::endl;/**/
  }


  // This should prepare homotopyContinuationStartSystem and
  // startPolynomialHessian appropriately.
  void FixedScaleOneLoopPotential::PrepareHomotopyContinuationStartSystem()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: FixedScaleOneLoopPotential::"
    << "PrepareHomotopyContinuationStartSystem()";
    std::cout << std::endl;/**/
  }

  // This should prepare homotopyContinuationStartValues to be all the
  // solutions of homotopyContinuationStartSystem.
  void FixedScaleOneLoopPotential::PrepareHomotopyContinuationStartValues()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: FixedScaleOneLoopPotential::"
    << "PrepareHomotopyContinuationStartValues()";
    std::cout << std::endl;/**/
  }

} /* namespace VevaciousPlusPlus */
