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
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "FixedScaleOneLoopPotential::FixedScaleOneLoopPotential( \""
    << modelFilename << "\" )";
    std::cout << std::endl;/**/
  }

  FixedScaleOneLoopPotential::FixedScaleOneLoopPotential(
         PotentialFromPolynomialAndMasses& potentialFromPolynomialAndMasses ) :
    PotentialFromPolynomialAndMasses( potentialFromPolynomialAndMasses )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "FixedScaleOneLoopPotential::FixedScaleOneLoopPotential( [copy] )";
    std::cout << std::endl;/**/
  }

  FixedScaleOneLoopPotential::~FixedScaleOneLoopPotential()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "FixedScaleOneLoopPotential::~FixedScaleOneLoopPotential()";
    std::cout << std::endl;/**/
  }


  // This evaluates the target system and places the values in
  // destinationVector.
  void FixedScaleOneLoopPotential::HomotopyContinuationSystemValues(
                                   std::vector< double > solutionConfiguration,
                                     std::vector< double >& destinationVector )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "FixedScaleOneLoopPotential::HomotopyContinuationSystemValues( ... )";
    std::cout << std::endl;/**/
  }

  // This evaluates the derivatives of the target system and places the
  // values in destinationMatrix.
  void FixedScaleOneLoopPotential::HomotopyContinuationSystemGradients(
                                   std::vector< double > solutionConfiguration,
                      std::vector< std::vector< double > >& destinationMatrix )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: FixedScaleOneLoopPotential::"
    << "HomotopyContinuationSystemGradients( ... )";
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
