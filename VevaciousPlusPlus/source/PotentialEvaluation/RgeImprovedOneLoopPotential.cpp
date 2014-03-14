/*
 * RgeImprovedOneLoopPotential.cpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  RgeImprovedOneLoopPotential::RgeImprovedOneLoopPotential(
                                           std::string const& modelFilename ) :
    PotentialFromPolynomialAndMasses( modelFilename )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "RgeImprovedOneLoopPotential::RgeImprovedOneLoopPotential( \""
    << modelFilename << "\" )";
    std::cout << std::endl;/**/

  }

  RgeImprovedOneLoopPotential::~RgeImprovedOneLoopPotential()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "RgeImprovedOneLoopPotential::~RgeImprovedOneLoopPotential()";
    std::cout << std::endl;/**/
  }


  // This evaluates the target system and places the values in
  // destinationVector.
  void RgeImprovedOneLoopPotential::HomotopyContinuationSystemValues(
                                   std::vector< double > solutionConfiguration,
                                     std::vector< double >& destinationVector )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "RgeImprovedOneLoopPotential::HomotopyContinuationSystemValues( ... )";
    std::cout << std::endl;/**/
  }


  // This evaluates the derivatives of the target system and places the
  // values in destinationMatrix.
  void RgeImprovedOneLoopPotential::HomotopyContinuationSystemGradients(
                                   std::vector< double > solutionConfiguration,
                      std::vector< std::vector< double > >& destinationMatrix )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: RgeImprovedOneLoopPotential::"
    << "HomotopyContinuationSystemGradients( ... )";
    std::cout << std::endl;/**/
  }

  // This prepares a system of polynomials for the homotopy continuation
  // based on the current SLHA input data. Each polynomial term in the
  // tree-level potential generates its derivatives in its fields with the
  // coefficients fitted to a polynomial in the logarithm of the
  // renormalization scale, and then also a polynomial relating the logarithm
  // of the renormalization scale to minimumRenormalizationScaleSquared and
  // the field values is also prepared.
  void RgeImprovedOneLoopPotential::PrepareHomotopyContinuationPolynomials()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: RgeImprovedOneLoopPotential::"
    << "PrepareHomotopyContinuationPolynomials()";
    std::cout << std::endl;/**/
  }

} /* namespace VevaciousPlusPlus */
