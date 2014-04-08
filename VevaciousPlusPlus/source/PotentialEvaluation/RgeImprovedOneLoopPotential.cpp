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
                                              std::string const& modelFilename,
                           RunningParameterManager& runningParameterManager ) :
    PotentialFromPolynomialAndMasses( modelFilename,
                                      runningParameterManager )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "RgeImprovedOneLoopPotential::RgeImprovedOneLoopPotential( \""
    << modelFilename << "\" )";
    std::cout << std::endl;/**/
  }

  RgeImprovedOneLoopPotential::RgeImprovedOneLoopPotential(
         PotentialFromPolynomialAndMasses& potentialFromPolynomialAndMasses ) :
    PotentialFromPolynomialAndMasses( potentialFromPolynomialAndMasses )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "RgeImprovedOneLoopPotential::RgeImprovedOneLoopPotential( [copy] )";
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


  // This should prepare homotopyContinuationPotentialPolynomial
  // appropriately.
  void
  RgeImprovedOneLoopPotential::PrepareHomotopyContinuationPotentialPolynomial()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: RgeImprovedOneLoopPotential::"
    << "PrepareHomotopyContinuationPotentialPolynomial()";
    std::cout << std::endl;/**/
  }

  // This should prepare homotopyContinuationStartSystem and
  // startPolynomialHessian appropriately.
  void RgeImprovedOneLoopPotential::PrepareHomotopyContinuationStartSystem()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: RgeImprovedOneLoopPotential::"
    << "PrepareHomotopyContinuationStartSystem()";
    std::cout << std::endl;/**/
  }

  // This should prepare homotopyContinuationStartValues to be all the
  // solutions of homotopyContinuationStartSystem.
  void RgeImprovedOneLoopPotential::PrepareHomotopyContinuationStartValues()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: RgeImprovedOneLoopPotential::"
    << "PrepareHomotopyContinuationStartValues()";
    std::cout << std::endl;/**/
  }

} /* namespace VevaciousPlusPlus */
