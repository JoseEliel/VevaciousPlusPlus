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


  // This should return a vector of field values corresponding to the field
  // configuration as it should be passed to operator() for evaluating the
  // potential, given a vector of values that solve this instance's homotopy
  // continuation system. It should return an empty vector if
  // homotopyContinuatioConfiguration does not correspond to a valid field
  // configuration. (For example, RgeImprovedOneLoopPotential uses the
  // logarithm of the renormalization scale as an extra variable in the
  // homotopy continuation, so homotopyContinuatioConfiguration actually has
  // an extra entry compared to a valid field configuration. Also, it may
  // be that homotopyContinuatioConfiguration did not correspond to a valid
  // solution where the scale is close to the Euclidean length of the field
  // configuration, so it would not correspond to a valid solution.)
  std::vector< double >
  RgeImprovedOneLoopPotential::ValidFieldsFromHomotopyContinuation(
              std::vector< double > homotopyContinuatioConfiguration ) const
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "RgeImprovedOneLoopPotential::ValidFieldsFromHomotopyContinuation( ..."
    << " )";
    std::cout << std::endl;
    return std::vector< double >();/**/
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
