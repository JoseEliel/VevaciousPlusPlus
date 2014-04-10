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
                                      runningParameterManager ),
    inverseRenormalizationScaleSquared( NAN )
  {
    // This constructor is just an initialization list.
  }

  FixedScaleOneLoopPotential::FixedScaleOneLoopPotential(
         PotentialFromPolynomialAndMasses& potentialFromPolynomialAndMasses ) :
    PotentialFromPolynomialAndMasses( potentialFromPolynomialAndMasses ),
    inverseRenormalizationScaleSquared( NAN )
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


  // This returns the energy density in GeV^4 of the potential for a state
  // strongly peaked around expectation values (in GeV) for the fields given
  // by the values of fieldConfiguration and temperature in GeV given by
  // temperatureValue.
  double FixedScaleOneLoopPotential::operator()(
                               std::vector< double > const& fieldConfiguration,
                                          double const temperatureValue ) const
  {
    std::vector< DoubleVectorWithDouble > scalarMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( fieldConfiguration,
                                      scalarSquareMasses,
                                      scalarMassesSquaredWithFactors );
    std::vector< DoubleVectorWithDouble > fermionMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( fieldConfiguration,
                                      fermionMasses,
                                      fermionMassesSquaredWithFactors );
    AddMassesSquaredWithMultiplicity( fieldConfiguration,
                                      fermionMassSquareds,
                                      fermionMassesSquaredWithFactors );
    std::vector< DoubleVectorWithDouble > vectorMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( fieldConfiguration,
                                      vectorSquareMasses,
                                      vectorMassesSquaredWithFactors );
    return ( treeLevelPotential( fieldConfiguration )
             + polynomialLoopCorrections( fieldConfiguration )
             + LoopAndThermalCorrections( fieldConfiguration,
                                          scalarMassesSquaredWithFactors,
                                          fermionMassesSquaredWithFactors,
                                          vectorMassesSquaredWithFactors,
                                          inverseRenormalizationScaleSquared,
                                          temperatureValue ) );
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
