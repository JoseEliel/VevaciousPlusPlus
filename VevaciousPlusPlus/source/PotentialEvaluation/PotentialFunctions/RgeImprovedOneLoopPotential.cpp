/*
 * RgeImprovedOneLoopPotential.cpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  RgeImprovedOneLoopPotential::RgeImprovedOneLoopPotential(
                                              std::string const& modelFilename,
                           RunningParameterManager& runningParameterManager ) :
    PotentialFromPolynomialAndMasses( modelFilename,
                                      runningParameterManager ),
    logarithmOfMinimumRenormalizationScale( NAN ),
    logarithmOfMaximumRenormalizationScale( NAN )
  {
    // This constructor is just an initialization list.
  }

  RgeImprovedOneLoopPotential::RgeImprovedOneLoopPotential(
   PotentialFromPolynomialAndMasses const& potentialFromPolynomialAndMasses ) :
    PotentialFromPolynomialAndMasses( potentialFromPolynomialAndMasses ),
    logarithmOfMinimumRenormalizationScale( log(
                                        currentMinimumRenormalizationScale ) ),
    logarithmOfMaximumRenormalizationScale( log(
                                         currentMaximumRenormalizationScale ) )
  {
    // This constructor is just an initialization list.
  }

  RgeImprovedOneLoopPotential::~RgeImprovedOneLoopPotential()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "RgeImprovedOneLoopPotential::~RgeImprovedOneLoopPotential()";
    std::cout << std::endl;/**/
  }


  // This returns the energy density in GeV^4 of the potential for a state
  // strongly peaked around expectation values (in GeV) for the fields given
  // by the values of fieldConfiguration and temperature in GeV given by
  // temperatureValue.
  double RgeImprovedOneLoopPotential::operator()(
                               std::vector< double > const& fieldConfiguration,
                                          double const temperatureValue ) const
  {
    double renormalizationScaleSquared( RenormalizationScaleSquared(
                                                            fieldConfiguration,
                                                          temperatureValue ) );
    double logarithmOfScale( 0.5 * log( renormalizationScaleSquared ) );
    std::vector< DoubleVectorWithDouble > scalarMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( fieldConfiguration,
                                      scalarSquareMasses,
                                      scalarMassesSquaredWithFactors,
                                      logarithmOfScale );
    std::vector< DoubleVectorWithDouble > fermionMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( fieldConfiguration,
                                      fermionSquareMasses,
                                      fermionMassesSquaredWithFactors,
                                      logarithmOfScale );
    std::vector< DoubleVectorWithDouble > vectorMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( fieldConfiguration,
                                      vectorSquareMasses,
                                      vectorMassesSquaredWithFactors,
                                      logarithmOfScale );
    return ( treeLevelPotential( fieldConfiguration,
                                 logarithmOfScale )
             + polynomialLoopCorrections( fieldConfiguration,
                                          logarithmOfScale )
             + LoopAndThermalCorrections( fieldConfiguration,
                                          scalarMassesSquaredWithFactors,
                                          fermionMassesSquaredWithFactors,
                                          vectorMassesSquaredWithFactors,
                                         ( 1.0 / renormalizationScaleSquared ),
                                          temperatureValue ) );
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
