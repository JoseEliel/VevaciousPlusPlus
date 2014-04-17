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
    logarithmOfMaximumRenormalizationScale( NAN ),
    homotopyContinuationTargetSystem( treeLevelPotential,
                                      numberOfFields,
                                      *this )
  {
    // This constructor is just an initialization list.
  }

  RgeImprovedOneLoopPotential::RgeImprovedOneLoopPotential(
   PotentialFromPolynomialAndMasses const& potentialFromPolynomialAndMasses ) :
    PotentialFromPolynomialAndMasses( potentialFromPolynomialAndMasses ),
    logarithmOfMinimumRenormalizationScale( log(
                                        currentMinimumRenormalizationScale ) ),
    logarithmOfMaximumRenormalizationScale( log(
                                        currentMaximumRenormalizationScale ) ),
    homotopyContinuationTargetSystem( treeLevelPotential,
                                      numberOfFields,
                                      *this )
  {
    // This constructor is just an initialization list.
  }

  RgeImprovedOneLoopPotential::~RgeImprovedOneLoopPotential()
  {
    // This does nothing.
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

} /* namespace VevaciousPlusPlus */
