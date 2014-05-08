/*
 * FixedScaleOneLoopPotential.cpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  FixedScaleOneLoopPotential::FixedScaleOneLoopPotential(
                                              std::string const& modelFilename,
                           RunningParameterManager& runningParameterManager ) :
    PotentialFromPolynomialAndMasses( modelFilename,
                                      runningParameterManager ),
    inverseRenormalizationScaleSquared( NAN ),
    homotopyContinuationTargetSystem( treeLevelPotential,
                                      numberOfFields,
                                      *this )
  {
    // This constructor is just an initialization list.
  }

  FixedScaleOneLoopPotential::FixedScaleOneLoopPotential(
         PotentialFromPolynomialAndMasses& potentialFromPolynomialAndMasses ) :
    PotentialFromPolynomialAndMasses( potentialFromPolynomialAndMasses ),
    inverseRenormalizationScaleSquared( NAN ),
    homotopyContinuationTargetSystem( treeLevelPotential,
                                      numberOfFields,
                                      *this )
  {
    // This constructor is just an initialization list.
  }

  FixedScaleOneLoopPotential::~FixedScaleOneLoopPotential()
  {
    // This does nothing.
  }


  // This returns the energy density in GeV^4 of the potential for a state
  // strongly peaked around expectation values (in GeV) for the fields given
  // by the values of fieldConfiguration and temperature in GeV given by
  // temperatureValue.
  double FixedScaleOneLoopPotential::operator()(
                               std::vector< double > const& fieldConfiguration,
                                          double const temperatureValue ) const
  {
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "FixedScaleOneLoopPotential::operator()( {";
    for( std::vector< double >::const_iterator
         fieldValue( fieldConfiguration.begin() );
         fieldValue < fieldConfiguration.end();
         ++fieldValue )
    {
      std::cout << "  " << *fieldValue;
    }
    std::cout << "  }, " << temperatureValue
    << " ) called. treeLevelPotential( fieldConfiguration ) = "
    << treeLevelPotential( fieldConfiguration );
    std::cout << std::endl;*/

    std::vector< double > cappedFieldConfiguration( fieldConfiguration );
    double const squaredLengthBeyondCap( CapFieldConfiguration(
                                                  cappedFieldConfiguration ) );

    std::vector< DoubleVectorWithDouble > scalarMassesSquaredWithFactors;

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << " scalarSquareMasses";
    std::cout << std::endl;*/
    AddMassesSquaredWithMultiplicity( cappedFieldConfiguration,
                                      scalarSquareMasses,
                                      scalarMassesSquaredWithFactors );
    std::vector< DoubleVectorWithDouble > fermionMassesSquaredWithFactors;

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "FixedScaleOneLoopPotential::operator():"
    << " fermionMasses";
    std::cout << std::endl;*/
    AddMassesSquaredWithMultiplicity( cappedFieldConfiguration,
                                      fermionSquareMasses,
                                      fermionMassesSquaredWithFactors );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "FixedScaleOneLoopPotential::operator():"
    << " vectorSquareMasses";
    std::cout << std::endl;*/
    std::vector< DoubleVectorWithDouble > vectorMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( cappedFieldConfiguration,
                                      vectorSquareMasses,
                                      vectorMassesSquaredWithFactors );
    return ( ( squaredLengthBeyondCap * squaredLengthBeyondCap )
             + treeLevelPotential( cappedFieldConfiguration )
             + polynomialLoopCorrections( cappedFieldConfiguration )
             + LoopAndThermalCorrections( cappedFieldConfiguration,
                                          scalarMassesSquaredWithFactors,
                                          fermionMassesSquaredWithFactors,
                                          vectorMassesSquaredWithFactors,
                                          inverseRenormalizationScaleSquared,
                                          temperatureValue ) );
  }

} /* namespace VevaciousPlusPlus */
