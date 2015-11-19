/*
 * OldFixedScaleOneLoopPotential.cpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/PotentialFunctions/OldFixedScaleOneLoopPotential.hpp"

namespace VevaciousPlusPlus
{

  OldFixedScaleOneLoopPotential::OldFixedScaleOneLoopPotential(
                                              std::string const& modelFilename,
                                          double const scaleRangeMinimumFactor,
            bool const treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions,
                               double const assumedPositiveOrNegativeTolerance,
                           RunningParameterManager& runningParameterManager ) :
    OldPotentialFromPolynomialAndMasses( modelFilename,
                                      scaleRangeMinimumFactor,
                       treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions,
                                      assumedPositiveOrNegativeTolerance,
                                      runningParameterManager ),
    inverseRenormalizationScaleSquared( -1.0 ),
    homotopyContinuationTargetSystem( treeLevelPotential,
                                      numberOfFields,
                                      *this,
                                      fieldsAssumedPositive,
                                      fieldsAssumedNegative,
                       treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions,
                                      assumedPositiveOrNegativeTolerance )
  {
    // This constructor is just an initialization list.
  }

  OldFixedScaleOneLoopPotential::OldFixedScaleOneLoopPotential(
      OldPotentialFromPolynomialAndMasses& potentialFromPolynomialAndMasses ) :
    OldPotentialFromPolynomialAndMasses( potentialFromPolynomialAndMasses ),
    inverseRenormalizationScaleSquared( -1.0 ),
    homotopyContinuationTargetSystem( treeLevelPotential,
                                      numberOfFields,
                                      *this,
                                      fieldsAssumedPositive,
                                      fieldsAssumedNegative,
                      treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions,
                                      assumedPositiveOrNegativeTolerance )
  {
    // This constructor is just an initialization list.
  }

  OldFixedScaleOneLoopPotential::~OldFixedScaleOneLoopPotential()
  {
    // This does nothing.
  }


  // This returns the energy density in GeV^4 of the potential for a state
  // strongly peaked around expectation values (in GeV) for the fields given
  // by the values of fieldConfiguration and temperature in GeV given by
  // temperatureValue.
  double OldFixedScaleOneLoopPotential::operator()(
                               std::vector< double > const& fieldConfiguration,
                                          double const temperatureValue ) const
  {
    std::vector< double > cappedFieldConfiguration( fieldConfiguration );
    double const squaredLengthBeyondCap( CapFieldConfiguration(
                                                  cappedFieldConfiguration ) );

    std::vector< DoubleVectorWithDouble > scalarMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( cappedFieldConfiguration,
                                      scalarSquareMasses,
                                      scalarMassesSquaredWithFactors );
    std::vector< DoubleVectorWithDouble > fermionMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( cappedFieldConfiguration,
                                      fermionSquareMasses,
                                      fermionMassesSquaredWithFactors );
    std::vector< DoubleVectorWithDouble > vectorMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( cappedFieldConfiguration,
                                      vectorSquareMasses,
                                      vectorMassesSquaredWithFactors );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl << "OldFixedScaleOneLoopPotential::operator()( "
    << FieldConfigurationAsMathematica( fieldConfiguration )
    << ", T =" << temperatureValue
    << " ) called. squaredLengthBeyondCap = " << squaredLengthBeyondCap
    << ", cappedFieldConfiguration = "
    << FieldConfigurationAsMathematica( cappedFieldConfiguration )
    << std::endl << "scalarMassesSquaredWithFactors = {" << std::endl;
    for( std::vector< DoubleVectorWithDouble >::const_iterator
         massesWithFactor( scalarMassesSquaredWithFactors.begin() );
         massesWithFactor < scalarMassesSquaredWithFactors.end();
         ++massesWithFactor )
    {
      std::cout
      << "[ factor: " << massesWithFactor->second
      << ", masses-squared: { ";
      for( std::vector< double >::const_iterator
           massSquared( massesWithFactor->first.begin() );
           massSquared < massesWithFactor->first.end();
           ++massSquared )
      {
        std::cout << *massSquared << ",";
      }
      std::cout << " } ]" << std::endl;
    }
    std::cout << "}"
    << std::endl << "fermionMassesSquaredWithFactors = {" << std::endl;
    for( std::vector< DoubleVectorWithDouble >::const_iterator
         massesWithFactor( fermionMassesSquaredWithFactors.begin() );
         massesWithFactor < fermionMassesSquaredWithFactors.end();
         ++massesWithFactor )
    {
      std::cout
      << "[ factor: " << massesWithFactor->second
      << ", masses-squared: { ";
      for( std::vector< double >::const_iterator
           massSquared( massesWithFactor->first.begin() );
           massSquared < massesWithFactor->first.end();
           ++massSquared )
      {
        std::cout << *massSquared << ",";
      }
      std::cout << " } ]" << std::endl;
    }
    std::cout << "}"
    << std::endl << "vectorMassesSquaredWithFactors = {" << std::endl;
    for( std::vector< DoubleVectorWithDouble >::const_iterator
         massesWithFactor( vectorMassesSquaredWithFactors.begin() );
         massesWithFactor < vectorMassesSquaredWithFactors.end();
         ++massesWithFactor )
    {
      std::cout
      << "[ factor: " << massesWithFactor->second
      << ", masses-squared: { ";
      for( std::vector< double >::const_iterator
           massSquared( massesWithFactor->first.begin() );
           massSquared < massesWithFactor->first.end();
           ++massSquared )
      {
        std::cout << *massSquared << ",";
      }
      std::cout << " } ]" << std::endl;
    }
    std::cout << " }";
    std::cout << std::endl;
    std::cout << "treeLevelPotential = "
    << treeLevelPotential( cappedFieldConfiguration );
    std::cout << std::endl;
    std::cout << "polynomialLoopCorrections = "
    << polynomialLoopCorrections( cappedFieldConfiguration )
    // << std::endl << polynomialLoopCorrections.AsDebuggingString();
    << std::endl;
    std::cout << "LoopAndThermalCorrections = "
    << LoopAndThermalCorrections( cappedFieldConfiguration,
                                  scalarMassesSquaredWithFactors,
                                  fermionMassesSquaredWithFactors,
                                  vectorMassesSquaredWithFactors,
                                  inverseRenormalizationScaleSquared,
                                  temperatureValue );
    std::cout << std::endl;*/

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
