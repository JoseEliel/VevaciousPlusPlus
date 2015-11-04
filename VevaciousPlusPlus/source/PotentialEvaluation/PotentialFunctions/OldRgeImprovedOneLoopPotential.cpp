/*
 * OldRgeImprovedOneLoopPotential.cpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../../include/PotentialEvaluation/PotentialFunctions/OldRgeImprovedOneLoopPotential.hpp"

namespace VevaciousPlusPlus
{

  OldRgeImprovedOneLoopPotential::OldRgeImprovedOneLoopPotential(
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
    logarithmOfMinimumRenormalizationScale( 0.0 ),
    logarithmOfMaximumRenormalizationScale( 0.0 ),
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

  OldRgeImprovedOneLoopPotential::OldRgeImprovedOneLoopPotential(
   OldPotentialFromPolynomialAndMasses const& potentialFromPolynomialAndMasses ) :
    OldPotentialFromPolynomialAndMasses( potentialFromPolynomialAndMasses ),
    logarithmOfMinimumRenormalizationScale( log(
                                        currentMinimumRenormalizationScale ) ),
    logarithmOfMaximumRenormalizationScale( log(
                                        currentMaximumRenormalizationScale ) ),
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

  OldRgeImprovedOneLoopPotential::~OldRgeImprovedOneLoopPotential()
  {
    // This does nothing.
  }


  // This returns the energy density in GeV^4 of the potential for a state
  // strongly peaked around expectation values (in GeV) for the fields given
  // by the values of fieldConfiguration and temperature in GeV given by
  // temperatureValue.
  double OldRgeImprovedOneLoopPotential::operator()(
                               std::vector< double > const& fieldConfiguration,
                                          double const temperatureValue ) const
  {
    std::vector< double > cappedFieldConfiguration( fieldConfiguration );
    double const squaredLengthBeyondCap( CapFieldConfiguration(
                                                  cappedFieldConfiguration ) );
    double renormalizationScaleSquared( RenormalizationScaleSquared(
                                                      cappedFieldConfiguration,
                                                          temperatureValue ) );
    double logarithmOfScale( 0.5 * log( renormalizationScaleSquared ) );
    std::vector< DoubleVectorWithDouble > scalarMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( cappedFieldConfiguration,
                                      scalarSquareMasses,
                                      scalarMassesSquaredWithFactors,
                                      logarithmOfScale );
    std::vector< DoubleVectorWithDouble > fermionMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( cappedFieldConfiguration,
                                      fermionSquareMasses,
                                      fermionMassesSquaredWithFactors,
                                      logarithmOfScale );
    std::vector< DoubleVectorWithDouble > vectorMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( cappedFieldConfiguration,
                                      vectorSquareMasses,
                                      vectorMassesSquaredWithFactors,
                                      logarithmOfScale );
    return ( ( squaredLengthBeyondCap * squaredLengthBeyondCap )
             + treeLevelPotential( cappedFieldConfiguration,
                                 logarithmOfScale )
             + polynomialLoopCorrections( cappedFieldConfiguration,
                                          logarithmOfScale )
             + LoopAndThermalCorrections( cappedFieldConfiguration,
                                          scalarMassesSquaredWithFactors,
                                          fermionMassesSquaredWithFactors,
                                          vectorMassesSquaredWithFactors,
                                         ( 1.0 / renormalizationScaleSquared ),
                                          temperatureValue ) );
  }

  // This sets dsbFieldValueInputs based on the SLHA file just read in.
  void OldRgeImprovedOneLoopPotential::UpdateSelfForNewSlha(
                                               SlhaManager const& slhaManager )
  {
    currentMinimumRenormalizationScale = runningParameters.LowestBlockScale();
    squareOfMinimumRenormalizationScale = ( currentMinimumRenormalizationScale
                                        * currentMinimumRenormalizationScale );
    logarithmOfMinimumRenormalizationScale
    = log( currentMinimumRenormalizationScale );
    currentMaximumRenormalizationScale = runningParameters.HighestBlockScale();
    if( currentMaximumRenormalizationScale
        < ( scaleRangeMinimumFactor * currentMinimumRenormalizationScale ) )
    {
      currentMaximumRenormalizationScale
      = ( scaleRangeMinimumFactor * currentMinimumRenormalizationScale );
    }
    squareOfMaximumRenormalizationScale = ( currentMaximumRenormalizationScale
                                        * currentMaximumRenormalizationScale );
    logarithmOfMaximumRenormalizationScale
    = log( currentMaximumRenormalizationScale );
    std::vector< double > fieldOrigin( numberOfFields,
                                       0.0 );
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      dsbFieldValueInputs[ fieldIndex ]
      = dsbFieldValuePolynomials[ fieldIndex ]( fieldOrigin );
    }
  }

} /* namespace VevaciousPlusPlus */
