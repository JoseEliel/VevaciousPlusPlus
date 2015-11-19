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
                      OldPotentialFromPolynomialAndMasses const& copySource ) :
    OldPotentialFromPolynomialAndMasses( copySource ),
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

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl << "OldRgeImprovedOneLoopPotential::operator()( "
    << FieldConfigurationAsMathematica( fieldConfiguration )
    << ", T =" << temperatureValue
    << " ) called. squaredLengthBeyondCap = " << squaredLengthBeyondCap
    << ", cappedFieldConfiguration = "
    << FieldConfigurationAsMathematica( cappedFieldConfiguration )
    << ", renormalizationScaleSquared = "
    << renormalizationScaleSquared
    << ", sqrt( renormalizationScaleSquared ) = "
    << sqrt( renormalizationScaleSquared )
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
                                  ( 1.0 / renormalizationScaleSquared ),
                                  temperatureValue );
    std::cout << std::endl;*/

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

    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "At end of OldRgeImprovedOneLoopPotential::UpdateSelfForNewSlha(...),"
    << " currentMinimumRenormalizationScale = "
    << currentMinimumRenormalizationScale
    << ", squareOfMinimumRenormalizationScale = "
    << squareOfMinimumRenormalizationScale
    << ", logarithmOfMinimumRenormalizationScale = "
    << logarithmOfMinimumRenormalizationScale
    << ", currentMaximumRenormalizationScale = "
    << currentMaximumRenormalizationScale
    << ", squareOfMaximumRenormalizationScale = "
    << squareOfMaximumRenormalizationScale
    << ", logarithmOfMaximumRenormalizationScale = "
    << logarithmOfMaximumRenormalizationScale;
    std::cout << std::endl;/**/
  }

} /* namespace VevaciousPlusPlus */
