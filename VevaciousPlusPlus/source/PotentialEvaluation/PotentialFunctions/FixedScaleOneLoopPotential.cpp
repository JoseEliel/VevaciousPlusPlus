/*
 * FixedScaleOneLoopPotential.cpp
 *
 *  Created on: Nov 4, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/PotentialFunctions/FixedScaleOneLoopPotential.hpp"

namespace VevaciousPlusPlus
{

  FixedScaleOneLoopPotential::FixedScaleOneLoopPotential(
                                              std::string const& modelFilename,
                               double const assumedPositiveOrNegativeTolerance,
                     LagrangianParameterManager& lagrangianParameterManager ) :
    PotentialFromPolynomialWithMasses( modelFilename,
                                       assumedPositiveOrNegativeTolerance,
                                       lagrangianParameterManager ),
    BOL::BasicObserver(),
    renormalizationScale( -1.0 ),
    inverseRenormalizationScaleSquared( -1.0 )
  {
    lagrangianParameterManager.registerObserver( this );
  }

  FixedScaleOneLoopPotential::FixedScaleOneLoopPotential(
                   PotentialFromPolynomialWithMasses const& potentialToCopy ) :
    PotentialFromPolynomialWithMasses( potentialToCopy ),
    BOL::BasicObserver(),
    renormalizationScale( -1.0 ),
    inverseRenormalizationScaleSquared( -1.0 )
  {
    lagrangianParameterManager.registerObserver( this );
  }

  FixedScaleOneLoopPotential::~FixedScaleOneLoopPotential()
  {
    // This does nothing.
  }


  // This updates the scale used for the loop corrections based on the
  // appropriate scale from lagrangianParameterManager, and updates all
  // components used to evaluate the potential to use the Lagrangian
  // parameters evaluated at that scale.
  void FixedScaleOneLoopPotential::respondToObservedSignal()
  {
    renormalizationScale
    = lagrangianParameterManager.AppropriateSingleFixedScale();
    inverseRenormalizationScaleSquared
    = ( 1.0 / ( renormalizationScale * renormalizationScale ) );
    std::vector< double > fixedParameterValues;
    lagrangianParameterManager.ParameterValues( log( renormalizationScale ),
                                                fixedParameterValues );
    UpdateDsbValues( log( renormalizationScale ) );

    treeLevelPotential.UpdateForFixedScale( fixedParameterValues );
    polynomialLoopCorrections.UpdateForFixedScale( fixedParameterValues );
    for( std::vector< RealMassesSquaredMatrix >::iterator
         massMatrix( scalarMassSquaredMatrices.begin() );
         massMatrix < scalarMassSquaredMatrices.end();
         ++massMatrix )
    {
      massMatrix->UpdateForFixedScale( fixedParameterValues );
    }
    for( std::vector< SymmetricComplexMassMatrix >::iterator
         massMatrix( fermionMassMatrices.begin() );
         massMatrix < fermionMassMatrices.end();
         ++massMatrix )
    {
      massMatrix->UpdateForFixedScale( fixedParameterValues );
    }
    for( std::vector< ComplexMassSquaredMatrix >::iterator
         massMatrix( fermionMassSquaredMatrices.begin() );
         massMatrix < fermionMassSquaredMatrices.end();
         ++massMatrix )
    {
      massMatrix->UpdateForFixedScale( fixedParameterValues );
    }
    for( std::vector< RealMassesSquaredMatrix >::iterator
         massMatrix( vectorMassSquaredMatrices.begin() );
         massMatrix < vectorMassSquaredMatrices.end();
         ++massMatrix )
    {
      massMatrix->UpdateForFixedScale( fixedParameterValues );
    }
  }

  // This returns a string that is valid Python with no indentation to evaluate
  // the potential in three functions:
  // TreeLevelPotential( fv ), JustLoopCorrectedPotential( fv ), and
  // LoopAndThermallyCorrectedPotential( fv ).
   std::string FixedScaleOneLoopPotential::WriteActualPythonFunction() const
   {
     std::stringstream stringBuilder;
     stringBuilder << std::setprecision( 12 );
     stringBuilder
     << "fixedScaleInverseSquare = " << inverseRenormalizationScaleSquared
     << "\n"
     "\n"
     "def TreeLevelPotential( fv ):\n"
     "  return TreeLevelContribution( fv,\n"
     "                                fixedScaleLagrangianParameters )\n"
     "\n"
     "def JustLoopCorrectedPotential( fv ):\n"
     "  return ( TreeLevelContribution( fv,\n"
     "                                  fixedScaleLagrangianParameters )\n"
     "           + PolynomialLoopCorrections( fv,\n"
     "                                      fixedScaleLagrangianParameters )\n"
     "           + JustLoopCorrections( fv,\n"
     "                                  fixedScaleLagrangianParameters,\n"
     "                                  fixedScaleInverseSquare ) )\n"
     "\n"
     "def LoopAndThermallyCorrectedPotential( fv ):\n"
     "  return ( TreeLevelContribution( fv,\n"
     "                                  fixedScaleLagrangianParameters )\n"
     "           + PolynomialLoopCorrections( fv,\n"
     "                                      fixedScaleLagrangianParameters )\n"
     "           + LoopAndThermalCorrections( fv,\n"
     "                                       fixedScaleLagrangianParameters,\n"
     "                                        fixedScaleInverseSquare,\n"
     "                                        temperatureInverseSquare ) )\n";
     return stringBuilder.str();
   }

   // This is for debugging.
   std::string FixedScaleOneLoopPotential::PrintEvaluation(
                               std::vector< double > const& fieldConfiguration,
                                          double const temperatureValue ) const
   {
     std::stringstream stringBuilder;
     std::vector< DoubleVectorWithDouble > scalarMassesSquaredWithFactors;
     AddMassesSquaredWithMultiplicity( fieldConfiguration,
                                       scalarSquareMasses,
                                       scalarMassesSquaredWithFactors );
     std::vector< DoubleVectorWithDouble > fermionMassesSquaredWithFactors;
     AddMassesSquaredWithMultiplicity( fieldConfiguration,
                                       fermionSquareMasses,
                                       fermionMassesSquaredWithFactors );
     std::vector< DoubleVectorWithDouble > vectorMassesSquaredWithFactors;
     AddMassesSquaredWithMultiplicity( fieldConfiguration,
                                       vectorSquareMasses,
                                       vectorMassesSquaredWithFactors );

     stringBuilder << "FixedScaleOneLoopPotential::PrintEvaluation( "
     << FieldConfigurationAsMathematica( fieldConfiguration )
     << ", T =" << temperatureValue << " ) called."
     << std::endl << "scalarMassesSquaredWithFactors = {" << std::endl;
     for( std::vector< DoubleVectorWithDouble >::const_iterator
          massesWithFactor( scalarMassesSquaredWithFactors.begin() );
          massesWithFactor < scalarMassesSquaredWithFactors.end();
          ++massesWithFactor )
     {
       stringBuilder
       << "[ factor: " << massesWithFactor->second
       << ", masses-squared: { ";
       for( std::vector< double >::const_iterator
            massSquared( massesWithFactor->first.begin() );
            massSquared < massesWithFactor->first.end();
            ++massSquared )
       {
         stringBuilder << *massSquared << ",";
       }
       stringBuilder << " } ]" << std::endl;
     }
     stringBuilder << "}"
     << std::endl << "fermionMassesSquaredWithFactors = {" << std::endl;
     for( std::vector< DoubleVectorWithDouble >::const_iterator
          massesWithFactor( fermionMassesSquaredWithFactors.begin() );
          massesWithFactor < fermionMassesSquaredWithFactors.end();
          ++massesWithFactor )
     {
       stringBuilder
       << "[ factor: " << massesWithFactor->second
       << ", masses-squared: { ";
       for( std::vector< double >::const_iterator
            massSquared( massesWithFactor->first.begin() );
            massSquared < massesWithFactor->first.end();
            ++massSquared )
       {
         stringBuilder << *massSquared << ",";
       }
       stringBuilder << " } ]" << std::endl;
     }
     stringBuilder << "}"
     << std::endl << "vectorMassesSquaredWithFactors = {" << std::endl;
     for( std::vector< DoubleVectorWithDouble >::const_iterator
          massesWithFactor( vectorMassesSquaredWithFactors.begin() );
          massesWithFactor < vectorMassesSquaredWithFactors.end();
          ++massesWithFactor )
     {
       stringBuilder
       << "[ factor: " << massesWithFactor->second
       << ", masses-squared: { ";
       for( std::vector< double >::const_iterator
            massSquared( massesWithFactor->first.begin() );
            massSquared < massesWithFactor->first.end();
            ++massSquared )
       {
         stringBuilder << *massSquared << ",";
       }
       stringBuilder << " } ]" << std::endl;
     }
     stringBuilder << " }";
     stringBuilder << std::endl;
     stringBuilder << "treeLevelPotential = "
     << treeLevelPotential( fieldConfiguration );
     stringBuilder << std::endl;
     stringBuilder << "polynomialLoopCorrections = "
     << polynomialLoopCorrections( fieldConfiguration )
     << std::endl;
     stringBuilder << "LoopAndThermalCorrections = "
     << LoopAndThermalCorrections( scalarMassesSquaredWithFactors,
                                   fermionMassesSquaredWithFactors,
                                   vectorMassesSquaredWithFactors,
                                   inverseRenormalizationScaleSquared,
                                   temperatureValue );
     stringBuilder << std::endl;

     stringBuilder
     << "Total value = " << ( treeLevelPotential( fieldConfiguration )
                              + polynomialLoopCorrections( fieldConfiguration )
                   + LoopAndThermalCorrections( scalarMassesSquaredWithFactors,
                                               fermionMassesSquaredWithFactors,
                                                vectorMassesSquaredWithFactors,
                                            inverseRenormalizationScaleSquared,
                                                temperatureValue ) );
     return stringBuilder.str();
   }

} /* namespace VevaciousPlusPlus */
