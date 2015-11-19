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
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "FixedScaleOneLoopPotential::respondToObservedSignal() called.";
    std::cout << std::endl;*/

    renormalizationScale
    = lagrangianParameterManager.AppropriateFixedScaleForParameterPoint();
    inverseRenormalizationScaleSquared
    = ( 1.0 / ( renormalizationScale * renormalizationScale ) );
    std::vector< double > const
    fixedParameterValues( lagrangianParameterManager.ParameterValues(
                                               log( renormalizationScale ) ) );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "FixedScaleOneLoopPotential::respondToObservedSignal(): before updating"
    << "treeLevelPotential, treeLevelPotential = "
    << treeLevelPotential.AsDebuggingString();
    std::cout << std::endl;*/

    treeLevelPotential.UpdateForFixedScale( fixedParameterValues );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "fixedParameterValues = { ";
    for( std::vector< double >::const_iterator
         parameterValue( fixedParameterValues.begin() );
         parameterValue < fixedParameterValues.end();
         ++parameterValue )
    {
      std::cout << *parameterValue << ", ";
    }
    std::cout << "}, now treeLevelPotential = "
    << treeLevelPotential.AsDebuggingString();
    std::cout << std::endl;*/

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

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "renormalizationScale = " << renormalizationScale
    << ", inverseRenormalizationScaleSquared = "
    << inverseRenormalizationScaleSquared
    << ", log( renormalizationScale ) = "
    << log( renormalizationScale )
    << std::endl
    << "fixedParameterValues = { ";
    for( std::vector< double >::const_iterator
         parameterValue( fixedParameterValues.begin() );
         parameterValue < fixedParameterValues.end();
         ++parameterValue )
    {
      std::cout << *parameterValue << ",";
    }
    std::cout << " }" << std::endl;
    std::vector< double > fieldConfiguration( FieldValuesOrigin() );
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
    std::cout << "For origin:"
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
    std::cout << "}";
    std::cout << std::endl;*/
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
     "                                        fixedScaleInverseSquare\n"
     "                                        temperatureInverseSquare ) )\n";
     return stringBuilder.str();
   }

} /* namespace VevaciousPlusPlus */
