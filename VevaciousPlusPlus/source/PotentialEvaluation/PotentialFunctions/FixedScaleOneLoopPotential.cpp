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
                       ParameterUpdatePropagator& parameterUpdatePropagator ) :
    PotentialFromPolynomialWithMasses( modelFilename,
                                       assumedPositiveOrNegativeTolerance,
                   parameterUpdatePropagator.GetLagrangianParameterManager() ),
    ParameterUpdatePropagator( parameterUpdatePropagator ),
    renormalizationScale( -1.0 ),
    inverseRenormalizationScaleSquared( -1.0 )
  {
    // This constructor is just an initialization list.
  }

  FixedScaleOneLoopPotential::FixedScaleOneLoopPotential(
                    PotentialFromPolynomialWithMasses const& potentialToCopy,
                      ParameterUpdatePropagator& parameterUpdatePropagator ) :
    PotentialFromPolynomialWithMasses( potentialToCopy ),
    ParameterUpdatePropagator( parameterUpdatePropagator ),
    renormalizationScale( -1.0 ),
    inverseRenormalizationScaleSquared( -1.0 )
  {
    // This constructor is just an initialization list.
  }

  FixedScaleOneLoopPotential::~FixedScaleOneLoopPotential()
  {
    // This does nothing.
  }


  // This updates the scale used for the loop corrections based on the
  // appropriate scale from lagrangianParameterManager, and updates all
  // components used to evaluate the potential to use the Lagrangian
  // parameters evaluated at that scale.
  void FixedScaleOneLoopPotential::UpdateSelfForNewParameterPoint(
                 LagrangianParameterManager const& lagrangianParameterManager )
  {
    renormalizationScale
    = lagrangianParameterManager.AppropriateFixedScaleForParameterPoint();
    inverseRenormalizationScaleSquared
    = ( 1.0 / ( renormalizationScale * renormalizationScale ) );
    std::vector< double > const
    fixedParameterValues( lagrangianParameterManager.ParameterValues(
                                               log( renormalizationScale ) ) );
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
     "                                        fixedScaleInverseSquare\n"
     "                                        temperatureInverseSquare ) )\n";
     return stringBuilder.str();
   }

} /* namespace VevaciousPlusPlus */
