/*
 * TreeLevelPotential.cpp
 *
 *  Created on: Oct 12, 2020
 *      Author: José Eliel Camargo-Molina (Eliel@camargo-molina.com) 
 */

#include "PotentialEvaluation/PotentialFunctions/TreeLevelPotential.hpp"

namespace VevaciousPlusPlus
{

  TreeLevelPotential::TreeLevelPotential(
                                              std::string const& modelFilename,
                               double const assumedPositiveOrNegativeTolerance,
                     LagrangianParameterManager& lagrangianParameterManager ) :
    PotentialFromPolynomialWithMasses( modelFilename,
                                       assumedPositiveOrNegativeTolerance,
                                       lagrangianParameterManager ),
    LHPC::BasicObserver(),
    renormalizationScale( -1.0 ),
    inverseRenormalizationScaleSquared( -1.0 )
  {
    lagrangianParameterManager.RegisterObserver( this );
  }

  TreeLevelPotential::TreeLevelPotential(
                   PotentialFromPolynomialWithMasses const& potentialToCopy ) :
    PotentialFromPolynomialWithMasses( potentialToCopy ),
    LHPC::BasicObserver(),
    renormalizationScale( -1.0 ),
    inverseRenormalizationScaleSquared( -1.0 )
  {
    lagrangianParameterManager.RegisterObserver( this );
  }

  TreeLevelPotential::~TreeLevelPotential()
  {
    // This does nothing.
  }

  void TreeLevelPotential::WriteAsPython(
                                      std::string const& pythonFilename ) const
  {
    throw std::runtime_error(
            "PotentialFunction::WriteAsPython(..) does not work with TreeLevelPotential."
            "CosmoTransitions integration is not implemented for this potential class."
            "Use FixedScaleOneLoopPotential Instead." );
  }

  std::string TreeLevelPotential::WriteActualPythonFunction() const
   {
     throw std::runtime_error(
            "PotentialFunction::WriteActualPythonFunction(..) does not work with TreeLevelPotential."
            "CosmoTransitions integration is not implemented for this potential class."
            "Use FixedScaleOneLoopPotential Instead." );

   }



  // This updates the scale used for the loop corrections based on the
  // appropriate scale from lagrangianParameterManager, and updates all
  // components used to evaluate the potential to use the Lagrangian
  // parameters evaluated at that scale.
  void TreeLevelPotential::RespondToObservedSignal()
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

} /* namespace VevaciousPlusPlus */
