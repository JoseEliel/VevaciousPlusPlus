/*
 * TreeLevelPotential.hpp
 *
 *  Created on: Oct 12, 2020
 *      Author: José Eliel Camargo-Molina (Eliel@camargo-molina.com) 
 */

#ifndef TreeLevelPotential_HPP_
#define TreeLevelPotential_HPP_

#include "PotentialFromPolynomialWithMasses.hpp"
#include "LHPC/Utilities/BasicObserverPattern.hpp"
#include <string>
#include "LagrangianParameterManagement/LagrangianParameterManager.hpp"
#include <vector>
#include "PotentialMinimization/PotentialMinimum.hpp"
#include <cmath>
#include "PotentialEvaluation/MassesSquaredCalculators/RealMassesSquaredMatrix.hpp"
#include "PotentialEvaluation/MassesSquaredCalculators/SymmetricComplexMassMatrix.hpp"
#include "PotentialEvaluation/MassesSquaredCalculators/ComplexMassSquaredMatrix.hpp"
#include <sstream>
#include <iomanip>

namespace VevaciousPlusPlus
{

  class TreeLevelPotential : public PotentialFromPolynomialWithMasses,
                                     public LHPC::BasicObserver
  {
  public:
    TreeLevelPotential( std::string const& modelFilename,
                               double const assumedPositiveOrNegativeTolerance,
                      LagrangianParameterManager& lagrangianParameterManager );
    TreeLevelPotential(
                    PotentialFromPolynomialWithMasses const& potentialToCopy );
    virtual ~TreeLevelPotential();


    // This returns the energy density in GeV^4 of the potential for a state
    // strongly peaked around expectation values (in GeV) for the fields given
    // by the values of fieldConfiguration and temperature in GeV given by
    // temperatureValue.
    virtual double
    operator()( std::vector< double > const& fieldConfiguration,
                double const temperatureValue = 0.0 ) const;

    // This returns the square of the current renormalization scale.
    virtual double
    ScaleSquaredRelevantToTunneling( PotentialMinimum const& falseVacuum,
                                     PotentialMinimum const& trueVacuum ) const
    { return ( renormalizationScale * renormalizationScale ); }

    // This updates the scale used for the loop corrections based on the
    // appropriate scale from lagrangianParameterManager, and updates all
    // components used to evaluate the potential to use the Lagrangian
    // parameters evaluated at that scale.
    virtual void RespondToObservedSignal();


    // This is for debugging.
    std::string
    PrintEvaluation( std::vector< double > const& fieldConfiguration,
                     double const temperatureValue = 0.0 ) const;

    // Not Implemented, no CosmoTransitions support for tree level potential. 
    virtual std::string WriteActualPythonFunction() const;
    virtual void WriteAsPython( std::string const& pythonFilename ) const;


  protected:
    double renormalizationScale;
    double inverseRenormalizationScaleSquared;
  };




  // This returns the energy density in GeV^4 of the potential for a state
  // strongly peaked around expectation values (in GeV) for the fields given
  // by the values of fieldConfiguration and temperature in GeV given by
  // temperatureValue.
  inline double TreeLevelPotential::operator()(
                               std::vector< double > const& fieldConfiguration,
                                          double const temperatureValue ) const
  {
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
    return ( treeLevelPotential( fieldConfiguration )
             + JustThermalCorrections( scalarMassesSquaredWithFactors,
                                          fermionMassesSquaredWithFactors,
                                          vectorMassesSquaredWithFactors,
                                          inverseRenormalizationScaleSquared,
                                          temperatureValue ) );
  }

} /* namespace VevaciousPlusPlus */

#endif /* TreeLevelPotential_HPP_ */
