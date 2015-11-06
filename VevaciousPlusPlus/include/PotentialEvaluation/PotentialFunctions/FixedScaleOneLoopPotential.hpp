/*
 * FixedScaleOneLoopPotential.hpp
 *
 *  Created on: Nov 4, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef FIXEDSCALEONELOOPPOTENTIAL_HPP_
#define FIXEDSCALEONELOOPPOTENTIAL_HPP_

#include "CommonIncludes.hpp"
#include "PotentialFromPolynomialWithMasses.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"

namespace VevaciousPlusPlus
{

  class FixedScaleOneLoopPotential : public PotentialFromPolynomialWithMasses,
                     public BOL::PushedToObserver< LagrangianParameterManager >
  {
  public:
    FixedScaleOneLoopPotential( std::string const& modelFilename,
                               double const assumedPositiveOrNegativeTolerance,
                      LagrangianParameterManager& lagrangianParameterManager );
    FixedScaleOneLoopPotential(
                    PotentialFromPolynomialWithMasses const& potentialToCopy );
    virtual ~FixedScaleOneLoopPotential();


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
    virtual void respondToPush(
                LagrangianParameterManager const& lagrangianParameterManager );

    // This returns a string that is valid Python with no indentation to
    // evaluate the potential in three functions:
    // TreeLevelPotential( fv ), JustLoopCorrectedPotential( fv ), and
    // LoopAndThermallyCorrectedPotential( fv ).
    virtual std::string WriteActualPythonFunction() const;


  protected:
    double renormalizationScale;
    double inverseRenormalizationScaleSquared;
  };





  // This returns the energy density in GeV^4 of the potential for a state
  // strongly peaked around expectation values (in GeV) for the fields given
  // by the values of fieldConfiguration and temperature in GeV given by
  // temperatureValue.
  inline double FixedScaleOneLoopPotential::operator()(
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
             + polynomialLoopCorrections( fieldConfiguration )
             + LoopAndThermalCorrections( scalarMassesSquaredWithFactors,
                                          fermionMassesSquaredWithFactors,
                                          vectorMassesSquaredWithFactors,
                                          inverseRenormalizationScaleSquared,
                                          temperatureValue ) );
  }

} /* namespace VevaciousPlusPlus */

#endif /* FIXEDSCALEONELOOPPOTENTIAL_HPP_ */
