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
#include "LagrangianParameterManagement/ParameterUpdatePropagator.hpp"

namespace VevaciousPlusPlus
{

  class FixedScaleOneLoopPotential : public PotentialFromPolynomialWithMasses,
                                     public ParameterUpdatePropagator
  {
  public:
    FixedScaleOneLoopPotential( std::string const& modelFilename,
                               double const assumedPositiveOrNegativeTolerance,
                        ParameterUpdatePropagator& parameterUpdatePropagator );
    FixedScaleOneLoopPotential(
                      PotentialFromPolynomialWithMasses const& potentialToCopy,
                        ParameterUpdatePropagator& parameterUpdatePropagator );
    virtual ~FixedScaleOneLoopPotential();


    // This returns the energy density in GeV^4 of the potential for a state
    // strongly peaked around expectation values (in GeV) for the fields given
    // by the values of fieldConfiguration and temperature in GeV given by
    // temperatureValue.
    virtual double
    operator()( std::vector< double > const& fieldConfiguration,
                double const temperatureValue = 0.0 ) const;

    // This returns the tree-level evaluation, ignoring the temperature.
    virtual double
    QuickApproximation( std::vector< double > const& fieldConfiguration,
                        double const temperatureValue = 0.0 )
    { return treeLevelPotential( fieldConfiguration ); }

    // This returns the square of the current renormalization scale.
    virtual double
    ScaleSquaredRelevantToTunneling( PotentialMinimum const& falseVacuum,
                                PotentialMinimum const& trueVacuum ) const
    { return ( renormalizationScale * renormalizationScale ); }

    // This updates the scale used for the loop corrections based on the
    // appropriate scale from lagrangianParameterManager.
    virtual void UpdateSelfForNewParameterPoint(
                LagrangianParameterManager const& lagrangianParameterManager );


  protected:
    double renormalizationScale;
    double inverseRenormalizationScaleSquared;
  };





  // This updates the scale used for the loop corrections based on the
  // appropriate scale from lagrangianParameterManager.
  inline void FixedScaleOneLoopPotential::UpdateSelfForNewParameterPoint(
                 LagrangianParameterManager const& lagrangianParameterManager )
  {
    renormalizationScale
    = lagrangianParameterManager.AppropriateFixedScaleForParameterPoint();
    inverseRenormalizationScaleSquared
    = ( 1.0 / ( renormalizationScale * renormalizationScale ) );
  }

} /* namespace VevaciousPlusPlus */

#endif /* FIXEDSCALEONELOOPPOTENTIAL_HPP_ */
