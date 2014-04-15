/*
 * FixedScaleOneLoopPotential.hpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef FIXEDSCALEONELOOPPOTENTIAL_HPP_
#define FIXEDSCALEONELOOPPOTENTIAL_HPP_

#include "../../StandardIncludes.hpp"
#include "PotentialFromPolynomialAndMasses.hpp"

namespace VevaciousPlusPlus
{

  class FixedScaleOneLoopPotential : public PotentialFromPolynomialAndMasses
  {
  public:
    FixedScaleOneLoopPotential( std::string const& modelFilename,
                            RunningParameterManager& runningParameterManager );
    FixedScaleOneLoopPotential(
          PotentialFromPolynomialAndMasses& potentialFromPolynomialAndMasses );
    virtual
    ~FixedScaleOneLoopPotential();


    // This returns the energy density in GeV^4 of the potential for a state
    // strongly peaked around expectation values (in GeV) for the fields given
    // by the values of fieldConfiguration and temperature in GeV given by
    // temperatureValue.
    virtual double
    operator()( std::vector< double > const& fieldConfiguration,
                double const temperatureValue = 0.0 ) const;

    // This returns the tree-level potential energy density evaluated at the
    // correct scale.
    virtual double
    QuickApproximation( std::vector< double > const& fieldConfiguration,
                        double const temperatureValue = 0.0 );

    // This returns the square of the renormalization scale.
    virtual double
    ScaleSquaredRelevantToTunneling( PotentialMinimum const& falseVacuum,
                                     PotentialMinimum const& trueVacuum ) const
    { return ( currentMinimumRenormalizationScale
               * currentMinimumRenormalizationScale ); }

    // This should return a vector of field values corresponding to the field
    // configuration as it should be passed to operator() for evaluating the
    // potential, given a vector of values that solve this instance's homotopy
    // continuation system. It should return an empty vector if
    // homotopyContinuatioConfiguration does not correspond to a valid field
    // configuration. (For example, RgeImprovedOneLoopPotential uses the
    // logarithm of the renormalization scale as an extra variable in the
    // homotopy continuation, so homotopyContinuatioConfiguration actually has
    // an extra entry compared to a valid field configuration. Also, it may
    // be that homotopyContinuatioConfiguration did not correspond to a valid
    // solution where the scale is close to the Euclidean length of the field
    // configuration, so it would not correspond to a valid solution.)
    virtual std::vector< double > ValidFieldsFromHomotopyContinuation(
                std::vector< double > homotopyContinuatioConfiguration ) const;


  protected:
    double inverseRenormalizationScaleSquared;

    // This sets dsbFieldValueInputs based on the SLHA file just read in.
    virtual void EvaluateDsbInputAndSetScale();

    // This should prepare homotopyContinuationStartSystem appropriately.
    virtual void PrepareHomotopyContinuationStartSystem();

    // This should prepare homotopyContinuationStartValues to be all the
    // solutions of homotopyContinuationStartSystem.
    virtual void PrepareHomotopyContinuationStartValues();
  };




  // This returns the tree-level potential energy density evaluated at the
  // correct scale.
  inline double FixedScaleOneLoopPotential::QuickApproximation(
                               std::vector< double > const& fieldConfiguration,
                                                double const temperatureValue )
  {
    return treeLevelPotential( fieldConfiguration );
  }

  // FixedScaleOneLoopPotential prepares the homotopy continuation polynomials
  // in 1-to-1 correspondence with how they should be sent to operator().
  inline std::vector< double >
  FixedScaleOneLoopPotential::ValidFieldsFromHomotopyContinuation(
                 std::vector< double > homotopyContinuatioConfiguration ) const
  {
    return homotopyContinuatioConfiguration;
  }

  // This sets dsbFieldValueInputs based on the SLHA file just read in.
  inline void FixedScaleOneLoopPotential::EvaluateDsbInputAndSetScale()
  {
    currentMaximumRenormalizationScale = runningParameters.HighestBlockScale();
    currentMinimumRenormalizationScale = runningParameters.LowestBlockScale();
    squareOfMinimumRenormalizationScale = ( currentMinimumRenormalizationScale
                                        * currentMinimumRenormalizationScale );
    inverseRenormalizationScaleSquared
    = ( 1.0 / squareOfMinimumRenormalizationScale );
    runningParameters.UpdateRunningParameters(
                                          currentMinimumRenormalizationScale );
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
#endif /* FIXEDSCALEONELOOPPOTENTIAL_HPP_ */
