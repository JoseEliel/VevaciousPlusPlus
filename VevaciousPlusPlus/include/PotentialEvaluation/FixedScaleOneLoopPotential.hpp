/*
 * FixedScaleOneLoopPotential.hpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef FIXEDSCALEONELOOPPOTENTIAL_HPP_
#define FIXEDSCALEONELOOPPOTENTIAL_HPP_

#include "../StandardIncludes.hpp"
#include "PotentialFromPolynomialAndMasses.hpp"

namespace VevaciousPlusPlus
{

  class FixedScaleOneLoopPotential : public PotentialFromPolynomialAndMasses
  {
  public:
    FixedScaleOneLoopPotential( std::string const& modelFilename );
    virtual
    ~FixedScaleOneLoopPotential();


    // This returns the energy density in GeV^4 of the potential for a state
    // strongly peaked around expectation values (in GeV) for the fields given
    // by the values of fieldConfiguration and temperature in GeV given by
    // temperatureValue.
    virtual double
    operator()( std::vector< double > const& fieldConfiguration,
                double const temperatureValue = 0.0 );

    // This returns the tree-level potential energy density evaluated at the
    // correct scale.
    virtual double
    QuickApproximation( std::vector< double > const& fieldConfiguration,
                        double const temperatureValue = 0.0 );

    // This returns the square of the Euclidean distance between the two vacua.
    virtual double
    ScaleSquaredRelevantToTunneling( PotentialMinimum const& falseVacuum,
                                     PotentialMinimum const& trueVacuum );

    // This evaluates the target system and places the values in
    // destinationVector.
    virtual void
    HomotopyContinuationSystemValues(
                                   std::vector< double > solutionConfiguration,
                                    std::vector< double >& destinationVector );

    // This evaluates the derivatives of the target system and places the
    // values in destinationMatrix.
    virtual void
    HomotopyContinuationSystemGradients(
                                   std::vector< double > solutionConfiguration,
                     std::vector< std::vector< double > >& destinationMatrix );


  protected:
    // This sets dsbFieldValueInputs based on the SLHA file just read in.
    virtual void EvaluateDsbInputAndSetScale();

    // This prepares a system of polynomials for the homotopy continuation
    // based on the current SLHA input data. Each polynomial term in the
    // tree-level potential generates its derivatives in its fields with the
    // coefficients fitted to a polynomial in the logarithm of the
    // renormalization scale, and then also a polynomial relating the logarithm
    // of the renormalization scale to minimumRenormalizationScaleSquared and
    // the field values is also prepared.
    virtual void PreparePolynomialHomotopyContinuation();
  };




  // This returns the energy density in GeV^4 of the potential for a state
  // strongly peaked around expectation values (in GeV) for the fields given
  // by the values of fieldConfiguration and temperature in GeV given by
  // temperatureValue.
  inline double FixedScaleOneLoopPotential::operator()(
                               std::vector< double > const& fieldConfiguration,
                                                double const temperatureValue )
  {
    return LoopAndThermallyCorrectedPotential( fieldConfiguration,
                                               temperatureValue );
  }

  // This returns the tree-level potential energy density evaluated at the
  // correct scale.
  inline double FixedScaleOneLoopPotential::QuickApproximation(
                               std::vector< double > const& fieldConfiguration,
                                                double const temperatureValue )
  {
    return treeLevelPotential( fieldConfiguration );
  }

  // This returns the square of the Euclidean distance between the two vacua.
  inline double FixedScaleOneLoopPotential::ScaleSquaredRelevantToTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    return renormalizationScaleSquared;
  }

  // This sets dsbFieldValueInputs based on the SLHA file just read in.
  inline void FixedScaleOneLoopPotential::EvaluateDsbInputAndSetScale()
  {
    renormalizationScaleSquared = runningParameters.LowestBlockScale();
    renormalizationScaleSquared *= renormalizationScaleSquared;
    renormalizationScaleSquared = std::max( renormalizationScaleSquared,
                                          minimumRenormalizationScaleSquared );
    runningParameters.UpdateRunningParameters(
                                         sqrt( renormalizationScaleSquared ) );
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
    << "FixedScaleOneLoopPotential::EvaluateDsbInputAndSetScale() set"
    << " renormalizationScaleSquared to " << renormalizationScaleSquared
    << " and DSB field values to:" << std::endl;
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      std::cout << fieldNames[ fieldIndex ] << " -> "
      << dsbFieldValueInputs[ fieldIndex ] << std::endl;
    }
    std::cout << std::endl;/**/
  }

} /* namespace VevaciousPlusPlus */
#endif /* FIXEDSCALEONELOOPPOTENTIAL_HPP_ */
