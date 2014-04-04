/*
 * RgeImprovedOneLoopPotential.hpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef RGEIMPROVEDONELOOPPOTENTIAL_HPP_
#define RGEIMPROVEDONELOOPPOTENTIAL_HPP_

#include "../StandardIncludes.hpp"
#include "PotentialFromPolynomialAndMasses.hpp"

namespace VevaciousPlusPlus
{

  class RgeImprovedOneLoopPotential : public PotentialFromPolynomialAndMasses
  {
  public:
    RgeImprovedOneLoopPotential( std::string const& modelFilename );
    virtual
    ~RgeImprovedOneLoopPotential();


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

    // This sets the renormalization scale and broadcasts it to the running
    // parameters.
    void UpdateRenormalizationScale(
                               std::vector< double > const& fieldConfiguration,
                                     double const evaluationTemperature );

    // This prepares a system of polynomials for the homotopy continuation
    // based on the current SLHA input data. Each polynomial term in the
    // tree-level potential generates its derivatives in its fields with the
    // coefficients fitted to a polynomial in the logarithm of the
    // renormalization scale, and then also a polynomial relating the logarithm
    // of the renormalization scale to minimumRenormalizationScaleSquared and
    // the field values is also prepared.
    virtual void PreparePolynomialHomotopyContinuation();

    // This should prepare homotopyContinuationPotentialPolynomial
    // appropriately.
    virtual void PrepareHomotopyContinuationPotentialPolynomial() = 0;
  };




  // This returns the energy density in GeV^4 of the potential for a state
  // strongly peaked around expectation values (in GeV) for the fields given
  // by the values of fieldConfiguration and temperature in GeV given by
  // temperatureValue.
  inline double RgeImprovedOneLoopPotential::operator()(
                               std::vector< double > const& fieldConfiguration,
                                                double const temperatureValue )
  {
    UpdateRenormalizationScale( fieldConfiguration,
                                temperatureValue );
    return LoopAndThermallyCorrectedPotential( fieldConfiguration,
                                               temperatureValue );
  }

  // This returns the tree-level potential energy density evaluated at the
  // correct scale.
  inline double RgeImprovedOneLoopPotential::QuickApproximation(
                               std::vector< double > const& fieldConfiguration,
                                                double const temperatureValue )
  {
    UpdateRenormalizationScale( fieldConfiguration,
                                temperatureValue );
    return treeLevelPotential( fieldConfiguration );
  }

  // This returns the square of the Euclidean distance between the two vacua.
  inline double RgeImprovedOneLoopPotential::ScaleSquaredRelevantToTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    return falseVacuum.SquareDistanceTo( trueVacuum );
  }

  // This sets dsbFieldValueInputs based on the SLHA file just read in.
  inline void RgeImprovedOneLoopPotential::EvaluateDsbInputAndSetScale()
  {
    UpdateRenormalizationScale( fieldOrigin,
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
    << "RgeImprovedOneLoopPotential::EvaluateDsbInputAndSetScale() set DSB"
    << " field values to:" << std::endl;
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      std::cout << fieldNames[ fieldIndex ] << " -> "
      << dsbFieldValueInputs[ fieldIndex ] << std::endl;
    }
    std::cout << std::endl;/**/
  }

  // This sets the renormalization scale and broadcasts it to the running
  // parameters.
  inline void RgeImprovedOneLoopPotential::UpdateRenormalizationScale(
                               std::vector< double > const& fieldConfiguration,
                                           double const evaluationTemperature )
  {
    renormalizationScaleSquared
    = ( minimumRenormalizationScaleSquared
        + ( evaluationTemperature * evaluationTemperature ) );
    for( std::vector< double >::const_iterator
         whichField( fieldConfiguration.begin() );
         whichField < fieldConfiguration.end();
         ++whichField )
    {
      renormalizationScaleSquared += ( (*whichField) * (*whichField) );
    }
    runningParameters.UpdateRunningParameters(
                                         sqrt( renormalizationScaleSquared ) );
  }

  // This fills polynomialGradient, homotopyContinuationStartSystem, and
  // polynomialHessian appropriately.
  inline void
  RgeImprovedOneLoopPotential::PreparePolynomialHomotopyContinuation()
  {
    PrepareHomotopyContinuationPotentialPolynomial();
    PrepareHomotopyContinuationGradient();
    // Now we add in the constraint on the log of the renormalization scale.
    PolynomialSum scaleConstraint;
    std::vector< PolynomialTerm >&
    constraintTerms( scaleConstraint.PolynomialTerms() );
    PolynomialTerm constraintTerm;
    constraintTerm.RaiseFieldPower( numberOfFields,
                                    2 );
    constraintTerms.push_back( constraintTerm );
    polynomialGradient.push_back( scaleConstraint );
    PrepareHomotopyContinuationHessian();
    PrepareHomotopyContinuationStartSystem();
    PrepareHomotopyContinuationStartValues();
  }

} /* namespace VevaciousPlusPlus */
#endif /* RGEIMPROVEDONELOOPPOTENTIAL_HPP_ */
