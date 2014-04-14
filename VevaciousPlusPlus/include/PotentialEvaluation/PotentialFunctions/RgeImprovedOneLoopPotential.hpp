/*
 * RgeImprovedOneLoopPotential.hpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef RGEIMPROVEDONELOOPPOTENTIAL_HPP_
#define RGEIMPROVEDONELOOPPOTENTIAL_HPP_

#include "../../StandardIncludes.hpp"
#include "PotentialFromPolynomialAndMasses.hpp"

namespace VevaciousPlusPlus
{

  class RgeImprovedOneLoopPotential : public PotentialFromPolynomialAndMasses
  {
  public:
    RgeImprovedOneLoopPotential( std::string const& modelFilename,
                            RunningParameterManager& runningParameterManager );
    RgeImprovedOneLoopPotential(
    PotentialFromPolynomialAndMasses const& potentialFromPolynomialAndMasses );
    virtual
    ~RgeImprovedOneLoopPotential();


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

    // This returns the square of the Euclidean distance between the two vacua.
    virtual double
    ScaleSquaredRelevantToTunneling( PotentialMinimum const& falseVacuum,
                                    PotentialMinimum const& trueVacuum ) const;

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
    double logarithmOfMinimumRenormalizationScale;
    double logarithmOfMaximumRenormalizationScale;

    // This sets dsbFieldValueInputs based on the SLHA file just read in.
    virtual void EvaluateDsbInputAndSetScale();

    // This returns the square of an appropriate renormalization scale.
    double RenormalizationScaleSquared(
                               std::vector< double > const& fieldConfiguration,
                                   double const temperatureValue = 0.0 ) const;

    // This appends the masses-squared and multiplicity from each
    // MassesSquaredFromMatrix in massSquaredMatrices to
    // massSquaredMatrices, with all functionoids evaluated at the natural
    // exponent of logarithmOfScale.
    void AddMassesSquaredWithMultiplicity(
                               std::vector< double > const& fieldConfiguration,
            std::vector< MassesSquaredCalculator* > const& massSquaredMatrices,
                std::vector< DoubleVectorWithDouble >& massesSquaredWithFactors,
                                         double const logarithmOfScale ) const;

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
    virtual void PrepareHomotopyContinuationPotentialPolynomial();

    // This should prepare homotopyContinuationStartSystem appropriately.
    virtual void PrepareHomotopyContinuationStartSystem();

    // This should prepare homotopyContinuationStartValues to be all the
    // solutions of homotopyContinuationStartSystem.
    virtual void PrepareHomotopyContinuationStartValues();
  };




  // This returns the square of an appropriate renormalization scale.
  inline double RgeImprovedOneLoopPotential::RenormalizationScaleSquared(
                               std::vector< double > const& fieldConfiguration,
                                          double const temperatureValue ) const
  {
    double renormalizationScaleSquared( squareOfMinimumRenormalizationScale
                                   + ( temperatureValue * temperatureValue ) );
    for( std::vector< double >::const_iterator
         whichField( fieldConfiguration.begin() );
         whichField < fieldConfiguration.end();
         ++whichField )
    {
      renormalizationScaleSquared += ( (*whichField) * (*whichField) );
    }
    return renormalizationScaleSquared;
  }

  // This appends the masses-squared and multiplicity from each
  // MassesSquaredFromMatrix in massSquaredMatrices to
  // massSquaredMatrices, with all functionoids evaluated at the last scale
  // which was used to update them.
  inline void
  RgeImprovedOneLoopPotential::AddMassesSquaredWithMultiplicity(
                               std::vector< double > const& fieldConfiguration,
            std::vector< MassesSquaredCalculator* > const& massSquaredMatrices,
               std::vector< DoubleVectorWithDouble >& massesSquaredWithFactors,
                                          double const logarithmOfScale ) const
  {
    for( std::vector< MassesSquaredCalculator* >::const_iterator
         whichMatrix( massSquaredMatrices.begin() );
         whichMatrix < massSquaredMatrices.end();
         ++whichMatrix )
    {
      massesSquaredWithFactors.push_back(
             std::make_pair( (*whichMatrix)->MassesSquared( fieldConfiguration,
                                                            logarithmOfScale ),
                                      (*whichMatrix)->MultiplicityFactor() ) );
    }
  }

  // This returns the tree-level potential energy density evaluated at the
  // correct scale.
  inline double RgeImprovedOneLoopPotential::QuickApproximation(
                               std::vector< double > const& fieldConfiguration,
                                                double const temperatureValue )
  {
    runningParameters.UpdateRunningParameters( sqrt(
                               RenormalizationScaleSquared( fieldConfiguration,
                                                        temperatureValue ) ) );
    return treeLevelPotential( fieldConfiguration );
  }

  // This returns the square of the Euclidean distance between the two vacua.
  inline double RgeImprovedOneLoopPotential::ScaleSquaredRelevantToTunneling(
                                           PotentialMinimum const& falseVacuum,
                                     PotentialMinimum const& trueVacuum ) const
  {
    return falseVacuum.SquareDistanceTo( trueVacuum );
  }

  // This prepares a system of polynomials for the homotopy continuation
  // based on the current SLHA input data. Each polynomial term in the
  // tree-level potential generates its derivatives in its fields with the
  // coefficients fitted to a polynomial in the logarithm of the
  // renormalization scale, and then also a polynomial relating the logarithm
  // of the renormalization scale to minimumRenormalizationScaleSquared and
  // the field values is also prepared.
  inline void
  RgeImprovedOneLoopPotential::PreparePolynomialHomotopyContinuation()
  {
    PrepareHomotopyContinuationPotentialPolynomial();
    PreparePolynomialGradient();
    // Now we add in the constraint on the log of the renormalization scale.
    PolynomialSum scaleConstraint;
    std::vector< PolynomialTerm >&
    constraintTerms( scaleConstraint.PolynomialTerms() );
    PolynomialTerm constraintTerm;
    constraintTerm.RaiseFieldPower( numberOfFields,
                                    2 );
    constraintTerms.push_back( constraintTerm );
    targetPolynomialGradient.push_back( scaleConstraint );
    PreparePolynomialHessian();
    PrepareHomotopyContinuationStartSystem();
    PrepareHomotopyContinuationStartValues();
  }

} /* namespace VevaciousPlusPlus */
#endif /* RGEIMPROVEDONELOOPPOTENTIAL_HPP_ */
