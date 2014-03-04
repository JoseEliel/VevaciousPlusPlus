/*
 * PotentialFromPolynomialAndMasses.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POTENTIALFROMPOLYNOMIALANDMASSES_HPP_
#define POTENTIALFROMPOLYNOMIALANDMASSES_HPP_

#include "../StandardIncludes.hpp"
#include "HomotopyContinuationReadyPotential.hpp"
#include "PolynomialTerm.hpp"
#include "PolynomialSum.hpp"
#include "MassSquaredMatrix.hpp"
#include "RunningParameterManager.hpp"

namespace VevaciousPlusPlus
{
  class PotentialFromPolynomialAndMasses :
                                      public HomotopyContinuationReadyPotential
  {
  public:
    PotentialFromPolynomialAndMasses( std::string const& modelFilename );
    virtual
    ~PotentialFromPolynomialAndMasses();


    // This returns the energy density in GeV^4 of the potential for a state
    // strongly peaked around expectation values (in GeV) for the fields given
    // by the values of fieldConfiguration and temperature in GeV given by
    // temperatureValue.
    virtual double
    operator()( std::vector< double > const& fieldConfiguration,
                double const temperatureValue = 0.0 );

    // This updates all the parameters of the potential that are not field
    // values based on the values that appear in blocks in the SLHA format in
    // the file given by slhaFilename.
    virtual void UpdateParameters( std::string const& slhaFilename );

    // This returns the tree-level potential energy density.
    virtual double
    QuickApproximation( std::vector< double > const& fieldConfiguration,
                        double const temperatureValue = 0.0 );

    // This returns the square of the scale (in GeV^2) relevant to tunneling
    // between the given minima for this potential.
    virtual double
    ScaleSquaredRelevantToTunneling( PotentialMinimum const& falseVacuum,
                                     PotentialMinimum const& trueVacuum );

    // This evaluates the target system and places the values in
    // destinationVector.
    virtual void
    HomotopyContinuationSystemValues(
                           std::vector< double > fieldConfigurationWithScale,
                                    std::vector< double >& destinationVector );

    // This evaluates the derivatives of the target system and places the
    // values in destinationMatrix.
    virtual void
    HomotopyContinuationSystemGradients(
                             std::vector< double > fieldConfigurationWithScale,
                     std::vector< std::vector< double > >& destinationMatrix );


  protected:
    static std::string const digitChars;
    static std::string const dotAndDigits;
    static std::string const allowedVariableInitials;
    static std::string const allowedVariableChars;

    RunningParameterManager runningParameters;
    double renormalizationScaleSquared;
    double minimumRenormalizationScaleSquared;
    PolynomialSum treeLevelPotential;
    PolynomialSum polynomialLoopCorrections;
    std::vector< MassSquaredMatrix > massSquaredMatrices;
    double vectorMassCorrectionConstant;
    std::vector< PolynomialSum > polynomialGradient;
    std::vector< std::vector< PolynomialSum > > polynomialHessian;
    std::vector< PolynomialSum > scaleSlopeOfGradient;


    PotentialFromPolynomialAndMasses();

    // This puts all index brackets into a consistent form.
    std::string FormatVariable( std::string const& unformattedVariable ) const;

    // This interprets stringToParse as a sum of polynomial terms and sets
    // polynomialSum accordingly.
    void ParseSumOfPolynomialTerms( std::string const& stringToParse,
                                    PolynomialSum& polynomialSum );

    // This reads in a whole number or variable (including possible raising to
    // a power), applies the correct operation to polynomialTerm, and then
    // returns the position of the character just after the interpreted word.
    size_t
    PutNextNumberOrVariableIntoPolynomial( std::string const& stringToParse,
                                           size_t wordStart,
                                           PolynomialTerm& polynomialTerm );

    // This sets the renormalization scale and broadcasts it to the running
    // parameters.
    void UpdateRenormalizationScale(
                                    std::vector< double > fieldConfiguration,
                                        double const evaluationTemperature );

    // This evaluates the sum of corrections for a set of real scalar degrees
    // of freedom with masses-squared given by massesSquared at a temperature
    // given by evaluationTemperature.
    double ScalarBosonCorrection( std::vector< double > const& massesSquared,
                                  double const evaluationTemperature );

    // This evaluates the sum of corrections for a set of Weyl fermion degrees
    // of freedom with masses-squared given by massesSquared at a temperature
    // given by evaluationTemperature.
    double WeylFermionCorrection( std::vector< double > const& massesSquared,
                                  double const evaluationTemperature );

    // This evaluates the sum of corrections for a set of vector gauge boson
    // degrees of freedom (transverse modes) with masses-squared given by
    // massesSquared at a temperature given by evaluationTemperature.
    double GaugeBosonCorrection( std::vector< double > const& massesSquared,
                                 double const evaluationTemperature );

    // This prepares system of polynomials for the homotopy continuation as a
    // set of polynomials in the field variables with coefficients from the
    // polynomials of the tree-level potential and the polynomial loop
    // corrections. The coefficient of each polynomial term is fitted to a
    // polynomial of degree powerOfScale in the logarithm of the scale. After
    // this, the differentials of the system are derived from these polynomials
    // in the fields and the logarithm of the scale, and also for the
    // constraint relating the scale to the field variables.
    void PrepareHomotopyContinuationPolynomials( int const powerOfScale );
  };



  inline double PotentialFromPolynomialAndMasses::QuickApproximation(
                               std::vector< double > const& fieldConfiguration,
                                                double const temperatureValue )
  {
    UpdateRenormalizationScale( fieldConfiguration,
                                temperatureValue );
    return treeLevelPotential( fieldConfiguration );
  }

  // This puts all index brackets into a consistent form.
  inline std::string PotentialFromPolynomialAndMasses::FormatVariable(
                                 std::string const& unformattedVariable ) const
  {
    return RunningParameterManager::FormatVariable( unformattedVariable );
  }

} /* namespace VevaciousPlusPlus */
#endif /* POTENTIALFROMPOLYNOMIALANDMASSES_HPP_ */
