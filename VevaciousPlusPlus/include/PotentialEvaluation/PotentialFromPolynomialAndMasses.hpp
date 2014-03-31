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
#include "RealMassesSquaredMatrix.hpp"
#include "ComplexMassMatrix.hpp"
#include "ComplexMassSquaredMatrix.hpp"
#include "MassesSquaredFromPolynomials.hpp"
#include "RunningParameterManager.hpp"
#include "ThermalFunctions.hpp"

namespace VevaciousPlusPlus
{
  class PotentialFromPolynomialAndMasses :
                                      public HomotopyContinuationReadyPotential
  {
  public:
    PotentialFromPolynomialAndMasses( std::string const& modelFilename );
    virtual
    ~PotentialFromPolynomialAndMasses();


    // This updates all the parameters of the potential that are not field
    // values based on the values that appear in blocks in the SLHA format in
    // the file given by slhaFilename.
    virtual void UpdateParameters( std::string const& slhaFilename );


  protected:
    static std::string const digitChars;
    static std::string const dotAndDigits;
    static std::string const allowedVariableInitials;
    static std::string const allowedVariableChars;
    static double const loopFactor;
    static double const thermalFactor;

    RunningParameterManager runningParameters;
    std::vector< PolynomialSum > dsbFieldValuePolynomials;
    double renormalizationScaleSquared;
    double minimumRenormalizationScaleSquared;
    PolynomialSum treeLevelPotential;
    PolynomialSum polynomialLoopCorrections;
    std::vector< RealMassesSquaredMatrix > scalarSquareMasses;
    std::vector< ComplexMassMatrix > fermionMasses;
    std::vector< ComplexMassSquaredMatrix > fermionMassSquareds;
    std::vector< RealMassesSquaredMatrix > vectorSquareMasses;
    double vectorMassCorrectionConstant;
    std::vector< PolynomialSum > polynomialGradient;
    std::vector< std::vector< PolynomialSum > > polynomialHessian;
    std::vector< PolynomialSum > scaleSlopeOfGradient;


    PotentialFromPolynomialAndMasses();


    // This should set dsbFieldValueInputs based on the SLHA file just read in.
    virtual void EvaluateDsbInputAndSetScale() = 0;

    // This evaluates the one-loop potential with thermal corrections assuming
    // that the scale has been set correctly.
    inline double LoopAndThermallyCorrectedPotential(
                               std::vector< double > const& fieldConfiguration,
                                               double const temperatureValue );

    // This puts all index brackets into a consistent form.
    std::string FormatVariable( std::string const& unformattedVariable ) const;

    // This interprets stringToParse as a sum of complex polynomial terms and
    // sets polynomialSum accordingly.
    void ParseSumOfPolynomialTerms( std::string const& stringToParse,
                   std::pair< PolynomialSum, PolynomialSum >& polynomialSums );

    // This interprets stringToParse as a sum of real polynomial terms and sets
    // polynomialSum accordingly.
    void ParseSumOfPolynomialTerms( std::string const& stringToParse,
                                    PolynomialSum& polynomialSum );

    // This reads in a whole number or variable (including possible raising to
    // a power), applies the correct operation to polynomialTerm, and then
    // returns the position of the character just after the interpreted word.
    // If there was a factor of "i", "I", "j", or "J", imaginaryTerm is set to
    // true.
    size_t
    PutNextNumberOrVariableIntoPolynomial( std::string const& stringToParse,
                                           size_t wordStart,
                                           PolynomialTerm& polynomialTerm,
                                           bool& imaginaryTerm );

    // This evaluates the sum of corrections for the real scalar degrees
    // of freedom with masses-squared given by scalarSquareMasses evaluated for
    // the field configuration given by fieldConfiguration at a temperature T
    // given by inverseTemperatureSquared for T^(-2) and temperatureFourthed
    // for T^4.
    double
    ScalarBosonCorrections( std::vector< double > const& fieldConfiguration,
                            double const inverseTemperatureSquared,
                            double const temperatureFourthed );

    // This evaluates the sum of corrections for a set of Weyl fermion degrees
    // of freedom with masses-squared given by fermionMasses evaluated for
    // the field configuration given by fieldConfiguration at a temperature T
    // given by inverseTemperatureSquared for T^(-2) and temperatureFourthed
    // for T^4.
    double
    WeylFermionCorrections( std::vector< double > const& fieldConfiguration,
                            double const inverseTemperatureSquared,
                            double const temperatureFourthed );

    // This evaluates the sum of corrections for the real scalar degrees
    // of freedom with masses-squared given by vectorSquareMasses evaluated for
    // the field configuration given by fieldConfiguration at a temperature T
    // given by inverseTemperatureSquared for T^(-2) and temperatureFourthed
    // for T^4.
    double
    GaugeBosonCorrections( std::vector< double > const& fieldConfiguration,
                           double const inverseTemperatureSquared,
                           double const temperatureFourthed );

    // This returns the J function thermal correction for a bosonic degree of
    // freedom based on a lookup table.
    double bosonThermalFunction(
                         double const massSquaredOverTemperatureSquared ) const
    { return ThermalFunctions::BosonicJ( massSquaredOverTemperatureSquared ); }

    // This returns the J function thermal correction for a fermionic degree of
    // freedom based on a lookup table.
    double fermionThermalFunction(
                         double const massSquaredOverTemperatureSquared ) const
    { return
      ThermalFunctions::FermionicJ( massSquaredOverTemperatureSquared ); }

    // This should prepare a system of polynomials for the homotopy
    // continuation based on the current SLHA input data.
    virtual void PrepareHomotopyContinuationPolynomials() = 0;
  };




  // This updates all the parameters of the potential that are not field
  // values based on the values that appear in blocks in the SLHA format in
  // the file given by slhaFilename.
  inline void PotentialFromPolynomialAndMasses::UpdateParameters(
                                              std::string const& slhaFilename )
  {
    runningParameters.UpdateSlhaParameters( slhaFilename );
    PrepareHomotopyContinuationPolynomials();
    EvaluateDsbInputAndSetScale();
  }

  // This evaluates the one-loop potential with thermal corrections assuming
  // that the scale has been set correctly.
  inline double
  PotentialFromPolynomialAndMasses::LoopAndThermallyCorrectedPotential(
                             std::vector< double > const& fieldConfiguration,
                                                double const temperatureValue )
  {
    double inverseTemperatureSquared( 1.0 );
    if( temperatureValue > 0.0 )
    {
      inverseTemperatureSquared
      = ( 1.0 / ( temperatureValue * temperatureValue ) );
    }
    double const temperatureFourthed( temperatureValue * temperatureValue
                                      * temperatureValue * temperatureValue );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "tree = " << treeLevelPotential( fieldConfiguration )
    << std::endl
    << "polynomial corrections = "
    << polynomialLoopCorrections( fieldConfiguration )
    << std::endl
    << "1-loop = " << ( treeLevelPotential( fieldConfiguration )
        + polynomialLoopCorrections( fieldConfiguration )
        + ScalarBosonCorrections( fieldConfiguration,
                                  inverseTemperatureSquared,
                                  temperatureFourthed )
        + WeylFermionCorrections( fieldConfiguration,
                                  inverseTemperatureSquared,
                                  temperatureFourthed )
        + GaugeBosonCorrections( fieldConfiguration,
                                 inverseTemperatureSquared,
                                 temperatureFourthed ) );
    std::cout << std::endl;*/


    return ( treeLevelPotential( fieldConfiguration )
             + polynomialLoopCorrections( fieldConfiguration )
             + ScalarBosonCorrections( fieldConfiguration,
                                       inverseTemperatureSquared,
                                       temperatureFourthed )
             + WeylFermionCorrections( fieldConfiguration,
                                       inverseTemperatureSquared,
                                       temperatureFourthed )
             + GaugeBosonCorrections( fieldConfiguration,
                                      inverseTemperatureSquared,
                                      temperatureFourthed ) );
  }

  // This puts all index brackets into a consistent form.
  inline std::string PotentialFromPolynomialAndMasses::FormatVariable(
                                 std::string const& unformattedVariable ) const
  {
    return RunningParameterManager::FormatVariable( unformattedVariable );
  }

  // This interprets stringToParse as a sum of real polynomial terms and sets
  // polynomialSum accordingly.
  inline void PotentialFromPolynomialAndMasses::ParseSumOfPolynomialTerms(
                                              std::string const& stringToParse,
                                                 PolynomialSum& polynomialSum )
  {
    std::pair< PolynomialSum, PolynomialSum > complexSum;
    ParseSumOfPolynomialTerms( stringToParse,
                               complexSum );
    if( !(complexSum.second.PolynomialTerms().empty()) )
    {
      throw std::runtime_error(
                      "Polynomial that should be real has imaginary factor!" );
    }
    polynomialSum.PolynomialTerms() = complexSum.first.PolynomialTerms();
  }

} /* namespace VevaciousPlusPlus */
#endif /* POTENTIALFROMPOLYNOMIALANDMASSES_HPP_ */
