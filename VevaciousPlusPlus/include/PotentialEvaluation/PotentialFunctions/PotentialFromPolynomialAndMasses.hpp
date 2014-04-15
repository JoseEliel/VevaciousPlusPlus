/*
 * PotentialFromPolynomialAndMasses.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POTENTIALFROMPOLYNOMIALANDMASSES_HPP_
#define POTENTIALFROMPOLYNOMIALANDMASSES_HPP_

#include "../../StandardIncludes.hpp"
#include "HomotopyContinuationReadyPolynomial.hpp"
#include "../PolynomialTerm.hpp"
#include "../PolynomialSum.hpp"
#include "../MassesSquaredCalculators/MassesSquaredCalculators.hpp"
#include "../RunningParameterManager.hpp"
#include "../ThermalFunctions.hpp"

namespace VevaciousPlusPlus
{
  class PotentialFromPolynomialAndMasses :
                                     public HomotopyContinuationReadyPolynomial
  {
  public:
    PotentialFromPolynomialAndMasses( std::string const& modelFilename,
                            RunningParameterManager& runningParameterManager );
    virtual
    ~PotentialFromPolynomialAndMasses();


    // This updates all the parameters of the potential that are not field
    // values based on the values that appear in blocks in the SLHA format in
    // the file managed by runningParameters.
    virtual void UpdateParameters();


  protected:
    typedef std::pair< std::vector< double >, double > DoubleVectorWithDouble;
    static std::string const digitChars;
    static std::string const dotAndDigits;
    static std::string const allowedVariableInitials;
    static std::string const allowedVariableChars;
    static double const loopFactor;
    static double const thermalFactor;

    RunningParameterManager& runningParameters;
    std::vector< PolynomialSum > dsbFieldValuePolynomials;
    double currentMinimumRenormalizationScale;
    double squareOfMinimumRenormalizationScale;
    double currentMaximumRenormalizationScale;
    PolynomialSum treeLevelPotential;
    PolynomialSum polynomialLoopCorrections;
    std::vector< MassesSquaredCalculator* > scalarSquareMasses;
    std::vector< MassesSquaredCalculator* > fermionSquareMasses;
    std::vector< MassesSquaredCalculator* > vectorSquareMasses;
    std::vector< RealMassesSquaredMatrix > scalarMassSquaredMatrices;
    std::vector< SymmetricComplexMassMatrix > fermionMassMatrices;
    std::vector< ComplexMassSquaredMatrix > fermionMassSquaredMatrices;
    std::vector< RealMassesSquaredMatrix > vectorMassSquaredMatrices;
    double vectorMassCorrectionConstant;


    // This is just for derived classes.
    PotentialFromPolynomialAndMasses(
                            RunningParameterManager& runningParameterManager );

    // This is just for derived classes.
    PotentialFromPolynomialAndMasses(
                          PotentialFromPolynomialAndMasses const& copySource );


    // This should set dsbFieldValueInputs based on the SLHA file just read in.
    virtual void EvaluateDsbInputAndSetScale() = 0;

    // This evaluates the one-loop potential with thermal corrections assuming
    // that the scale has been set correctly.
    double
    LoopAndThermalCorrections( std::vector< double > const& fieldConfiguration,
          std::vector< DoubleVectorWithDouble > scalarMassesSquaredWithFactors,
         std::vector< DoubleVectorWithDouble > fermionMassesSquaredWithFactors,
          std::vector< DoubleVectorWithDouble > vectorMassesSquaredWithFactors,
                                              double const inverseScaleSquared,
                                         double const temperatureValue ) const;

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

    // This appends the masses-squared and multiplicity from each
    // MassesSquaredFromMatrix in massSquaredMatrices to
    // massSquaredMatrices, with all functionoids evaluated at the last scale
    // which was used to update them.
    void AddMassesSquaredWithMultiplicity(
                               std::vector< double > const& fieldConfiguration,
            std::vector< MassesSquaredCalculator* > const& massSquaredMatrices,
       std::vector< DoubleVectorWithDouble >& massesSquaredWithFactors ) const;

    // This evaluates the sum of corrections for the degrees of freedom with
    // masses-squared given by massesSquaredWithFactors with
    // subtractFromLogarithm as the constant to subtract from the logarithm of
    // the ratio of mass-squared to square of renormalization scale, at a
    // temperature given by inverseTemperatureSquared^(-1/2) using a thermal
    // correction function given by ThermalFunction, and adds them to
    // cumulativeQuantumCorrection and cumulativeThermalCorrection.
    void AddToCorrections(
         std::vector< DoubleVectorWithDouble > const& massesSquaredWithFactors,
                           double const inverseScaleSquared,
                           double const inverseTemperatureSquared,
                           double const subtractFromLogarithm,
                           double (*ThermalFunction)( double const ),
                           double& cumulativeQuantumCorrection,
                           double& cumulativeThermalCorrection ) const;

    // This sets homotopyContinuationPotentialPolynomial to be a copy of
    // treeLevelPotential.
    virtual void PrepareHomotopyContinuationPotentialPolynomial()
    { homotopyContinuationPotentialPolynomial = treeLevelPotential; }
  };




  // This updates all the parameters of the potential that are not field
  // values based on the values that appear in blocks in the SLHA format in
  // the file managed by runningParameters.
  inline void PotentialFromPolynomialAndMasses::UpdateParameters()
  {
    EvaluateDsbInputAndSetScale();
    PreparePolynomialHomotopyContinuation();
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

  // This appends the masses-squared and multiplicity from each
  // MassesSquaredFromMatrix in massSquaredMatrices to
  // massSquaredMatrices, with all functionoids evaluated at the last scale
  // which was used to update them.
  inline void
  PotentialFromPolynomialAndMasses::AddMassesSquaredWithMultiplicity(
                               std::vector< double > const& fieldConfiguration,
            std::vector< MassesSquaredCalculator* > const& massSquaredMatrices,
        std::vector< DoubleVectorWithDouble >& massesSquaredWithFactors ) const
  {
    for( std::vector< MassesSquaredCalculator* >::const_iterator
         whichMatrix( massSquaredMatrices.begin() );
         whichMatrix < massSquaredMatrices.end();
         ++whichMatrix )
    {
      massesSquaredWithFactors.push_back(
           std::make_pair( (*whichMatrix)->MassesSquared( fieldConfiguration ),
                           (*whichMatrix)->MultiplicityFactor() ) );
    }
  }

} /* namespace VevaciousPlusPlus */
#endif /* POTENTIALFROMPOLYNOMIALANDMASSES_HPP_ */
