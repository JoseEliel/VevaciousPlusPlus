/*
 * PotentialFromPolynomialAndMasses.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POTENTIALFROMPOLYNOMIALANDMASSES_HPP_
#define POTENTIALFROMPOLYNOMIALANDMASSES_HPP_

#include "../../CommonIncludes.hpp"
#include "PotentialFunction.hpp"
#include "../RunningParameterManager.hpp"
#include "../PolynomialSum.hpp"
#include "../MassesSquaredCalculators.hpp"
#include "../ThermalFunctions.hpp"
#include "../../PotentialMinimization/HomotopyContinuation.hpp"
#include "IWritesPythonPotential.hpp"

namespace VevaciousPlusPlus
{
  class PotentialFromPolynomialAndMasses : public PotentialFunction,
                                           public IWritesPythonPotential
  {
  public:
    PotentialFromPolynomialAndMasses( std::string const& modelFilename,
                              RunningParameterManager& runningParameterManager,
                                 double const scaleRangeMinimumFactor = 10.0 );
    virtual
    ~PotentialFromPolynomialAndMasses();


    // This writes the potential as
    // def PotentialFunction( fv ): return ...
    // in pythonFilename for fv being an array of floating-point numbers in the
    // same order as they are for the field configurations as internal to this
    // C++ code. It uses the purely virtual function SetScaleInPythonFunction.
    virtual void WriteAsPython( std::string const pythonFilename ) const;

    // This gives out a PolynomialGradientTargetSystem reference while setting
    // needToUpdateHomotopyContinuation to true.
    PolynomialGradientTargetSystem& HomotopyContinuationTargetSystem();

    // This is for debugging.
    std::string AsDebuggingString() const;


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
    double squareOfMaximumRenormalizationScale;
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
    bool needToUpdateHomotopyContinuation;
    double const scaleRangeMinimumFactor;


    // This is just for derived classes.
    PotentialFromPolynomialAndMasses(
                            RunningParameterManager& runningParameterManager );

    // This is just for derived classes.
    PotentialFromPolynomialAndMasses(
                          PotentialFromPolynomialAndMasses const& copySource );

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

    // This should provide the lines of Python to set the scale according to
    // the field values given by fv in the code generated by WriteAsPython
    // above, in the PotentialFunction function.
    virtual std::string SetScaleInPythonFunction() const = 0;

    // This should give out a PolynomialGradientTargetSystem reference
    // appropriate to this potential.
    virtual PolynomialGradientTargetSystem&
    GetHomotopyContinuationTargetSystem() = 0;

    // This should return a string that is valid Python indented by 4 spaces to
    // set the scale Q, given by Q^(-2) called invQSq, for an evaluation of the
    // potential assuming that the field configuration is given as an array
    // called "fv" and the temperature by T, given by T^(-2) called invTSq, given
    // the rest of the Python code written by WriteAsPython. By default Q is left
    // as the lowest scale given by the blocks of the SLHA file.
    virtual std::string SetScaleForPythonPotentialCall() const
    { return std::string( "" ); }

    // This scales down the values of cappedFieldConfiguration if the sum of
    // the squares of the elements is larger than
    // squareOfMaximumRenormalizationScale so that the sum is equal to it, and
    // returns the difference of the original sum of squares from
    // squareOfMaximumRenormalizationScale.
    double CapFieldConfiguration(
                       std::vector< double >& cappedFieldConfiguration ) const;
  };




  // This gives out a PolynomialGradientTargetSystem reference while setting
  // needToUpdateHomotopyContinuation to true.
  inline PolynomialGradientTargetSystem&
  PotentialFromPolynomialAndMasses::HomotopyContinuationTargetSystem()
  {
    needToUpdateHomotopyContinuation = true;
    return GetHomotopyContinuationTargetSystem();
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

  // This scales down the values of cappedFieldConfiguration if the sum of
  // the squares of the elements is larger than
  // squareOfMaximumRenormalizationScale so that the sum is equal to it, and
  // returns the difference of the original sum of squares from
  // squareOfMaximumRenormalizationScale.
  inline double PotentialFromPolynomialAndMasses::CapFieldConfiguration(
                        std::vector< double >& cappedFieldConfiguration ) const
  {
    double lengthSquared( 0.0 );
    for( std::vector< double >::const_iterator
         fieldValue( cappedFieldConfiguration.begin() );
         fieldValue < cappedFieldConfiguration.end();
         ++fieldValue )
    {
      lengthSquared += ( (*fieldValue) * (*fieldValue) );
    }
    if( lengthSquared <= squareOfMaximumRenormalizationScale )
    {
      return 0.0;
    }

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "cappedFieldConfiguration = "
    << FieldConfigurationAsMathematica( cappedFieldConfiguration )
    << std::endl << "lengthSquared = " << lengthSquared
    << ", squareOfMaximumRenormalizationScale = "
    << squareOfMaximumRenormalizationScale;
    std::cout << std::endl;*/

    double const
    scaleFactor( sqrt( squareOfMaximumRenormalizationScale / lengthSquared ) );
    for( unsigned int fieldIndex( 0 );
         fieldIndex < cappedFieldConfiguration.size();
         ++fieldIndex )
    {
      cappedFieldConfiguration[ fieldIndex ] *= scaleFactor;
    }

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "cappedFieldConfiguration -> "
    << FieldConfigurationAsMathematica( cappedFieldConfiguration )
    << std::endl << "returning "
    << ( lengthSquared - squareOfMaximumRenormalizationScale );
    std::cout << std::endl;*/

    return ( lengthSquared - squareOfMaximumRenormalizationScale );
  }

} /* namespace VevaciousPlusPlus */
#endif /* POTENTIALFROMPOLYNOMIALANDMASSES_HPP_ */
