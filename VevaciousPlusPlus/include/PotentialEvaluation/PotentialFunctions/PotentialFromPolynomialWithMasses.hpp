/*
 * PotentialFromPolynomialWithMasses.hpp
 *
 *  Created on: Oct 09, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POTENTIALFROMPOLYNOMIALWITHMASSES_HPP_
#define POTENTIALFROMPOLYNOMIALWITHMASSES_HPP_

#include "LagrangianParameterManagement/LagrangianParameterManager.hpp"
#include "CommonIncludes.hpp"
#include "boost/math/constants/constants.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "BasicFunctions/ParametersAndFieldsProductSum.hpp"
#include "PotentialEvaluation/MassesSquaredCalculators/MassesSquaredFromMatrix.hpp"
#include "PotentialEvaluation/MassesSquaredCalculators/ComplexMassSquaredMatrix.hpp"
#include "PotentialEvaluation/MassesSquaredCalculators/RealMassesSquaredMatrix.hpp"
#include "PotentialEvaluation/MassesSquaredCalculators/SymmetricComplexMassMatrix.hpp"
#include "PotentialEvaluation/MassesSquaredCalculator.hpp"
#include "PotentialEvaluation/ThermalFunctions.hpp"
#include "VersionInformation.hpp"

namespace VevaciousPlusPlus
{
  class PotentialFromPolynomialWithMasses : public PotentialFunction
  {
  public:
    PotentialFromPolynomialWithMasses( std::string const& modelFilename,
                               double const assumedPositiveOrNegativeTolerance,
                      LagrangianParameterManager& lagrangianParameterManager );
    virtual ~PotentialFromPolynomialWithMasses();


    // This returns the appropriate polynomial of Lagrangian terms which best
    // matches the full potential. This base class allows for the
    // Liebler-Porod-Staub scheme where the polynomial part of the 1-loop
    // effective potential is split into a tree-level part and a 1-loop
    // polynomial correction part (something like (m^2_tree + delta m^2) phi^2
    // for example).
    virtual ParametersAndFieldsProductSum const&
    PolynomialApproximation() const
    { return treeLevelPotential; }

    // This writes the potential as
    // def PotentialFunction( fv ): return ...
    // in pythonFilename for fv being an array of floating-point numbers in the
    // same order as they are for the field configurations as internal to this
    // C++ code. It uses the virtual function WriteActualPythonFunction.
    virtual void WriteAsPython( std::string const& pythonFilename ) const;

    // This is for debugging.
    std::string AsDebuggingString() const;


  protected:
    typedef std::pair< std::vector< double >, double > DoubleVectorWithDouble;
    typedef
    std::pair< ParametersAndFieldsProductSum, ParametersAndFieldsProductSum >
    ComplexParametersAndFieldsProductSum;
    static std::string const digitChars;
    static std::string const dotAndDigits;
    static std::string const allowedVariableInitials;
    static std::string const allowedVariableChars;
    static double const piSquared;
    static double const loopFactor;
    static double const thermalFactor;
    static std::string const positiveByConvention;
    static std::string const negativeByConvention;

    ParametersAndFieldsProductSum treeLevelPotential;
    ParametersAndFieldsProductSum polynomialLoopCorrections;
    std::vector< MassesSquaredCalculator* > scalarSquareMasses;
    std::vector< MassesSquaredCalculator* > fermionSquareMasses;
    std::vector< MassesSquaredCalculator* > vectorSquareMasses;
    std::vector< RealMassesSquaredMatrix > scalarMassSquaredMatrices;
    std::vector< SymmetricComplexMassMatrix > fermionMassMatrices;
    std::vector< ComplexMassSquaredMatrix > fermionMassSquaredMatrices;
    std::vector< RealMassesSquaredMatrix > vectorMassSquaredMatrices;
    double vectorMassCorrectionConstant;
    std::vector< size_t > fieldsAssumedPositive;
    std::vector< size_t > fieldsAssumedNegative;
    double const assumedPositiveOrNegativeTolerance;


    // This is just for derived classes.
    PotentialFromPolynomialWithMasses(
                      LagrangianParameterManager& lagrangianParameterManager );

    // This is just for derived classes.
    PotentialFromPolynomialWithMasses(
                         PotentialFromPolynomialWithMasses const& copySource );

    // This evaluates the one-loop potential with thermal corrections assuming
    // that the squared masses were evaluated at the given scale correctly.
    double
    LoopAndThermalCorrections(
          std::vector< DoubleVectorWithDouble > scalarMassesSquaredWithFactors,
         std::vector< DoubleVectorWithDouble > fermionMassesSquaredWithFactors,
          std::vector< DoubleVectorWithDouble > vectorMassesSquaredWithFactors,
                                              double const inverseScaleSquared,
                                         double const temperatureValue ) const;

    // This interprets stringToParse as a sum of complex polynomial terms and
    // sets polynomialSum accordingly.
    void ParseSumOfPolynomialTerms( std::string const& stringToParse,
                         ComplexParametersAndFieldsProductSum& polynomialSum );

    // This interprets stringToParse as a sum of real polynomial terms and sets
    // polynomialSum accordingly.
    void ParseSumOfPolynomialTerms( std::string const& stringToParse,
                                ParametersAndFieldsProductSum& polynomialSum );

    // This reads in a whole number or variable (including possible raising to
    // a power), applies the correct operation to polynomialTerm, and then
    // returns the position of the character just after the interpreted word.
    // If there was a factor of "i", "I", "j", or "J", imaginaryTerm is set to
    // true.
    size_t
    PutNextNumberOrVariableIntoPolynomial( std::string const& stringToParse,
                                           size_t wordStart,
                                    ParametersAndFieldsProduct& polynomialTerm,
                                           bool& imaginaryTerm );

    // This appends the masses-squared and multiplicity from each
    // MassesSquaredFromMatrix in massSquaredMatrices to massSquaredMatrices,
    // with the values of the Lagrangian parameters given in parameterValues.
    void AddMassesSquaredWithMultiplicity(
                                  std::vector< double > const& parameterValues,
                               std::vector< double > const& fieldConfiguration,
            std::vector< MassesSquaredCalculator* > const& massSquaredMatrices,
       std::vector< DoubleVectorWithDouble >& massesSquaredWithFactors ) const;

    // This appends the masses-squared and multiplicity from each
    // MassesSquaredFromMatrix in massSquaredMatrices to massSquaredMatrices,
    // with all Lagrangian parameters evaluated at the last scale which was
    // used to update them.
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
                           bool const temperatureGreaterThanZero,
                           double const inverseTemperatureSquared,
                           double const subtractFromLogarithm,
                           double (*ThermalFunction)( double const ),
                           double& cumulativeQuantumCorrection,
                           double& cumulativeThermalCorrection ) const;

    // This should return a string that is valid Python with no indentedation
    // to evaluate the potential in three functions:
    // TreeLevelPotential( fv ), JustLoopCorrectedPotential( fv ), and
    // LoopAndThermallyCorrectedPotential( fv ).
    virtual std::string WriteActualPythonFunction() const = 0;
  };





  // This interprets stringToParse as a sum of real polynomial terms and sets
  // polynomialSum accordingly.
  inline void PotentialFromPolynomialWithMasses::ParseSumOfPolynomialTerms(
                                              std::string const& stringToParse,
                                 ParametersAndFieldsProductSum& polynomialSum )
  {
    ComplexParametersAndFieldsProductSum complexSum;
    ParseSumOfPolynomialTerms( stringToParse,
                               complexSum );
    if( !(complexSum.second.ParametersAndFieldsProducts().empty()) )
    {
      throw std::runtime_error(
                      "Polynomial that should be real has imaginary factor!" );
    }
    polynomialSum.ParametersAndFieldsProducts()
    = complexSum.first.ParametersAndFieldsProducts();
  }

  // This appends the masses-squared and multiplicity from each
  // MassesSquaredFromMatrix in massSquaredMatrices to massSquaredMatrices,
  // with the values of the Lagrangian parameters given in parameterValues.
  inline void
  PotentialFromPolynomialWithMasses::AddMassesSquaredWithMultiplicity(
                                  std::vector< double > const& parameterValues,
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
           std::make_pair( (*whichMatrix)->MassesSquared( parameterValues,
                                                          fieldConfiguration ),
                           (*whichMatrix)->MultiplicityFactor() ) );
    }
  }

  // This appends the masses-squared and multiplicity from each
  // MassesSquaredFromMatrix in massSquaredMatrices to massSquaredMatrices,
  // with all Lagrangian parameters evaluated at the last scale which was
  // used to update them.
  inline void
  PotentialFromPolynomialWithMasses::AddMassesSquaredWithMultiplicity(
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
#endif /* POTENTIALFROMPOLYNOMIALWITHMASSES_HPP_ */
