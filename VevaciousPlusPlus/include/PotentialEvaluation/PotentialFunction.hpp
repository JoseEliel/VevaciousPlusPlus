/*
 * PotentialFunction.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POTENTIALFUNCTION_HPP_
#define POTENTIALFUNCTION_HPP_

#include "../LagrangianParameterManagement/ParameterUpdatePropagator.hpp"
#include "CommonIncludes.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include "PotentialMinimization/HomotopyContinuation/PolynomialGradientTargetSystem.hpp"

namespace VevaciousPlusPlus
{
  class PotentialFunction : public ParameterUpdatePropagator
  {
  public:
    PotentialFunction( ParameterUpdatePropagator& parameterUpdatePropagator );
    PotentialFunction( PotentialFunction const& copySource );
    virtual ~PotentialFunction();


    size_t NumberOfFieldVariables() const
    { return numberOfFields; }

    std::vector< std::string > const& FieldNames() const{ return fieldNames; }

    std::string FieldConfigurationAsMathematica(
                       std::vector< double > const& fieldConfiguration ) const;

    // This should return the energy density in GeV^4 of the potential for a
    // state strongly peaked around expectation values (in GeV) for the fields
    // given by the values of fieldConfiguration and temperature in GeV given
    // by temperatureValue.
    virtual double
    operator()( std::vector< double > const& fieldConfiguration,
                double const temperatureValue = 0.0 ) const = 0;

    // This returns operator(), but could be over-ridden so that a partial
    // result from a derived class could be returned, such as say the
    // tree-level part of a potential that builds on it with loop
    // corrections.
    virtual double
    QuickApproximation( std::vector< double > const& fieldConfiguration,
                        double const temperatureValue = 0.0 )
    { return (*this)( fieldConfiguration,
                      temperatureValue ); }

    // This numerically evaluates the gradient at fieldConfiguration based on
    // steps of numericalStepSize GeV in each field direction and places the
    // gradient in gradientVector. Derived classes which can analytically
    // evaluate the gradient can over-write this function.
    virtual void SetAsGradientAt( std::vector< double >& gradientVector,
                               std::vector< double > const& fieldConfiguration,
                                  double const numericalStepSize = 1.0 ) const;

    // This should return the square of the scale (in GeV^2) relevant to
    // tunneling between the given minima for this potential.
    virtual double
    ScaleSquaredRelevantToTunneling( PotentialMinimum const& falseVacuum,
                                PotentialMinimum const& trueVacuum ) const = 0;

    // This is for ease of getting the index of a field of a given name. It
    // returns the largest possible unsigned int (-1 should tick over to that)
    // if fieldName was not found. Hence calling code can check that the return
    // from this function is less than NumberOfFieldVariables().
    size_t FieldIndex( std::string const& fieldName ) const;

    std::vector< double > const& DsbFieldValues() const
    { return dsbFieldValueInputs; }

    std::vector< double > FieldValuesOrigin() const
    { return std::vector< double >( numberOfFields,
                                    0.0 ); }


    // There are a few functions that return objects which the derived class
    // might not have but might be asked for by code anyway. It's not ideal to
    // have code dedicated to "this might or might not be able to do that", but
    // it's a lot better than dynamic casting at runtime.

    // This should give out a pointer to the PolynomialGradientTargetSystem
    // appropriate to this potential, if it has one.
    virtual PolynomialGradientTargetSystem* HomotopyContinuationTargetSystem()
    { return NULL; }


  protected:
    ParameterUpdatePropagator& parameterUpdatePropagator;
    std::vector< std::string > fieldNames;
    size_t numberOfFields;
    std::vector< double > dsbFieldValueInputs;
  };




  inline std::string PotentialFunction::FieldConfigurationAsMathematica(
                        std::vector< double > const& fieldConfiguration ) const
  {
    std::stringstream stringBuilder;
    stringBuilder << "{";
    for( size_t fieldIndex( 0 );
         fieldIndex < fieldConfiguration.size();
         ++fieldIndex )
    {
      if( fieldIndex != 0 )
      {
        stringBuilder << ",";
      }
      stringBuilder << " " << fieldNames[ fieldIndex ] << " -> "
      << fieldConfiguration[ fieldIndex ];
    }
    stringBuilder << " }";
    return stringBuilder.str();
  }

  // This numerically evaluates the gradient at fieldConfiguration based on
  // steps of numericalStepSize GeV in each field direction, places the
  // the gradient in gradientVector, and returns the potential evaluated at
  // fieldConfiguration. Derived classes which can analytically evaluate the
  // gradient can over-write this function, as well as SetAsGradientAt( ... ).
  inline void
  PotentialFunction::SetAsGradientAt( std::vector< double >& gradientVector,
                               std::vector< double > const& fieldConfiguration,
                                         double const numericalStepSize ) const
  {
    gradientVector.resize( numberOfFields );
    double const potentialValue( (*this)( fieldConfiguration ) );
    std::vector< double > displacementVector( fieldConfiguration );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      displacementVector[ fieldIndex ] += numericalStepSize;
      gradientVector[ fieldIndex ]
      = ( ( (*this)( displacementVector ) - potentialValue )
          / numericalStepSize );
      displacementVector[ fieldIndex ] -= numericalStepSize;
    }
  }

  // This is for ease of getting the index of a field of a given name. It
  // returns the largest possible unsigned int (-1 should tick over to that) if
  // fieldName was not found. Hence calling code can check that the return from
  // this function is less than NumberOfFieldVariables().
  inline size_t
  PotentialFunction::FieldIndex( std::string const& fieldName ) const
  {
    for( size_t fieldIndex( 0 );
         fieldIndex < fieldNames.size();
         ++fieldIndex )
    {
      if( fieldName.compare( fieldNames[ fieldIndex ] ) == 0 )
      {
        return fieldIndex;
      }
    }
    return -1;
  }

} /* namespace VevaciousPlusPlus */
#endif /* POTENTIALFUNCTION_HPP_ */
