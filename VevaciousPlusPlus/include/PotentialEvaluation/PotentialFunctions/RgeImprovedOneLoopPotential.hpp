/*
 * RgeImprovedOneLoopPotential.hpp
 *
 *  Created on: Nov 5, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef RGEIMPROVEDONELOOPPOTENTIAL_HPP_
#define RGEIMPROVEDONELOOPPOTENTIAL_HPP_

#include "CommonIncludes.hpp"
#include "PotentialFromPolynomialWithMasses.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include "LagrangianParameterManagement/LagrangianParameterManager.hpp"

namespace VevaciousPlusPlus
{

  class RgeImprovedOneLoopPotential : public PotentialFromPolynomialWithMasses,
                                      public BOL::BasicObserver
  {
  public:
    RgeImprovedOneLoopPotential( std::string const& modelFilename,
                               double const assumedPositiveOrNegativeTolerance,
                      LagrangianParameterManager& lagrangianParameterManager );
    RgeImprovedOneLoopPotential(
                    PotentialFromPolynomialWithMasses const& potentialToCopy );
    virtual ~RgeImprovedOneLoopPotential();


    // This returns the energy density in GeV^4 of the potential for a state
    // strongly peaked around expectation values (in GeV) for the fields given
    // by the values of fieldConfiguration and temperature in GeV given by
    // temperatureValue.
    virtual double
    operator()( std::vector< double > const& fieldConfiguration,
                double const temperatureValue = 0.0 ) const;

    // This returns the square of the Euclidean distance between the given
    // vacua in field space.
    virtual double
    ScaleSquaredRelevantToTunneling( PotentialMinimum const& falseVacuum,
                                     PotentialMinimum const& trueVacuum ) const
    { return trueVacuum.SquareDistanceTo( falseVacuum ); }

    // This updates the minimum scale to use when evaluating Lagrangian
    // parameters (as just using the Euclidean length of the field
    // configuration would lead to taking the logarithm of 0 when evaluating
    // the potential at the field origin).
    virtual void respondToObservedSignal();

    // This returns a string that is valid Python with no indentation to
    // evaluate the potential in three functions:
    // TreeLevelPotential( fv ), JustLoopCorrectedPotential( fv ), and
    // LoopAndThermallyCorrectedPotential( fv ).
    virtual std::string WriteActualPythonFunction() const;


  protected:
    double minimumScaleSquared;
    double maximumScaleSquared;
  };





  // This updates the minimum scale to use when evaluating Lagrangian
  // parameters (as just using the Euclidean length of the field
  // configuration would lead to taking the logarithm of 0 when evaluating
  // the potential at the field origin).
  inline void RgeImprovedOneLoopPotential::respondToObservedSignal()
  {
    double const
    minimumScale( lagrangianParameterManager.MinimumEvaluationScale() );
    minimumScaleSquared = ( minimumScale * minimumScale );
    double const
    maximumScale( lagrangianParameterManager.MaximumEvaluationScale() );
    maximumScaleSquared = ( maximumScale * maximumScale );
    if( minimumScale > maximumScale )
    {
      std::stringstream errorBuilder;
      errorBuilder
      << "Somehow minimum allowed scale (" << minimumScale
      << " GeV) is greater than the maximum allowed scale (" << maximumScale
      << ") for this parameter point.";
      throw std::runtime_error( errorBuilder.str() );
    }
  }

} /* namespace VevaciousPlusPlus */

#endif /* RGEIMPROVEDONELOOPPOTENTIAL_HPP_ */
