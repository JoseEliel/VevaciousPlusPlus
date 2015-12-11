/*
 * PolynomialAtFixedScalesSolver.hpp
 *
 *  Created on: Nov 20, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALATFIXEDSCALESSOLVER_HPP_
#define POLYNOMIALATFIXEDSCALESSOLVER_HPP_

#include "PotentialMinimization/StartingPointFinder.hpp"
#include "PotentialEvaluation/BuildingBlocks/ParametersAndFieldsProductSum.hpp"
#include "LagrangianParameterManagement/LagrangianParameterManager.hpp"
#include "PolynomialSystemSolver.hpp"
#include <vector>
#include "PotentialEvaluation/BuildingBlocks/ParametersAndFieldsProductTerm.hpp"
#include <stdexcept>
#include <cmath>
#include "Utilities/VectorUtilities.hpp"
#include <sstream>
#include "Eigen/Dense"


namespace VevaciousPlusPlus
{

  class PolynomialAtFixedScalesSolver : public StartingPointFinder
  {
  public:
    PolynomialAtFixedScalesSolver(
                    ParametersAndFieldsProductSum const& polynomialToExtremize,
                  LagrangianParameterManager const& lagrangianParameterManager,
                          PolynomialSystemSolver* const polynomialSystemSolver,
                                   unsigned int const numberOfScales,
                                   bool const returnOnlyPolynomialMinima,
                                   unsigned int const numberOfFields );
    virtual ~PolynomialAtFixedScalesSolver();


    // This uses polynomialSystemSolver repeatedly to find all the solutions of
    // the polynomial minimization conditions at a set of fixed scales, and
    // puts them into startingPoints. The scale used and whether solutions are
    // used or discarded depend on numberOfScales.
    // If numberOfScales is 1, the scale used is
    // lagrangianParameterManager.AppropriateSingleFixedScale().
    // If numberOfScales is 2 or more, the lowest scale used is
    // lagrangianParameterManager.MinimumEvaluationScale()
    // and the highest scale used is
    // lagrangianParameterManager.MaximumEvaluationScale().
    // Intermediate scales form a geometric sequence between these scales.
    // For example, if numberOfScales is 5, the lowest scale is 10^3, and the
    // highest is 10^11, the scales used would be 10^3, 10^5, 10^7, 10^9, and
    // 10^11.
    // For each scale, a set of solutions is found, and from that set, any
    // solutions with Euclidean length lower than the scale below the current
    // scale or larger than the scale above the current scale are discarded.
    // Using the above example, from the solution set found at 10^5, any
    // solutions with length less than 10^3 or greater than 10^7 are discarded.
    // The lowest scale does not discard any solution for being too small, and
    // the highest scale does not discard any solutions for being too large.
    virtual void
    operator()( std::vector< std::vector< double > >& startingPoints ) const;


  protected:
    std::vector< std::vector< ParametersAndFieldsProductTerm > >
    minimizationConditions;
    std::vector< std::vector< ParametersAndFieldsProductSum > >
    polynomialHermitian;
    LagrangianParameterManager const& lagrangianParameterManager;
    PolynomialSystemSolver* polynomialSystemSolver;
    unsigned int const numberOfScales;
    bool const returnOnlyPolynomialMinima;
    unsigned int const numberOfFields;

    // This uses polynomialSystemSolver to solve the system at the scale given
    // by exp(logCurrentScale), discarding solutions with Euclidean length
    // smaller than lowerSolutionLengthBound or greater than
    // upperSolutionLengthBound. (If upperSolutionLengthBound is <= 0.0, then
    // the upper limit is not applied.)
    void AddSolutions( std::vector< std::vector< double > >& startingPoints,
                       double const logCurrentScale,
                       double const lowerSolutionLengthBound,
                       double const upperSolutionLengthBound ) const;

    // This puts the polynomial constraint system for the given scale into
    // polynomialConstraints.
    void
    PolynomialConstraints( std::vector< double > const& lagrangianParameters,
                           double const logCurrentScale,
       PolynomialSystemSolver::ConstraintSystem& polynomialConstraints ) const;

    // This appends the solutions from candidateSolutions which have positive
    // semi-definite eigenvalues for their Hessian matrices.
    void AddOnlyPolynomialMinima(
                             std::vector< double > const& lagrangianParameters,
                std::vector< std::vector< double > > const& candidateSolutions,
                  std::vector< std::vector< double > >& startingPoints ) const;
  };

} /* namespace VevaciousPlusPlus */

#endif /* POLYNOMIALATFIXEDSCALESSOLVER_HPP_ */
