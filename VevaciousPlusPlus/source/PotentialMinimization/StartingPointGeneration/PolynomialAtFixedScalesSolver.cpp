/*
 * PolynomialAtFixedScalesSolver.cpp
 *
 *  Created on: Nov 20, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialMinimization/StartingPointGeneration/PolynomialAtFixedScalesSolver.hpp"

namespace VevaciousPlusPlus
{

  PolynomialAtFixedScalesSolver::PolynomialAtFixedScalesSolver(
                    ParametersAndFieldsProductSum const& polynomialToExtremize,
                  LagrangianParameterManager const& lagrangianParameterManager,
                          std::unique_ptr<PolynomialSystemSolver> polynomialSystemSolver,
                                             unsigned int const numberOfScales,
                                         bool const returnOnlyPolynomialMinima,
                                                size_t const numberOfFields ) :
    StartingPointFinder(),
    minimizationConditions( numberOfFields,
                            std::vector< ParametersAndFieldsProductTerm >() ),
    polynomialHermitian(),
    lagrangianParameterManager( lagrangianParameterManager ),
    polynomialSystemSolver( std::move(polynomialSystemSolver) ),
    numberOfScales( numberOfScales ),
    returnOnlyPolynomialMinima( returnOnlyPolynomialMinima ),
    numberOfFields( numberOfFields )
  {
    if( returnOnlyPolynomialMinima )
    {
      polynomialHermitian
      = std::vector< std::vector< ParametersAndFieldsProductSum > >(
                                                                numberOfFields,
                  std::vector< ParametersAndFieldsProductSum >( numberOfFields,
                                           ParametersAndFieldsProductSum() ) );
    }

    // The minimization conditions are the set of partial derivatives of the
    // polynomial with respect to the fields.
    for( std::vector< ParametersAndFieldsProductTerm >::const_iterator
         polynomialTerm(
                 polynomialToExtremize.ParametersAndFieldsProducts().begin() );
         polynomialTerm
         != polynomialToExtremize.ParametersAndFieldsProducts().end();
         ++polynomialTerm )
    {
      for( size_t fieldIndex( 0 );
           fieldIndex < numberOfFields;
           ++fieldIndex )
      {
        if( polynomialTerm->NonZeroDerivative( fieldIndex ) )
        {
          minimizationConditions[ fieldIndex ].push_back(
                             polynomialTerm->PartialDerivative( fieldIndex ) );

          if( returnOnlyPolynomialMinima )
          {
            // We take advantage of the Hessian being real and symmetric to
            // avoid bothering to calculate half of the off-diagonal second
            // derivatives.
            ParametersAndFieldsProductTerm const&
            firstDerivative( minimizationConditions[ fieldIndex ].back() );
            std::vector< ParametersAndFieldsProductSum >&
            fieldRow( polynomialHermitian[ fieldIndex ] );
            for( size_t secondFieldIndex( 0 );
                 secondFieldIndex <= fieldIndex;
                 ++secondFieldIndex )
            {
              if( firstDerivative.NonZeroDerivative( secondFieldIndex ) )
              {
                fieldRow[ secondFieldIndex ].ParametersAndFieldsProducts(
                                                                   ).push_back(
                                             firstDerivative.PartialDerivative(
                                                          secondFieldIndex ) );
              }
            }
          }
        }
      }
    }
  }

  PolynomialAtFixedScalesSolver::~PolynomialAtFixedScalesSolver()
  {
  }


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
  void PolynomialAtFixedScalesSolver::operator()(
                   std::vector< std::vector< double > >& startingPoints ) const
  {
    if( numberOfScales == 0 )
    {
      throw std::runtime_error(
                "PolynomialAtFixedScalesSolver::numberOfScales cannot be 0!" );
    }
    else if( numberOfScales == 1 )
    {
      AddSolutions( startingPoints,
               log( lagrangianParameterManager.AppropriateSingleFixedScale() ),
                    0.0,
                    -1.0 );
    }
    else
    {
      double const logLowestScale(
                  log( lagrangianParameterManager.MinimumEvaluationScale() ) );
      double const logHighestScale(
                  log( lagrangianParameterManager.MaximumEvaluationScale() ) );
      double const
      logStep( ( logHighestScale - logLowestScale ) / numberOfScales );

      // Take solutions at the lowest scale, with the allowed solution length
      // range being [ 0.0, exp( logLowestScale + logStep )].
      AddSolutions( startingPoints,
                    logLowestScale,
                    0.0,
                    exp( logLowestScale + logStep ) );

      // Take solutions from every intermediate scale, with the allowed range
      // being exp(+/- logStep) around each scale.
      double logCurrentScale( logLowestScale );
      for( unsigned int scaleStep( 1 );
           scaleStep < ( numberOfScales - 1 );
           ++scaleStep )
      {
        logCurrentScale += logStep;
        AddSolutions( startingPoints,
                      logCurrentScale,
                      exp( logCurrentScale - logStep ),
                      exp( logCurrentScale + logStep ) );
      }

      // Take solutions at the highest scale, with the allowed solution length
      // lower bound being exp(logCurrentScale) which should be
      // exp( logHighestScale - logStep ), and there is no upper bound.
      AddSolutions( startingPoints,
                    logHighestScale,
                    exp( logCurrentScale ),
                    -1.0 );
    }
  }

  // This uses polynomialSystemSolver to solve the system at the scale given
  // by exp(logCurrentScale), discarding solutions with Euclidean length
  // smaller than lowerSolutionLengthBound or greater than
  // upperSolutionLengthBound. (If upperSolutionLengthBound is <= 0.0, then the
  // upper limit is not applied.)
  void PolynomialAtFixedScalesSolver::AddSolutions(
                          std::vector< std::vector< double > >& startingPoints,
                                                  double const logCurrentScale,
                                         double const lowerSolutionLengthBound,
                                  double const upperSolutionLengthBound ) const
  {
    std::vector< double > lagrangianParameters;
    lagrangianParameterManager.ParameterValues( logCurrentScale,
                                                lagrangianParameters );
    std::vector< std::vector< double > > solutionSet;
    PolynomialSystemSolver::ConstraintSystem
    polynomialConstraints( numberOfFields );
    PolynomialConstraints( lagrangianParameters,
                           logCurrentScale,
                           polynomialConstraints );

    // Now polynomialSystemSolver does its job.
    (*polynomialSystemSolver)( polynomialConstraints,
                               solutionSet );

    std::vector< std::vector< double > > solutionsInRange;
    for( std::vector< std::vector< double > >::const_iterator
         solutionVector( solutionSet.begin() );
         solutionVector != solutionSet.end();
         ++solutionVector )
    {
      double const solutionLengthSquared( VectorUtilities::LengthSquared(
                                                           *solutionVector ) );
      if( ( solutionLengthSquared
            >= ( lowerSolutionLengthBound * lowerSolutionLengthBound ) )
          &&
          ( ( upperSolutionLengthBound <= 0.0 )
            ||
            ( solutionLengthSquared
              <= ( upperSolutionLengthBound * upperSolutionLengthBound ) ) ) )
      {
        solutionsInRange.push_back( *solutionVector );
      }
    }
    if( returnOnlyPolynomialMinima )
    {
      AddOnlyPolynomialMinima( lagrangianParameters,
                               solutionsInRange,
                               startingPoints );
    }
    else
    {
      startingPoints.insert( startingPoints.end(),
                             solutionsInRange.begin(),
                             solutionsInRange.end() );
    }
  }

  // This puts the polynomial constraint system for the given scale into
  // polynomialConstraints.
  void PolynomialAtFixedScalesSolver::PolynomialConstraints(
                             std::vector< double > const& lagrangianParameters,
                                                  double const logCurrentScale,
        PolynomialSystemSolver::ConstraintSystem& polynomialConstraints ) const
  {
    double polynomialCoefficient( 0.0 );
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      for( std::vector< ParametersAndFieldsProductTerm >::const_iterator
           polynomialTerm( minimizationConditions[ fieldIndex ].begin() );
           polynomialTerm != minimizationConditions[ fieldIndex ].end();
           ++polynomialTerm )
      {
        polynomialCoefficient
        = polynomialTerm->CoefficientFactor( lagrangianParameters );
        if( polynomialCoefficient != 0.0 )
        {
          polynomialConstraints[ fieldIndex ].push_back(
                                      PolynomialSystemSolver::FactorWithPowers(
                                                         polynomialCoefficient,
                                      polynomialTerm->FieldPowersByIndex() ) );
        }
      }
      if( polynomialConstraints[ fieldIndex ].empty() )
      {
        std::stringstream errorBuilder;
        errorBuilder
        << "At scale " << exp( logCurrentScale )
        << " the minimization condition for field " << fieldIndex
        << " had no non-zero terms.";
        throw std::runtime_error( errorBuilder.str() );
      }
    }
  }

  // This appends the solutions from candidateSolutions which have positive
  // semi-definite eigenvalues for their Hessian matrices.
  void PolynomialAtFixedScalesSolver::AddOnlyPolynomialMinima(
                             std::vector< double > const& lagrangianParameters,
                std::vector< std::vector< double > > const& candidateSolutions,
                  std::vector< std::vector< double > >& startingPoints ) const
  {
    Eigen::MatrixXd valuesMatrix( numberOfFields,
                                  numberOfFields );

    for( std::vector< std::vector< double > >::const_iterator
         candidateSolution( candidateSolutions.begin() );
         candidateSolution != candidateSolutions.end();
         ++candidateSolution )
    {
      // Each solution needs its Hessian constructed out of the solution with
      // the Lagrangian parameter values.
      for( size_t firstFieldIndex( 0 );
           firstFieldIndex < numberOfFields;
           ++firstFieldIndex )
      {
        // We take advantage of the Hessian being real and symmetric to
        // avoid bothering to calculate half of the off-diagonal second
        // derivatives.
        valuesMatrix( firstFieldIndex,
                      firstFieldIndex )
        = polynomialHermitian[ firstFieldIndex ][ firstFieldIndex ](
                                                          lagrangianParameters,
                                                          *candidateSolution );
        for( size_t secondFieldIndex( 0 );
             secondFieldIndex < firstFieldIndex;
             ++secondFieldIndex )
        {
          valuesMatrix( firstFieldIndex,
                        secondFieldIndex ) = valuesMatrix( secondFieldIndex,
                                                           firstFieldIndex )
          = polynomialHermitian[ firstFieldIndex ][ secondFieldIndex ](
                                                          lagrangianParameters,
                                                          *candidateSolution );
        }
      }

      Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd >
      eigenvalueFinder( valuesMatrix,
                        Eigen::EigenvaluesOnly );

      // If any of the eigenvalues are negative, the solution is not appended
      // to startingPoints.
      bool noNegativeEigenvalues( true );
      for( size_t eigenvalueIndex( 0 );
           eigenvalueIndex < numberOfFields;
           ++eigenvalueIndex )
      {
        if( eigenvalueFinder.eigenvalues()( eigenvalueIndex ) < 0.0 )
        {
          noNegativeEigenvalues = false;
        }
      }
      if( noNegativeEigenvalues )
      {
        startingPoints.push_back( *candidateSolution );
      }
    }
  }

} /* namespace VevaciousPlusPlus */
