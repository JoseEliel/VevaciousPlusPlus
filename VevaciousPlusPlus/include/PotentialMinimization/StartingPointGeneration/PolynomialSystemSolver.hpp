/*
 * PolynomialSystemSolver.hpp
 *
 *  Created on: Nov 20, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALSYSTEMSOLVER_HPP_
#define POLYNOMIALSYSTEMSOLVER_HPP_

#include <utility>
#include <vector>
#include "Utilities/VectorUtilities.hpp"

namespace VevaciousPlusPlus
{

  class PolynomialSystemSolver
  {
  public:
    typedef std::pair< double, std::vector< std::vector::size_type > >
            FactorWithPowers;
    typedef std::vector< FactorWithPowers > PolynomialConstraint;
    typedef std::vector< PolynomialConstraint > ConstraintSystem;

    PolynomialSystemSolver();
    virtual ~PolynomialSystemSolver();


    // This should solve the system of polynomial constraints and put all
    // the solutions into systemSolutions.
    virtual void operator()(
                      std::vector< PolynomialConstraint > const& systemToSolve,
             std::vector< std::vector< double > >& systemSolutions ) const = 0;


  protected:
    // This is useful if the main method somehow manages to miss some vaild
    // solutions where some of the fields have flipped sign, as we observed
    // would happen sometimes with HOM4PS2. It finds all the possible solutions
    // which have one or more of solutionConfiguration's elements flipped in
    // sign, and appends solutionConfiguration and all of the sign-flip
    // variations which pass IsValidSolution(...) to the end of solutionSet.
    static void AppendSolutionAndValidSignFlips(
                            std::vector< double > const& solutionConfiguration,
                             std::vector< std::vector< double > >& solutionSet,
                      std::vector< PolynomialConstraint > const& systemToSolve,
                                                 double const resolutionSize );

    // This puts solutionToFlip and every possible variant where at least one
    // element has had a sign flip into signFlips.
    static void AllSignFlips( std::vector< double > solutionToFlip,
                             std::vector< std::vector< double > >& signFlips );

    // This returns true if givenSolution is within a hypercube of side
    // resolutionSize with any of the solutions in solutionSet.
    static bool
    AlreadyInSolutionSet( std::vector< double > const& givenSolution,
                          std::vector< std::vector< double > >& solutionSet,
                          double const resolutionSize );

    // This appends any solution in candidateSolutions which is not already in
    // destinationSolutions and which also solves systemToSolve, within
    // resolutionSize.
    static void AppendValidSolutions(
                std::vector< std::vector< double > > const& candidateSolutions,
                    std::vector< std::vector< double > >& destinationSolutions,
                      std::vector< PolynomialConstraint > const& systemToSolve,
                                      double const resolutionSize );

    // This goes through each of the constraints in systemToSolve, stepping
    // resolutionSize either side of givenSolution in the field appropriate to
    // the constraint, and returns false if any of the constraints do not
    // change sign in stepping from one side of givenSolution to the other in
    // any of the fields.
    static bool IsValidSolution( std::vector< double > givenSolution,
                      std::vector< PolynomialConstraint > const& systemToSolve,
                                 double const resolutionSize );

    // This returns the value of fieldConstraint for the field values given in
    // fieldConfiguration.
    static double PartialSlope( PolynomialConstraint const& fieldConstraint,
                                std::vector< double > fieldConfiguration );
  };





  // This is useful if the main method somehow manages to miss some vaild
  // solutions where some of the fields have flipped sign, as we observed
  // would happen sometimes with HOM4PS2. It finds all the possible solutions
  // which have one or more of solutionConfiguration's elements flipped in
  // sign, and appends solutionConfiguration and all of the sign-flip
  // variations which pass IsValidSolution(...) to the end of solutionSet.
  inline void PolynomialSystemSolver::AppendSolutionAndValidSignFlips(
                            std::vector< double > const& solutionConfiguration,
                             std::vector< std::vector< double > >& solutionSet,
                      std::vector< PolynomialConstraint > const& systemToSolve,
                                                  double const resolutionSize )
  {
    if( !(AlreadyInSolutionSet( solutionConfiguration,
                                solutionSet,
                                resolutionSize )) )
    {
      solutionSet.push_back( solutionConfiguration );
      std::vector< std::vector< double > > signFlips;
      AllSignFlips( solutionConfiguration,
                    signFlips );
      AppendValidSolutions( signFlips,
                            solutionSet,
                            systemToSolve,
                            resolutionSize );
    }
  }

  // This puts solutionToFlip and every possible variant where at least one
  // element has had a sign flip into signFlips.
  inline void
  PolynomialSystemSolver::AllSignFlips( std::vector< double > solutionToFlip,
                              std::vector< std::vector< double > >& signFlips )
  {
    signFlips.push_back( solutionToFlip );
    std::vector::size_type const numberOfElements( solutionToFlip.size() );
    for( std::vector::size_type flipIndex( 0 );
         flipIndex < numberOfElements;
         ++flipIndex )
    {
      std::vector::size_type const vectorSize( signFlips.size() );
      for( std::vector::size_type vectorIndex( 0 );
           vectorIndex < vectorSize;
           ++vectorIndex )
      {
        if( signFlips[ vectorIndex ][ flipIndex ] != 0.0 )
        {
          solutionToFlip = signFlips[ vectorIndex ];
          solutionToFlip[ flipIndex ] = -(solutionToFlip[ flipIndex ]);
          signFlips.push_back( solutionToFlip );
        }
      }
    }
  }

  // This appends any solution in candidateSolutions which is not already in
  // destinationSolutions and which also solves systemToSolve, within
  // resolutionSize.
  inline void PolynomialSystemSolver::AppendValidSolutions(
                std::vector< std::vector< double > > const& candidateSolutions,
                    std::vector< std::vector< double > >& destinationSolutions,
                      std::vector< PolynomialConstraint > const& systemToSolve,
                                                  double const resolutionSize )
  {
    for( std::vector< std::vector< double > >::const_iterator
         candidateSolution( candidateSolutions.begin() );
         candidateSolution != candidateSolutions.end();
         ++candidateSolution )
    {
      if( !(AlreadyInSolutionSet( *candidateSolution,
                                  destinationSolutions,
                                  resolutionSize ))
          &&
          IsValidSolution( *candidateSolution,
                           systemToSolve,
                           resolutionSize ) )
      {
        destinationSolutions.push_back( *candidateSolution );
      }
    }
  }

  // This returns true if givenSolution is within a hypercube of side
  // resolutionSize with any of the solutions in solutionSet.
  inline bool PolynomialSystemSolver::AlreadyInSolutionSet(
                                    std::vector< double > const& givenSolution,
                             std::vector< std::vector< double > >& solutionSet,
                                                  double const resolutionSize )
  {
    for( std::vector< std::vector< double > >::const_iterator
         existingSolution( solutionSet.begin() );
         existingSolution < solutionSet.end();
         ++existingSolution )
    {
      if( VectorUtilities::DifferenceIsWithinHypercube( givenSolution,
                                                        *existingSolution,
                                                        resolutionSize ) )
      {
        return true;
      }
    }
    return false;
  }

} /* namespace VevaciousPlusPlus */

#endif /* POLYNOMIALSYSTEMSOLVER_HPP_ */
