/*
 * PHCRunner.hpp
 *
 *  Created on: Nov 22, 2017
 *      Author: Simon Geisler (simon.geisler94@gmail.com)
 */

#ifndef PHCRUNNER_HPP_
#define PHCRUNNER_HPP_

#include "PotentialMinimization/StartingPointGeneration/PolynomialSystemSolver.hpp"
#include <string>
#include <vector>
#include <map>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <climits>
#include <unistd.h>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <streambuf>
#include "LHPC/Utilities/ParsingUtilities.hpp"
#include <complex>
#include <cmath>
#include <regex>
#include <sys/stat.h>
#include <chrono>
namespace VevaciousPlusPlus
{

  class PHCRunner : public PolynomialSystemSolver
  {
  public:
    PHCRunner( std::string const& pathToPHC,
                   double const resolutionSize );
    virtual ~PHCRunner();


    // This uses PHC to fill systemSolutions with all the solutions of
    // systemToSolve.
    virtual void
    operator()( std::vector< PolynomialConstraint > const& systemToSolve,
                std::vector< std::vector< double > >& systemSolutions ) const;


  protected:
    static std::string const fieldNamePrefix;

    std::string const pathToPHC;
    double const resolutionSize;

    // This sets up the variable names in variableNames and nameToIndexMap,
    // then writes systemToSolve using these names in the correct form for
    // PHC in a file with name PHCInputFilename.
    void
    WritePHCInput( std::vector< PolynomialConstraint > const& systemToSolve,
                      std::vector< std::string >& variableNames,
                      std::map< std::string, size_t >& nameToIndexMap,
                      std::string const& PHCInputFilename ) const;

    // This returns the constraint as a string of terms joined by '+' or '-'
    // appropriately, where each term is of the form
    // coefficient " * " variableName[ fieldIndex ] "^" appropriate power
    // (without writing any power part if the power is only 1, and without
    // writing the field name at all if its power is 0).
    std::string WritePHCConstraint( PolynomialConstraint const& constraintToWrite,
                       std::vector< std::string > const& variableNames ) const;

	// This takes the solutioncontainers, which are appended to the Inputfile of PHC
	// and fills the purely real Solutionvectors.
    void ParsePHCOutput( std::string const& PHCInputFilename,
                  std::vector< std::vector< double > >& purelyRealSolutionSets,
                             std::vector< std::string > const& variableNames,
                         std::map< std::string, size_t > const& nameToIndexMap,
              std::vector< PolynomialConstraint > const& systemToSolve ) const;
  };

} /* namespace VevaciousPlusPlus */

#endif /* PHCRUNNER_HPP_ */
