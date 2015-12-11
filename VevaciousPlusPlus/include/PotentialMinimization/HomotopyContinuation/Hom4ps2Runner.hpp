/*
 * Hom4ps2Runner.hpp
 *
 *  Created on: Nov 23, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef HOM4PS2RUNNER_HPP_
#define HOM4PS2RUNNER_HPP_

#include "PotentialMinimization/StartingPointGeneration/PolynomialSystemSolver.hpp"
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <stdexcept>
#include <climits>
#include <unistd.h>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include "LHPC/Utilities/ParsingUtilities.hpp"
#include <cmath>

namespace VevaciousPlusPlus
{

  class Hom4ps2Runner : public PolynomialSystemSolver
  {
  public:
    Hom4ps2Runner( std::string const& pathToHom4ps2,
                   std::string const& homotopyType,
                   double const resolutionSize );
    virtual ~Hom4ps2Runner();


    // This uses HOM4PS2 to fill systemSolutions with all the solutions of
    // systemToSolve.
    virtual void operator()(
                      std::vector< PolynomialConstraint > const& systemToSolve,
                 std::vector< std::vector< double > >& systemSolutions ) const;


  protected:
    static std::string const fieldNamePrefix;

    std::string const pathToHom4ps2;
    std::string const homotopyType;
    double const resolutionSize;

    // This sets up the variable names in variableNames and nameToIndexMap,
    // then writes systemToSolve using these names in the correct form for
    // HOM4PS2 in a file with name hom4ps2InputFilename.
    void
    WriteHom4p2Input( std::vector< PolynomialConstraint > const& systemToSolve,
                      std::vector< std::string >& variableNames,
                      std::map< std::string, size_t >& nameToIndexMap,
                      std::string const& hom4ps2InputFilename ) const;

    // This returns the constraint as a string of terms joined by '+' or '-'
    // appropriately, where each term is of the form
    // coefficient " * " variableName[ fieldIndex ] "^" appropriate power
    // (without writing any power part if the power is only 1, and without
    // writing the field name at all if its power is 0).
    std::string WriteConstraint( PolynomialConstraint const& constraintToWrite,
                       std::vector< std::string > const& variableNames ) const;

    void ParseHom4ps2Output( std::string const& hom4ps2OutputFilename,
                  std::vector< std::vector< double > >& purelyRealSolutionSets,
                               std::vector< std::string > const& variableNames,
                 std::map< std::string, size_t > const& nameToIndexMap,
              std::vector< PolynomialConstraint > const& systemToSolve ) const;
  };

} /* namespace VevaciousPlusPlus */

#endif /* HOM4PS2RUNNER_HPP_ */
