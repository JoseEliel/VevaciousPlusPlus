/*
 * Hom4ps2Runner.hpp
 *
 *  Created on: Apr 14, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef HOM4PS2RUNNER_HPP_
#define HOM4PS2RUNNER_HPP_

#include "../StandardIncludes.hpp"
#include "HomotopyContinuationSolver.hpp"
#include "../PotentialEvaluation/PotentialFunctions/PotentialFunctions.hpp"


namespace VevaciousPlusPlus
{

  class Hom4ps2Runner : public HomotopyContinuationSolver
  {
  public:
    Hom4ps2Runner( HomotopyContinuationReadyPolynomial& polynomialPotential,
                   std::string const& pathToHom4ps2,
                   std::string const homotopyType = "1" );
    virtual
    ~Hom4ps2Runner();


    // This uses HOM4PS2 to fill purelyRealSolutionSets with all the extrema of
    // polynomialPotential.TargetPolynomialGradient().
    virtual void FindTreeLevelExtrema(
                std::vector< std::vector< double > >& purelyRealSolutionSets );


  protected:
    HomotopyContinuationReadyPolynomial& polynomialPotential;
    std::string const pathToHom4ps2;
    std::string const homotopyType;
    BOL::StringParser variableNamer;
    std::vector< std::complex< long double > > complexSolutions;
    std::vector< std::string > variableNames;
    std::map< std::string, unsigned int > nameToIndexMap;

    void WriteHom4p2Input( std::string const& hom4ps2InputFilename );

    void ParseHom4ps2Output( std::string const& hom4ps2OutputFilename,
                std::vector< std::vector< double > >& purelyRealSolutionSets );
  };

} /* namespace VevaciousPlusPlus */
#endif /* HOM4PS2RUNNER_HPP_ */
