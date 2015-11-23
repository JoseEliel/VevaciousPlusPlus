/*
 * OldHom4ps2Runner.hpp
 *
 *  Created on: Apr 14, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef OLDHOM4PS2RUNNER_HPP_
#define OLDHOM4PS2RUNNER_HPP_

#include "CommonIncludes.hpp"
#include "HomotopyContinuationSolver.hpp"
#include "PolynomialGradientTargetSystem.hpp"


namespace VevaciousPlusPlus
{

  class OldHom4ps2Runner : public HomotopyContinuationSolver
  {
  public:
    OldHom4ps2Runner( PolynomialGradientTargetSystem& targetSystem,
                      std::string const& pathToHom4ps2,
                      std::string const& homotopyType );
    virtual ~OldHom4ps2Runner();


    // This uses HOM4PS2 to fill startingPoints with all the extrema of
    // targetSystem.TargetPolynomialGradient().
    virtual void
    operator()( std::vector< std::vector< double > >& startingPoints );

    // This is just a dummy to make things work while I have it side-by-side to
    // the new version for comparison.
    virtual void operator()(
                   std::vector< std::vector< double > >& startingPoints ) const
    { ; /* This does nothing. */ }

  protected:
    PolynomialGradientTargetSystem& targetSystem;
    std::string const pathToHom4ps2;
    std::string const homotopyType;
    BOL::StringParser variableNamer;
    std::vector< std::complex< long double > > complexSolutions;
    std::vector< std::string > variableNames;
    std::map< std::string, size_t > nameToIndexMap;

    void WriteHom4p2Input( std::string const& hom4ps2InputFilename );

    void ParseHom4ps2Output( std::string const& hom4ps2OutputFilename,
                std::vector< std::vector< double > >& purelyRealSolutionSets );
  };

} /* namespace VevaciousPlusPlus */
#endif /* OLDHOM4PS2RUNNER_HPP_ */
