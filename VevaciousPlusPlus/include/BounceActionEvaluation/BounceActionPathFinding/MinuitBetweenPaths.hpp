/*
 * MinuitBetweenPaths.hpp
 *
 *  Created on: Sep 2, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITBETWEENPATHS_HPP_
#define MINUITBETWEENPATHS_HPP_

#include "CommonIncludes.hpp"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include "PotentialMinimization/GradientBasedMinimization/MinuitMinimum.hpp"
#include "../BounceActionCalculator.hpp"
#include "../PathParameterization/TunnelPath.hpp"
#include "../PathParameterization/LinearSplineThroughNodes.hpp"

namespace VevaciousPlusPlus
{

  class MinuitBetweenPaths : public ROOT::Minuit2::FCNBase
  {
  public:
    MinuitBetweenPaths( BounceActionCalculator* const bounceActionCalculator,
                        unsigned int const minuitStrategy,
                        double const minuitToleranceFraction,
                        size_t const movesPerImprovement );
    virtual ~MinuitBetweenPaths();


    // This takes minuitParameters as a pair of weightings and creates a node
    // path of weighted averages of each node in curvedPath with the
    // corresponding node in straightPath, where the weighting for the
    // straightPath node is
    // minuitParameters[ 0 ] + ( minuitParameters[ 1 ] * the fraction along the
    // path of the node), and the weighting for the curvedPath node is
    // 1.0 - the other weighting, then it returns the bounce action along the
    // path from the composed nodes.
    virtual double
    operator()( std::vector< double > const& minuitParameters ) const;

    // This implements Up() for ROOT::Minuit2::FCNBase just to stick to a basic
    // value.
    virtual double Up() const { return 1.0; }

    // This calls Minuit2's MIGRAD function on this MinuitBetweenPaths object
    // (as a ROOT::Minuit2::FCNBase) and then returns the tunneling path
    // corresponding to the minimum found after movesPerImprovement function
    // calls.
    TunnelPath const* ImprovePath();

    // This sets up straightPath to go from curvedPath.front() to
    // curvedPath.back() in a straight line with as many nodes as curvedPath
    // has, as well as updating pathTemperature and currentMinuitTolerance.
    void UpdateNodes( std::vector< std::vector< double > > const& curvedPath,
                      double const pathTemperature );


  protected:
    std::vector< std::vector< double > > const* curvedPath;
    double pathTemperature;
    BounceActionCalculator* const bounceActionCalculator;
    std::vector< std::vector< double > > straightPath;
    size_t numberOfFields;
    size_t numberOfSegments;
    unsigned int const minuitStrategy;
    size_t const movesPerImprovement;
    double const minuitToleranceFraction;
    double currentMinuitTolerance;
    MinuitMinimum currentMinuitResult;


    // This should create a set of nodes based on minuitParameters, curvedPath,
    // and straightPath, and return the path of straight lines between the
    // composed nodes.
    virtual TunnelPath const* PathForParameterization(
                     std::vector< double > const& minuitParameters ) const = 0;

    // This should prepare currentMinuitResult based on the updated curvedPath.
    virtual void PrepareMinuitStartingPoint() = 0;
  };




  // This takes the tunneling path from minuitParameters according to
  // PathForParameterization( minuitParameters ) and returns the bounce action
  // along that path.
  inline double MinuitBetweenPaths::operator()(
                          std::vector< double > const& minuitParameters ) const
  {
    TunnelPath const*
    composedPath( PathForParameterization( minuitParameters ) );
    double const bounceAction( (*bounceActionCalculator)( *composedPath ) );
    delete composedPath;
    return bounceAction;
  }

  // This calls Minuit2's MIGRAD function on this MinuitBetweenPaths object
  // (as a ROOT::Minuit2::FCNBase) and then returns the tunneling path
  // corresponding to the minimum found.
  inline TunnelPath const* MinuitBetweenPaths::ImprovePath()
  {
    ROOT::Minuit2::MnMigrad mnMigrad( (*this),
                                      currentMinuitResult.VariableValues(),
                                      currentMinuitResult.VariableErrors(),
                                      minuitStrategy );
    MinuitMinimum minuitMinimum( currentMinuitResult.VariableValues().size(),
                                 mnMigrad( movesPerImprovement,
                                           currentMinuitTolerance ) );
    currentMinuitResult = minuitMinimum;
    return PathForParameterization( minuitMinimum.VariableValues() );
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINUITBETWEENPATHS_HPP_ */
