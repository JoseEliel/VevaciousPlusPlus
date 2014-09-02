/*
 * MinimizingPotentialOnHypersurfaces.hpp
 *
 *  Created on: Sep 2, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINIMIZINGPOTENTIALONHYPERSURFACES_HPP_
#define MINIMIZINGPOTENTIALONHYPERSURFACES_HPP_

#include "CommonIncludes.hpp"
#include "MinuitPathFinder.hpp"
#include "Minuit2/MnMigrad.h"
#include "Eigen/Dense"
#include "PotentialMinimization/GradientBasedMinimization/MinuitMinimum.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"
#include "../PathParameterization/LinearSplineThroughNodesFactory.hpp"

namespace VevaciousPlusPlus
{

  class MinimizingPotentialOnHypersurfaces : public MinuitPathFinder
  {
  public:
    MinimizingPotentialOnHypersurfaces(
                                    PotentialFunction const& potentialFunction,
                                        size_t const movesPerImprovement = 100,
                                        unsigned int const minuitStrategy = 1,
                                  double const minuitToleranceFraction = 0.5 );
    virtual ~MinimizingPotentialOnHypersurfaces();


    // This resets the BouncePathFinder so that it sets up currentPath as its
    // initial path between the given vacua. It also resets pathCanBeImproved
    // and sets pathTemperature appropriately.
    virtual TunnelPath const*
    SetInitialPath( PotentialMinimum const& falseVacuum,
                    PotentialMinimum const& trueVacuum,
                    TunnelPath const* startingPath = NULL,
                    double const pathTemperature = 0.0 );

    // This allows Minuit2 to adjust the full path a set number of times to try
    // to minimize the sum of potentials at a set of nodes or bounce action
    // along the adjusted path, and then sets the path.
    //virtual TunnelPath const* ImprovePath() = 0;

    // It may seem unwise to have this object call Minuit on itself, but really
    // it's just a handy way of keeping the minimization function within the
    // class that ends up finding its minimum. In this case,
    // nodeParameterization is just the parameterization of the node.
    //virtual double
    //operator()( std::vector< double > const& nodeParameterization ) const = 0;


  protected:
    PotentialFunction const& potentialFunction;
    size_t const numberOfFields;
    std::vector< std::vector< double > > pathNodes;
    Eigen::MatrixXd reflectionMatrix;


    // This sets up reflectionMatrix to be the Householder reflection matrix
    // which reflects the axis of field 0 to be parallel to the difference
    // between startNode and endNode.
    void SetUpHouseholderReflection( std::vector< double > const& startNode,
                                     std::vector< double > const& endNode );

    // This takes the numberOfFields-1-dimensional vector and prepends a 0 to
    // make an numberOfFields-dimensional Eigen::VectorXd.
    Eigen::VectorXd UntransformedNode(
                     std::vector< double > const& nodeParameterization ) const;


    // This sets up pathFactory, pathTemperature, currentMinuitTolerance, and
    // currentMinuitResult to be appropriate for the initial path.
    void SetUpPathFactoryAndMinuit( PotentialMinimum const& falseVacuum,
                                    PotentialMinimum const& trueVacuum,
                                    TunnelPath const* startingPath = NULL,
                                    double const pathTemperature = 0.0 );
  };

} /* namespace VevaciousPlusPlus */
#endif /* MINIMIZINGPOTENTIALONHYPERSURFACES_HPP_ */
