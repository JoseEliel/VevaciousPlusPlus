/*
 * MinimizingPotentialOnHypersurfaces.cpp
 *
 *  Created on: Sep 2, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "BounceActionEvaluation/BounceActionPathFinding/MinimizingPotentialOnHypersurfaces.hpp"

namespace VevaciousPlusPlus
{

  MinimizingPotentialOnHypersurfaces::MinimizingPotentialOnHypersurfaces(
                                    PotentialFunction const& potentialFunction,
                                              size_t const movesPerImprovement,
                                             unsigned int const minuitStrategy,
                                       double const minuitToleranceFraction ) :
    MinuitPathFinder( movesPerImprovement,
                      minuitStrategy,
                      minuitToleranceFraction ),
    potentialFunction( potentialFunction ),
    numberOfFields( potentialFunction.NumberOfFieldVariables() ),
    pathNodes(),
    reflectionMatrix( numberOfFields,
                      numberOfFields )
  {
    // This constructor is just an initialization list.
  }

  MinimizingPotentialOnHypersurfaces::~MinimizingPotentialOnHypersurfaces()
  {
    // This does nothing.
  }


  // This resets the BouncePathFinder so that it sets up currentPath as its
  // initial path between the given vacua. It also resets pathCanBeImproved
  // and sets pathTemperature appropriately.
  TunnelPath const* MinimizingPotentialOnHypersurfaces::SetInitialPath(
                                           PotentialMinimum const& falseVacuum,
                                            PotentialMinimum const& trueVacuum,
                                                TunnelPath const* startingPath,
                                                 double const pathTemperature )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MinimizingPotentialOnHypersurfaces::ImprovePath()";
    std::cout << std::endl;
    return NULL;/**/

    this->pathTemperature = pathTemperature;
  }


  // This sets up reflectionMatrix to be the Householder reflection matrix
  // which reflects the axis of field 0 to be parallel to the difference
  // between startNode and endNode.
  void MinimizingPotentialOnHypersurfaces::SetUpHouseholderReflection( std::vector< double > const& startNode,
                                   std::vector< double > const& endNode );

  // This takes the numberOfFields-1-dimensional vector and prepends a 0 to
  // make an numberOfFields-dimensional Eigen::VectorXd.
  Eigen::VectorXd MinimizingPotentialOnHypersurfaces::UntransformedNode(
                   std::vector< double > const& nodeParameterization ) const;


  // This sets up pathFactory, pathTemperature, currentMinuitTolerance, and
  // currentMinuitResult to be appropriate for the initial path.
  void MinimizingPotentialOnHypersurfaces::SetUpPathFactoryAndMinuit( PotentialMinimum const& falseVacuum,
                                  PotentialMinimum const& trueVacuum,
                                  TunnelPath const* startingPath = NULL,
                                  double const pathTemperature = 0.0 );
} /* namespace VevaciousPlusPlus */
