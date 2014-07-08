/*
 * QuadraticSplinePathThroughNodesFactory.hpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef QUADRATICSPLINEPATHTHROUGHNODESFACTORY_HPP_
#define QUADRATICSPLINEPATHTHROUGHNODESFACTORY_HPP_

#include "CommonIncludes.hpp"
#include "PathFromNodesFactory.hpp"
#include "QuadraticSplinePathThroughNodes.hpp"

namespace VevaciousPlusPlus
{

  class QuadraticSplinePathThroughNodesFactory
  {
  public:
    QuadraticSplinePathThroughNodesFactory(
                                      std::vector< double > const& falseVacuum,
                                       std::vector< double > const& trueVacuum,
                                            std::string const& xmlArguments );
    virtual ~QuadraticSplinePathThroughNodesFactory();


    // This returns a pointer to a new QuadraticSplinePathThroughNodes.
    // The calling code is responsible for memory management! (It'd be great to
    // return a std::unique_Ptr< TunnelPath >, but we're stubbornly sticking to
    // the requirements not including a C++11-compliant compiler.)
    virtual TunnelPath*
    operator()( std::vector< std::vector< double > > const& pathNodes,
                double const pathTemperature = 0.0 ) const
    { return new QuadraticSplinePathThroughNodes( pathNodes,
                                                  pathTemperature ); }
  };

} /* namespace VevaciousPlusPlus */
#endif /* QUADRATICSPLINEPATHTHROUGHNODESFACTORY_HPP_ */
