/*
 * QuadraticSplineThroughNodesFactory.hpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef QUADRATICSPLINETHROUGHNODESFACTORY_HPP_
#define QUADRATICSPLINETHROUGHNODESFACTORY_HPP_

#include "CommonIncludes.hpp"
#include "PathFromNodesFactory.hpp"
#include "QuadraticSplineThroughNodes.hpp"
#include "NodesFromParameterization.hpp"

namespace VevaciousPlusPlus
{

  class QuadraticSplineThroughNodesFactory : public PathFromNodesFactory
  {
  public:
    QuadraticSplineThroughNodesFactory( size_t const numberOfFields,
                  NodesFromParameterization* const nodesFromParameterization );
    virtual ~QuadraticSplineThroughNodesFactory();


    // This returns a pointer to a new QuadraticSplineThroughNodes.
    // The calling code is responsible for memory management! (It'd be great to
    // return a std::unique_Ptr< TunnelPath >, but we're stubbornly sticking to
    // the requirements not including a C++11-compliant compiler.)
    virtual TunnelPath*
    operator()( std::vector< std::vector< double > > const& pathNodes,
                std::vector< double > const& pathParameterization,
                double const pathTemperature = 0.0 ) const
    { return new QuadraticSplineThroughNodes( pathNodes,
                                              pathParameterization,
                                              pathTemperature ); }
  };

} /* namespace VevaciousPlusPlus */
#endif /* QUADRATICSPLINETHROUGHNODESFACTORY_HPP_ */
