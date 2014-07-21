/*
 * LinearSplineThroughNodesFactory.hpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef LINEARSPLINETHROUGHNODESFACTORY_HPP_
#define LINEARSPLINETHROUGHNODESFACTORY_HPP_

#include "CommonIncludes.hpp"
#include "PathFromNodesFactory.hpp"
#include "LinearSplineThroughNodes.hpp"
#include "NodesFromParameterization.hpp"

namespace VevaciousPlusPlus
{

  class LinearSplineThroughNodesFactory : public PathFromNodesFactory
  {
  public:
    LinearSplineThroughNodesFactory( size_t const numberOfFields,
                  NodesFromParameterization* const nodesFromParameterization );
    virtual ~LinearSplineThroughNodesFactory();


    // This returns a pointer to a new LinearSplineThroughNodes.
    // The calling code is responsible for memory management! (It'd be great to
    // return a std::unique_Ptr< TunnelPath >, but we're stubbornly sticking to
    // the requirements not including a C++11-compliant compiler.)
    virtual TunnelPath*
    operator()( std::vector< std::vector< double > > const& pathNodes,
                double const pathTemperature = 0.0 ) const
    { return new LinearSplineThroughNodes( pathNodes,
                                           pathTemperature ); }
  };

} /* namespace VevaciousPlusPlus */
#endif /* LINEARSPLINETHROUGHNODESFACTORY_HPP_ */
