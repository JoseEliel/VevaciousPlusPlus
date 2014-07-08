/*
 * LinearSplinePathThroughNodesFactory.hpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef LINEARSPLINEPATHTHROUGHNODESFACTORY_HPP_
#define LINEARSPLINEPATHTHROUGHNODESFACTORY_HPP_

#include "CommonIncludes.hpp"
#include "PathFromNodesFactory.hpp"
#include "LinearSplinePathThroughNodes.hpp"

namespace VevaciousPlusPlus
{

  class LinearSplinePathThroughNodesFactory
  {
  public:
    LinearSplinePathThroughNodesFactory(
                                      std::vector< double > const& falseVacuum,
                                       std::vector< double > const& trueVacuum,
                                         std::string const& xmlArguments );
    virtual ~LinearSplinePathThroughNodesFactory();


    // This returns a pointer to a new LinearSplinePathThroughNodes.
    // The calling code is responsible for memory management! (It'd be great to
    // return a std::unique_Ptr< TunnelPath >, but we're stubbornly sticking to
    // the requirements not including a C++11-compliant compiler.)
    virtual TunnelPath*
    operator()( std::vector< std::vector< double > > const& pathNodes,
                double const pathTemperature = 0.0 ) const
    { return new LinearSplinePathThroughNodes( pathNodes,
                                               pathTemperature ); }
  };

} /* namespace VevaciousPlusPlus */
#endif /* LINEARSPLINEPATHTHROUGHNODESFACTORY_HPP_ */
