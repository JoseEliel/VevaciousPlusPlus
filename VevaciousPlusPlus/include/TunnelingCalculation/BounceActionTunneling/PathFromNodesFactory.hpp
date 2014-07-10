/*
 * PathFromNodesFactory.hpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef PATHFROMNODESFACTORY_HPP_
#define PATHFROMNODESFACTORY_HPP_

#include "CommonIncludes.hpp"
#include "TunnelPathFactory.hpp"
#include "NodesFromParameterization.hpp"

namespace VevaciousPlusPlus
{

  class PathFromNodesFactory : public TunnelPathFactory
  {
  public:
    PathFromNodesFactory( std::string const& xmlArguments );
    virtual ~PathFromNodesFactory();


    // This does something.
    void SetVacua( PotentialMinimum const& falseVacuum,
                   PotentialMinimum const& trueVacuum,
                   TunnelPath const* startingPath = NULL,
                   double const pathTemperature = 0.0 );

    // This gets the nodes from nodesFromParameterization and then calls
    // operator()( std::vector< std::vector< double > > const& pathNodes ).
    virtual TunnelPath*
    operator()( std::vector< double > const& pathParameterization,
                double const pathTemperature = 0.0 ) const;


    // This should return a pointer to an appropriate derived class instance.
    // The calling code is responsible for memory management! (It'd be great to
    // return a std::unique_Ptr< TunnelPath >, but we're stubbornly sticking to
    // the requirements not including a C++11-compliant compiler.)
    virtual TunnelPath* operator()(
                         std::vector< std::vector< double > > const& pathNodes,
                                double const pathTemperature = 0.0 ) const = 0;

    NodesFromParameterization& NodesFromParameterization()
    { return *nodesFromParameterization; }


  protected:
    NodesFromParameterization* nodesFromParameterization;
  };




  // This gets the nodes from nodesFromParameterization and then calls
  // operator()( std::vector< std::vector< double > > const& pathNodes ).
  inline TunnelPath* PathFromNodesFactory::operator()(
                             std::vector< double > const& pathParameterization,
                                           double const pathTemperature ) const
  {
    std::vector< std::vector< double > > pathNodes;
    nodesFromParameterization->PathNodeSet( pathNodes,
                                            pathParameterization );
    return (*this)( pathNodes,
                    pathTemperature );
  }

} /* namespace VevaciousPlusPlus */
#endif /* PATHFROMNODESFACTORY_HPP_ */
