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
#include "PotentialMinimization/PotentialMinimum.hpp"

namespace VevaciousPlusPlus
{

  class PathFromNodesFactory : public TunnelPathFactory
  {
  public:
    PathFromNodesFactory( size_t const numberOfFields,
                        NodesFromParameterization* nodesFromParameterization );
    virtual ~PathFromNodesFactory();



    // This resets nodesFromParameterization so that it will produce
    // TunnelPath*s that parameterize the path between the given vacua.
    virtual void SetVacua( PotentialMinimum const& falseVacuum,
                           PotentialMinimum const& trueVacuum );

    // This gets the nodes from nodesFromParameterization and then calls
    // operator()( std::vector< std::vector< double > > const& pathNodes ).
    virtual TunnelPath*
    operator()( std::vector< double > const& pathParameterization,
                double const pathTemperature = 0.0 ) const;


    // This should return a pointer to an appropriate derived class instance.
    // The calling code is responsible for memory management! (It'd be great to
    // return a std::unique_Ptr< TunnelPath >, but we're stubbornly sticking to
    // the requirements not including a C++11-compliant compiler.)
    virtual TunnelPath*
    operator()( std::vector< std::vector< double > > const& pathNodes,
                std::vector< double > const& pathParameterization,
                double const pathTemperature = 0.0 ) const = 0;

    NodesFromParameterization& GetNodesFromParameterization()
    { return *nodesFromParameterization; }


  protected:
    NodesFromParameterization* nodesFromParameterization;
  };





  // This resets nodesFromParameterization so that it will produce
  // TunnelPath*s that parameterize the path between the given vacua.
  inline void
  PathFromNodesFactory::SetVacua( PotentialMinimum const& falseVacuum,
                                  PotentialMinimum const& trueVacuum )
  {
    nodesFromParameterization->SetVacua( falseVacuum,
                                         trueVacuum );
    initialStepSizes = nodesFromParameterization->InitialStepSizes();
  }

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
                    pathParameterization,
                    pathTemperature );
  }

} /* namespace VevaciousPlusPlus */
#endif /* PATHFROMNODESFACTORY_HPP_ */
