/*
 * PolynomialThroughNodesFactory.hpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALTHROUGHNODESFACTORY_HPP_
#define POLYNOMIALTHROUGHNODESFACTORY_HPP_

#include "CommonIncludes.hpp"
#include "PathFromNodesFactory.hpp"
#include "PolynomialThroughNodes.hpp"
#include "NodesFromParameterization.hpp"
#include "Eigen/Dense"

namespace VevaciousPlusPlus
{

  class PolynomialThroughNodesFactory : public PathFromNodesFactory
  {
  public:
    PolynomialThroughNodesFactory( size_t const numberOfFields,
                  NodesFromParameterization* const nodesFromParameterization );
    virtual ~PolynomialThroughNodesFactory();


    // This returns a pointer to a new PolynomialThroughNodes.
    // The calling code is responsible for memory management! (It'd be great to
    // return a std::unique_Ptr< TunnelPath >, but we're stubbornly sticking to
    // the requirements not including a C++11-compliant compiler.)
    virtual TunnelPath*
    operator()( std::vector< std::vector< double > > const& pathNodes,
                std::vector< double > const& pathParameterization,
                double const pathTemperature = 0.0 ) const
    { return new PolynomialThroughNodes( pathNodes,
                                         pathParameterization,
                                         pathTemperature,
                                         pathStepsInverse ); }


  protected:
    Eigen::MatrixXd pathStepsInverse;
  };

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALTHROUGHNODESFACTORY_HPP_ */
