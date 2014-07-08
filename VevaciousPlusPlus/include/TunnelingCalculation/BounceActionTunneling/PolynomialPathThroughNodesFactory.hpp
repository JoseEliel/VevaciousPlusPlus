/*
 * PolynomialPathThroughNodesFactory.hpp
 *
 *  Created on: Jul 8, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POLYNOMIALPATHTHROUGHNODESFACTORY_HPP_
#define POLYNOMIALPATHTHROUGHNODESFACTORY_HPP_

#include "CommonIncludes.hpp"
#include "PathFromNodesFactory.hpp"
#include "PolynomialPathThroughNodes.hpp"

namespace VevaciousPlusPlus
{

  class PolynomialPathThroughNodesFactory : public PathFromNodesFactory
  {
  public:
    PolynomialPathThroughNodesFactory(
                                      std::vector< double > const& falseVacuum,
                                       std::vector< double > const& trueVacuum,
                                       std::string const& xmlArguments );
    virtual ~PolynomialPathThroughNodesFactory();


    // This returns a pointer to a new PolynomialPathThroughNodes.
    // The calling code is responsible for memory management! (It'd be great to
    // return a std::unique_Ptr< TunnelPath >, but we're stubbornly sticking to
    // the requirements not including a C++11-compliant compiler.)
    virtual TunnelPath*
    operator()( std::vector< std::vector< double > > const& pathNodes,
                double const pathTemperature = 0.0 ) const
    { return new PolynomialPathThroughNodes( pathNodes,
                                             pathTemperature,
                                             pathStepsInverse ); }


  protected:
    Eigen::MatrixXd pathStepsInverse;
  };

} /* namespace VevaciousPlusPlus */
#endif /* POLYNOMIALPATHTHROUGHNODESFACTORY_HPP_ */
