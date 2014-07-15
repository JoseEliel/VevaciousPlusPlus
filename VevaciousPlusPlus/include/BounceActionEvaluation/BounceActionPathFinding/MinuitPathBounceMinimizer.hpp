/*
 * MinuitPathBounceMinimizer.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITPATHBOUNCEMINIMIZER_HPP_
#define MINUITPATHBOUNCEMINIMIZER_HPP_

#include "CommonIncludes.hpp"
#include "FullPathVaryingMinuit.hpp"
#include "../PathParameterization/TunnelPathFactory.hpp"
#include "../BounceActionCalculator.hpp"

namespace VevaciousPlusPlus
{

  class MinuitPathBounceMinimizer : public FullPathVaryingMinuit
  {
  public:
    MinuitPathBounceMinimizer( TunnelPathFactory* pathFactory,
                               BounceActionCalculator* bounceActionCalculator,
                               std::string const& xmlArguments );
    virtual ~MinuitPathBounceMinimizer();


    // It may seem unwise to have this object call Minuit on itself, but really
    // it's just a handy way of keeping the minimization function within the
    // class that ends up finding its minimum. In this case,
    // nodeParameterization is just the parameterization of the node.
    virtual double
    operator()( std::vector< double > const& pathParameterization ) const;


  protected:
    BounceActionCalculator* bounceActionCalculator;
  };




  // It may seem unwise to have this object call Minuit on itself, but really
  // it's just a handy way of keeping the minimization function within the
  // class that ends up finding its minimum. The node index is set by
  // ImprovePath before the Minuit minimization and the temperature is set by
  // SetInitialPath.
  inline double MinuitPathBounceMinimizer::operator()(
                      std::vector< double > const& pathParameterization ) const
  {
    return (*bounceActionCalculator)( (*pathFactory)( pathParameterization,
                                                      pathTemperature ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* MINUITPATHBOUNCEMINIMIZER_HPP_ */
