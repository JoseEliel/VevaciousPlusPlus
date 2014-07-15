/*
 * MinuitPathPotentialMinimizer.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef MINUITPATHPOTENTIALMINIMIZER_HPP_
#define MINUITPATHPOTENTIALMINIMIZER_HPP_

#include "CommonIncludes.hpp"
#include "FullPathVaryingMinuit.hpp"
#include "../PathParameterization/TunnelPathFactory.hpp"
#include "PotentialEvaluation/PotentialFunction.hpp"

namespace VevaciousPlusPlus
{

  class MinuitPathPotentialMinimizer : public FullPathVaryingMinuit
  {
  public:
    MinuitPathPotentialMinimizer( TunnelPathFactory* pathFactory,
                                  PotentialFunction const& potentialFunction,
                                  std::string const& xmlArguments );
    virtual ~MinuitPathPotentialMinimizer();


    // It may seem unwise to have this object call Minuit on itself, but really
    // it's just a handy way of keeping the minimization function within the
    // class that ends up finding its minimum. In this case,
    // nodeParameterization is just the parameterization of the node.
    virtual double
    operator()( std::vector< double > const& pathParameterization ) const;


  protected:
    PotentialFunction const& potentialFunction;
    size_t numberOfPotentialSamplePoints;
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
#endif /* MINUITPATHPOTENTIALMINIMIZER_HPP_ */
