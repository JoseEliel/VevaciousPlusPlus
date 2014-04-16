/*
 * PotentialMinimizer.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef POTENTIALMINIMIZER_HPP_
#define POTENTIALMINIMIZER_HPP_

#include "../StandardIncludes.hpp"
#include "../PotentialEvaluation.hpp"
#include "PotentialMinimum.hpp"
#include "BOLlib/include/BasicObserver.hpp"

namespace VevaciousPlusPlus
{
  class PotentialMinimizer
  {
  public:
    PotentialMinimizer( PotentialFunction& potentialFunction );
    virtual
    ~PotentialMinimizer();


    // This should find all the minima of the potential evaluated at a
    // temperature given by temperatureInGev, and record them in foundMinima.
    // It should also set dsbVacuum (which is assumed to not have been updated
    // directly by the SLHA manager yet; also it may be that the SLHA data only
    // gives the approximate position of the DSB minimum for potentialFunction
    // which may not be as high an order approximation as the calculation that
    // created the SLHA data, so the DSB input has to be taken as a starting
    // point for MINUIT, for example; furthermore, in general the DSB minimum
    // will change with temperature), and record the minima lower than
    // dsbVacuum in panicVacua, and of those, it should set panicVacuum to be
    // the prime candidate for tunneling out of dsbVacuum (by default, taken to
    // be the minimum in panicVacua closest to dsbVacuum).
    virtual void FindMinima( double const temperatureInGev = 0.0 ) = 0;

    // This should find the minimum at temperature temperatureInGev nearest to
    // minimumToAdjust (which is assumed to be a minimum of the potential at a
    // different temperature).
    virtual PotentialMinimum
    AdjustMinimumForTemperature( PotentialMinimum const& minimumToAdjust,
                                 double const temperatureInGev ) = 0;

    // This returns the set of minima of the potential which were found.
    std::vector< PotentialMinimum > const& FoundMinima() const
    { return foundMinima; }

    // This returns the minimum corresponding to the desired symmetry breaking.
    PotentialMinimum const& DsbVacuum() const{ return dsbVacuum; }

    // This returns the set of minima which were found to be lower than the DSB
    // minimum.
    std::vector< PotentialMinimum > const& PanicVacua() const
    { return panicVacua; }

    // This returns the minimum which is considered to be the prime candidate
    // for tunneling to from the DSB minimum, taken by default to be the
    // minimum from panicVacua which is closest to dsbVacuum.
    PotentialMinimum const& PanicVacuum() const { return panicVacuum; }

    bool DsbVacuumIsStable() const{ return panicVacua.empty(); }

    bool DsbVacuumIsMetaStable() const{ return !(panicVacua.empty()); }


  protected:
    PotentialFunction& potentialFunction;
    std::vector< PotentialMinimum > foundMinima;
    PotentialMinimum dsbVacuum;
    std::vector< PotentialMinimum > panicVacua;
    PotentialMinimum panicVacuum;

    // This returns false if the vector of the absolute values of the fields of
    // comparisonMinimum lies within a hypercube of side thresholdDistance
    // centered on dsbVacuum.
    bool IsNotPhaseRotationOfDsbVacuum(
                                     PotentialMinimum const& comparisonMinimum,
                                        double const thresholdDistance ) const;
  };




  // This returns false if the vector of the absolute values of the fields of
  // comparisonMinimum lies within a hypercube of side thresholdDistance
  // centered on dsbVacuum.
  inline bool PotentialMinimizer::IsNotPhaseRotationOfDsbVacuum(
                                     PotentialMinimum const& comparisonMinimum,
                                         double const thresholdDistance ) const
  {
    for( unsigned int fieldIndex( 0 );
         fieldIndex < dsbVacuum.FieldConfiguration().size();
         ++fieldIndex )
    {
      if( fabs( fabs( dsbVacuum.FieldConfiguration()[ fieldIndex ] )
               - fabs( comparisonMinimum.FieldConfiguration()[ fieldIndex ] ) )
          > thresholdDistance )
      {
        return true;
      }
    }
    return false;
  }

} /* namespace VevaciousPlusPlus */
#endif /* POTENTIALMINIMIZER_HPP_ */
