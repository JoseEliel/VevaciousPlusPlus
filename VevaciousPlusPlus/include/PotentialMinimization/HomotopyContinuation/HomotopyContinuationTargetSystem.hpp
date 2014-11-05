/*
 * HomotopyContinuationTargetSystem.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef HOMOTOPYCONTINUATIONTARGETSYSTEM_HPP_
#define HOMOTOPYCONTINUATIONTARGETSYSTEM_HPP_

#include "CommonIncludes.hpp"
#include "SlhaManagement/SlhaUpdatePropagator.hpp"

namespace VevaciousPlusPlus
{
  class HomotopyContinuationTargetSystem : public SlhaUpdatePropagator
  {
  public:
    HomotopyContinuationTargetSystem(
                                    SlhaUpdatePropagator& previousPropagator );
    HomotopyContinuationTargetSystem( SlhaManager& slhaManager );
    virtual ~HomotopyContinuationTargetSystem();


    // This should evaluate the target system and place the values in
    // destinationVector.
    virtual void HomotopyContinuationSystemValues(
                   std::vector< std::complex< double > > solutionConfiguration,
          std::vector< std::complex< double > >& destinationVector ) const = 0;

    // This should be the real-valued version of the above.
    virtual void HomotopyContinuationSystemValues(
                                   std::vector< double > solutionConfiguration,
                          std::vector< double >& destinationVector ) const = 0;

    // This should evaluate the derivatives of the target system and place the
    // values in destinationMatrix.
    virtual void HomotopyContinuationSystemGradients(
                   std::vector< std::complex< double > > solutionConfiguration,
        std::vector< std::vector< std::complex< double > > >& destinationMatrix
                                                                   ) const = 0;

    // This should return a vector of field values corresponding to the field
    // configuration as it should be passed to operator() for evaluating the
    // potential, given a vector of values that solve this instance's homotopy
    // continuation system. It should return an empty vector if
    // homotopyContinuatioConfiguration does not correspond to a valid field
    // configuration. (For example, RgeImprovedOneLoopPotential uses the
    // logarithm of the renormalization scale as an extra variable in the
    // homotopy continuation, so homotopyContinuatioConfiguration actually has
    // an extra entry compared to a valid field configuration. Also, it may
    // be that homotopyContinuatioConfiguration did not correspond to a valid
    // solution where the scale is close to the Euclidean length of the field
    // configuration, so it would not correspond to a valid solution.)
    virtual std::vector< double > ValidFieldsFromHomotopyContinuation(
                std::vector< double > homotopyContinuationConfiguration ) const
    { return homotopyContinuationConfiguration; }

    // This appends all sign-flip combinations of solutionConfiguration to
    // realSolutions if the combination is both purely real and solves the
    // homotopy continuation target system (checked by considering whether
    // each entry of the gradient changes sign when its variable is changed by
    // +- resolutionSize). It does not append any solutions that are too close
    // to those already in realSolutions, by being within a hypercube of side
    // resolutionSize centered on the existing entry.
    void AppendPureRealSolutionAndValidSignFlips(
            std::vector< std::complex< double > > const& solutionConfiguration,
                           std::vector< std::vector< double > >& realSolutions,
                                                 double const resolutionSize );

  protected:
    // This allows derived classes to veto particular solutions from being
    // appended to realSolutions in AppendPureRealSolutionAndValidSignFlips.
    // The numerical tolerance of whether the system is solved is given by
    // equationTolerance.
    virtual bool AllowedSolution(
                std::vector< double > const& solutionConfiguration,
                              double const equationTolerance = 1.0 ) const = 0;
  };

} /* namespace VevaciousPlusPlus */
#endif /* HOMOTOPYCONTINUATIONTARGETSYSTEM_HPP_ */
