/*
 * FieldPolynomialsWithScale.hpp
 *
 *  Created on: Apr 16, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef FIELDPOLYNOMIALSWITHSCALE_HPP_
#define FIELDPOLYNOMIALSWITHSCALE_HPP_

#include "../../LagrangianParameterManagement/ParameterUpdatePropagator.hpp"
#include "CommonIncludes.hpp"
#include "PolynomialGradientTargetSystem.hpp"
#include "BasicFunctions/PolynomialTerm.hpp"
#include "BasicFunctions/PolynomialSum.hpp"

namespace VevaciousPlusPlus
{

  class FieldPolynomialsWithScale : public PolynomialGradientTargetSystem
  {
  public:
    FieldPolynomialsWithScale( PolynomialSum const& potentialPolynomial,
                               size_t const numberOfVariables,
                               ParameterUpdatePropagator& previousPropagator,
                            std::vector< size_t > const& fieldsAssumedPositive,
                            std::vector< size_t > const& fieldsAssumedNegative,
                               bool const treeLevelMinimaOnlyAsValidSolutions,
                             double const assumedPositiveOrNegativeTolerance );
    virtual ~FieldPolynomialsWithScale();


    // This returns the first numberOfFields entries in
    // homotopyContinuationConfiguration as a vector (i.e. without the last
    // entry) if the scale variable (which should approximate the logarithm of
    // the scale), which is the last entry of
    // homotopyContinuationConfiguration, exponentiates to a value between
    // (0.5 * minimumScale) and (2.0 * maximumScale). Otherwise, an empty
    // vector is returned. (The values of minimumScale and maximumScale are set
    // by the last call of PrepareForHomotopyContinuation).
    virtual std::vector< double > ValidFieldsFromHomotopyContinuation(
               std::vector< double > homotopyContinuationConfiguration ) const;

    // This calls the base PolynomialGradientTargetSystem version then replaces
    // targetSystem[ numberOfFields ] (which should be the scale constraint)
    // with a polynomial approximation of the constraint l = ln( Q ) where l is
    // the variable with index numberOfFields and Q is the scale taken as the
    // square root of the sum of the squares of the field values. Hence it
    // inserts a polynomial approximation of
    // exp( 2 * l ) = minimumScale^2 + f[0]^2 + f[1]^2 + ...
    // for the field configuration f[]. The minimum renormalization scale
    // should be given by lowerEndOfStartValues, and upperEndOfStartValues
    // should give the maximum renormalization scale.
    void UpdateSelfForNewSlha( SlhaManager const& slhaManager );
  };




  // This returns the first numberOfFields entries in
  // homotopyContinuationConfiguration as a vector (i.e. without the last
  // entry) if the scale variable (which should approximate the logarithm of
  // the scale), which is the last entry of
  // homotopyContinuationConfiguration, exponentiates to a value between
  // (0.5 * minimumScale) and (2.0 * maximumScale). Otherwise, an empty vector
  // is returned. (The values of minimumScale and maximumScale are set by the
  // last call of PrepareForHomotopyContinuation).
  inline std::vector< double >
  FieldPolynomialsWithScale::ValidFieldsFromHomotopyContinuation(
                std::vector< double > homotopyContinuationConfiguration ) const
  {
    std::vector< double > returnVector;
    double const
    renormalizationScale( exp( homotopyContinuationConfiguration.back() ) );
    if( ( renormalizationScale > ( 0.5 * minimumScale ) )
        &&
        ( renormalizationScale < ( 2.0 * maximumScale ) ) )
    {
      returnVector.assign( homotopyContinuationConfiguration.begin(),
                           ( homotopyContinuationConfiguration.end() - 1 ) );
    }
    return returnVector;
  }

} /* namespace VevaciousPlusPlus */
#endif /* FIELDPOLYNOMIALSWITHSCALE_HPP_ */
