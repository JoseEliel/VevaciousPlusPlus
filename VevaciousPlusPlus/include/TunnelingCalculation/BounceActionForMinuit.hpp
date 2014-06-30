/*
 * BounceActionForMinuit.hpp
 *
 *  Created on: Jun 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BOUNCEACTIONFORMINUIT_HPP_
#define BOUNCEACTIONFORMINUIT_HPP_

#include "CommonIncludes.hpp"
#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "Eigen/Dense"
#include "boost/numeric/odeint/integrate/integrate.hpp"
#include "boost/math/constants/constants.hpp"
#include "PotentialEvaluation.hpp"
#include "PathFromNodes.hpp"
#include "PathFieldsAndPotential.hpp"
#include "SplinePotential.hpp"
#include "BubbleProfile.hpp"

namespace VevaciousPlusPlus
{

  class BounceActionForMinuit : public ROOT::Minuit2::FCNBase
  {
  public:
    BounceActionForMinuit( PotentialFunction const& potentialFunction,
                           size_t const numberOfVaryingPathNodes,
                           size_t const numberOfSplinesInPotential,
                           PotentialMinimum const& falseVacuum,
                           PotentialMinimum const& trueVacuum,
                           double const falseVacuumEvaporationTemperature,
                           size_t const undershootOvershootAttempts = 32,
                           size_t const maximumMultipleOfLongestLength = 16,
                         double const initialFractionOfShortestLength = 0.05,
                           double const minimumScaleSquared = 1.0,
                          double const shootingCloseEnoughThreshold = 0.01 );
    virtual
    ~BounceActionForMinuit();

    // The bounce action is calculated by the following process:
    // 1) The path in field space from the false vacuum to the true vacuum is
    //    decoded from pathParameterization to give the field configuration f
    //    as a function of a path auxiliary variable p, giving f(p) and df/dp.
    //    This is obtained from pathFromNodes.
    // 2) The potential is fitted as a polynomial in p of degree
    //    numberOfSplinesInPotential, giving V(p).
    // 3) The one-dimensional bubble equation of motion along p is integrated
    //    to get p as a function of the radial variable r, which is the spatial
    //    radius for three-dimensional actions, or the length of the
    //    four-dimensional Euclidean vector. This gets the correct bubble
    //    profile only if the transverse equations of motion in field space are
    //    also satisfied, or for a modified potential which has additional
    //    terms that raise the energy barrier on the side of the path that has
    //    a lower energy barrier. Hence this bubble profile gives an upper
    //    bound on the bounce action for the unmodified potential, as the
    //    modified potential (which has the same true vacuum and false vacuum
    //    and is exactly the same along the path) cannot have a lower bounce
    //    action.
    //    [The p equation of motion has a strange term giving a force
    //    proportional to (dp/dr)^2 which is not a friction term as the sign
    //    is not proportional to dp/dr, rather to d^2f/dp^2. This is because p
    //    is not linear in the distance in field space, so actually this term
    //    behaves a bit like extra kinetic energy because it rolled further
    //    "down the hill" from the true vacuum initially because there is more
    //    field space distance covered by a small change in p if the derivative
    //    is that way, or like extra friction because the approach to the
    //    false vacuum is actually longer in field space distance than in p.
    //    The alternative is to divide up the path into pathResolution segments
    //    and take the path as a series of straight lines at a constant
    //    "velocity" in field space, meaning that the weird pseudo-friction
    //    force in the bubble wall equation of motion in p disappears. This
    //    would require some contortions with linked lists (probably) of
    //    segments to get f(p(r)) and df/dp.]
    // 4) The bounce action along p(r) is numerically integrated and returned.
    virtual double
    operator()( std::vector< double > const& pathParameterization ) const;

    // This implements Up() for FCNBase just to stick to a basic value.
    virtual double Up() const { return 1.0; }

    // This sets up initialParameterization and initialStepSizes to be a
    // straight path in field space with initial step sizes of stepSizeFraction
    // times the differences in the fields between true vacuum and false
    // vacuum. After the node parameters are added to initialParameterization,
    // startingTemperature is appended if it is > 0.0. There is no return value
    // optimization because we cannot be sure that poor physicist users will
    // have access to a C++11-compliant compiler.
    void
    SetUpStraightPathForMinuit( std::vector< double >& initialParameterization,
                                std::vector< double >& initialStepSizes,
                                double const startingTemperature,
                                double const stepSizeFraction = 0.2 ) const;

    // This plots the fields as functions of the bubble radial value in a file
    // called plotFilename in .eps format, with each field plotted in the color
    // given by fieldColors: the field with index i is plotted in the color
    // given by fieldColors[ i ]. An empty string indicates that the field
    // should not be plotted.
    void PlotBubbleProfile( std::vector< double > const& pathParameterization,
                            std::vector< std::string > const& fieldColors,
                            std::string const& plotFilename ) const;


  protected:
    static double const radiusDifferenceThreshold;

    PotentialFunction const& potentialFunction;
    size_t const numberOfFields;
    size_t referenceFieldIndex;
    PathFromNodes pathFromNodes;
    size_t numberOfSplinesInPotential;
    PotentialMinimum const& falseVacuum;
    PotentialMinimum const& trueVacuum;
    // The vector in field space from falseVacuum to trueVacuum (which are
    // zero-temperature vacua) is given by zeroTemperatureStraightPath.
    std::vector< double > zeroTemperatureStraightPath;
    double zeroTemperatureStraightPathInverseLengthSquared;
    double const falseVacuumPotential;
    double const trueVacuumPotential;
    double const falseVacuumEvaporationTemperature;
    double const tunnelingScaleSquared;
    double shortestLength;
    double longestLength;
    size_t const undershootOvershootAttempts;
    double const initialFractionOfShortestLength;
    double const shootingThreshold;

    // This returns true if pathParameterization.size() is only
    // pathFromNodes.ParameterizationSize() (hence T = 0.0 is implicit) or if
    // pathParameterization.back() <= 0.0; otherwise false is returned, or an
    // exception is thrown if pathParameterization.size() is not
    // ( pathFromNodes.ParameterizationSize() + 1 ).
    bool ZeroTemperatureParameterization(
                     std::vector< double > const& pathParameterization ) const;

    // This turns pathParameterization into a PathFieldsAndPotential, by first
    // checking for a non-zero temperature, then setting up the straight path
    // in field space, and projecting the nodes extracted from
    // pathParameterization onto planes perpendicular to the straight path. A
    // few extra bits of information to do with the temperature are also
    // recorded in the PathFieldsAndPotential object returned.
    PathFieldsAndPotential DecodePathParameters(
                     std::vector< double > const& pathParameterization ) const;

    // This evaluates the bounce action density at the given point on the
    // bubble profile.
    double
    BounceActionDensity( PathFieldsAndPotential const& pathFieldsAndPotential,
                      BubbleRadialValueDescription const& profilePoint ) const;

    // This puts a polynomial approximation of the potential along the path
    // given by pathFieldsAndPotential into pathFieldsAndPotential.
    void
    PotentialAlongPath( PathFieldsAndPotential& pathFieldsAndPotential ) const;

    // This sets up the bubble profile, numerically integrates the bounce
    // action over it, and then returns the bounce action [S_4 or S_3(T)].
    double
    BounceAction( PathFieldsAndPotential const& pathFieldsAndPotential ) const;
  };




  // This returns either S_4 or S_3(T), where S_4 is the dimensionless
  // four-dimensional bounce action and S_3(T) the three-dimensional bounce
  // action in units of GeV for the given temperature T, also in GeV. If T is
  // greater than 0, S_3(T) is returned, otherwise S_4 is returned. If no
  // temperature is given, T is assumed to be 0.
  // The bounce action is calculated by the following process:
  // 1) The path in field space from the false vacuum to the true vacuum is
  //    decoded from pathParameterization to give the field configuration f as
  //    a function of a path auxiliary variable p, giving f(p) and df/dp. See
  //    the comment above DecodePathParameters for the format of
  //    pathParameterization.
  // 2) The potential is fitted as a polynomial in p of degree
  //    numberOfSplinesInPotential, giving V(p).
  // 3) The one-dimensional bubble equation of motion along p is integrated to
  //    get p as a function of the radial variable r, which is the spatial
  //    radius for three-dimensional actions, or the length of the
  //    four-dimensional Euclidean vector. This gets the correct bubble profile
  //    only if the transverse equations of motion in field space are also
  //    satisfied, or for a modified potential which has additional terms that
  //    raise the energy barrier on the side of the path that has a lower
  //    energy barrier. Hence this bubble profile gives an upper bound on the
  //    bounce action for the unmodified potential, as the modified potential
  //    (which has the same true vacuum and false vacuum and is exactly the
  //    same along the path) cannot have a lower bounce action.
  //    [The p equation of motion has a strange term giving a force
  //    proportional to (dp/dr)^2 which is not a friction term as the sign
  //    is not proportional to dp/dr, rather to d^2f/dp^2. This is because p
  //    is not linear in the distance in field space, so actually this term
  //    behaves a bit like extra kinetic energy because it rolled further
  //    "down the hill" from the true vacuum initially because there is more
  //    field space distance covered by a small change in p if the derivative
  //    is that way, or like extra friction because the approach to the
  //    false vacuum is actually longer in field space distance than in p.
  //    The alternative is to divide up the path into pathResolution segments
  //    and take the path as a series of straight lines at a constant
  //    "velocity" in field space, meaning that the weird pseudo-friction
  //    force in the bubble wall equation of motion in p disappears. This
  //    would require some contortions with linked lists (probably) of
  //    segments to get f(p(r)) and df/dp.]
  // 4) The bounce action along p(r) is numerically integrated and returned.
  inline double BounceActionForMinuit::operator()(
                      std::vector< double > const& pathParameterization ) const
  {
    PathFieldsAndPotential
    pathFieldsAndPotential( DecodePathParameters( pathParameterization ) );
    PotentialAlongPath( pathFieldsAndPotential );
    return BounceAction( pathFieldsAndPotential );
  }

  // This evaluates the bounce action density at the given point on the
  // bubble profile.
  inline double BounceActionForMinuit::BounceActionDensity(
                          PathFieldsAndPotential const& pathFieldsAndPotential,
                       BubbleRadialValueDescription const& profilePoint ) const
  {
    double const currentAuxiliary( profilePoint.auxiliaryValue );
    double kineticTerm( profilePoint.auxiliarySlope );
    kineticTerm *= ( 0.5 * kineticTerm
        * pathFieldsAndPotential.FieldDerivativesSquared( currentAuxiliary ) );
    return
    ( kineticTerm
      + pathFieldsAndPotential.PotentialApproximation( currentAuxiliary ) );
  }

  // This sets up initialParameterization and initialStepSizes to be a
  // straight path in field space with initial step sizes of stepSizeFraction
  // times the differences in the fields between true vacuum and false
  // vacuum. After the node parameters are added to initialParameterization,
  // startingTemperature is appended if it is > 0.0. There is no return value
  // optimization because we cannot be sure that poor physicist users will have
  // access to a C++11-compliant compiler.
  inline void BounceActionForMinuit::SetUpStraightPathForMinuit(
                                std::vector< double >& initialParameterization,
                                       std::vector< double >& initialStepSizes,
                                              double const startingTemperature,
                                          double const stepSizeFraction ) const
  {
    pathFromNodes.InitialStepsForMinuit( initialStepSizes,
                                         zeroTemperatureStraightPath,
                                         stepSizeFraction );
    initialParameterization.assign( pathFromNodes.ParameterizationSize(),
                                    0.0 );
    if( startingTemperature > 0.0 )
    {
      initialParameterization.push_back( startingTemperature );
      initialStepSizes.push_back( stepSizeFraction * startingTemperature );
    }
  }

  // This returns true if pathParameterization.size() is only
  // pathFromNodes.ParameterizationSize() (hence T = 0.0 is implicit) or if
  // pathParameterization.back() <= 0.0; otherwise false is returned, or an
  // exception is thrown if pathParameterization.size() is not
  // ( pathFromNodes.ParameterizationSize() + 1 ).
  inline bool BounceActionForMinuit::ZeroTemperatureParameterization(
                      std::vector< double > const& pathParameterization ) const
  {
    if( pathParameterization.size() == pathFromNodes.ParameterizationSize() )
    {
      return true;
    }
    if( pathParameterization.size()
        == ( pathFromNodes.ParameterizationSize() + 1 ) )
    {
      return ( pathParameterization.back() <= 0.0 );
    }
    std::stringstream errorStream;
    errorStream << "BounceActionForMinuit::operator() was given a wrong"
    << " number of parameters: " << pathParameterization.size()
    << "; it should have been given "
    << pathFromNodes.NumberOfVaryingPathNodes() << " nodes each of "
    << pathFromNodes.NumberOfParameterizationFields()
    << " fields, optionally 1 extra parameter for the temperature at the"
    << " end, so " << pathFromNodes.ParameterizationSize() << " or "
    << ( pathFromNodes.ParameterizationSize() + 1 ) << " parameters.";
    throw std::range_error( errorStream.str() );
  }

} /* namespace VevaciousPlusPlus */
#endif /* BOUNCEACTIONFORMINUIT_HPP_ */
