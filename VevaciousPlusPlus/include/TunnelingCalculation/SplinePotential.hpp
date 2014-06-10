/*
 * SplinePotential.hpp
 *
 *  Created on: Jun 10, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SPLINEPOTENTIAL_HPP_
#define SPLINEPOTENTIAL_HPP_

namespace VevaciousPlusPlus
{

  class SplinePotential
  {
  public:
    SplinePotential();
    SplinePotential( SplinePotential const& copySource );
    virtual
    ~SplinePotential();


    // This returns the value of the potential at auxiliaryValue, by finding
    // the correct spline and then returning its value at that point.
    virtual double operator()( double auxiliaryValue ) const;

    // This returns the value of the first derivative of the potential at
    // auxiliaryValue, by finding the correct spline and then returning its
    // slope at that point.
    virtual double FirstDerivative( double const auxiliaryValue ) const;

    // This adds another point for the spline, assuming that it goes after the
    // previously-added point by auxiliaryDifference from the previous point,
    // to a potential value of potentialValue relative to the false vacuum.
    void AddPoint( double const auxiliaryDifference,
                   double const potentialValue )
    { auxiliaryValues.push_back( auxiliaryDifference );
      potentialValues.push_back( potentialValue );
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "SplinePotential::AddPoint( " << auxiliaryDifference << ", "
    << potentialValue << " ) called.";
    std::cout << std::endl;/**/}

    // This sets up the splines based on auxiliaryValues and potentialValues,
    // ensuring that the potential reaches the correct values for the vacua,
    // and that the potential derivative vanishes at the vacua.
    void SetSplines( double const trueVacuumPotentialDifference );

    // This is for debugging.
    std::string AsDebuggingString() const;


  protected:
    // There are auxiliaryValues.size() normal splines, each represented by 4
    // values: the size of the spline in the auxiliary value, the potential
    // (relative to the false vacuum at zero auxiliary value), its first
    // derivative, and its second derivative. Within the spline, the potential
    // is approximated by
    // V(p_j + d) + d * V'(p_j) + d^2 * V''(p_j),
    // where the auxiliary value is p = p_j + d. There is also one final
    // spline, where for an auxiliary value = 1 - d, the potential is
    // approximated by
    // V(1.0) + d^2 * V''(1) + d^3 * V'''(1).
    std::vector< double > auxiliaryValues;
    std::vector< double > potentialValues;
    std::vector< double > firstDerivatives;
    std::vector< double > halfSecondDerivatives;
    double finalPotential;
    double halfFinalSecondDerivative;
    double finalCubicCoefficient;
  };

} /* namespace VevaciousPlusPlus */
#endif /* SPLINEPOTENTIAL_HPP_ */
