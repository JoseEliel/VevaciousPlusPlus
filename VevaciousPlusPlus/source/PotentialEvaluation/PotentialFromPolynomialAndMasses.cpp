/*
 * PotentialFromPolynomialAndMasses.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{

  PotentialFromPolynomialAndMasses::PotentialFromPolynomialAndMasses(
                                           std::string const& modelFilename ) :
    HomotopyContinuationReadyPotential(),
    runningParameters(),
    treeLevelPotential(),
    massSquaredMatrices(),
    vectorMassCorrectionConstant( NAN ),
    polynomialGradient(),
    polynomialHessian(),
    scaleSlopeOfGradient()
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MassCorrectedPotential::MassCorrectedPotential( \""
    << modelFilename << "\" )";
    std::cout << std::endl;/**/
  }

  PotentialFromPolynomialAndMasses::~PotentialFromPolynomialAndMasses()
  {
    // This does nothing.
  }


  double PotentialFromPolynomialAndMasses::operator()(
                               std::vector< double > const& fieldConfiguration,
                                             double const temperatureValue )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MassCorrectedPotential::operator()(...)";
    std::cout << std::endl;

    return 0.0;/**/
  }

  // This updates all the parameters of the potential that are not field
  // values based on the values that appear in blocks in the SLHA format in
  // the file given by slhaFilename.
  void PotentialFromPolynomialAndMasses::UpdateParameters(
                                              std::string const& slhaFilename )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "MassCorrectedPotential::UpdateParameters( \"" << slhaFilename
    << "\" )";
    std::cout << std::endl;/**/
  }

  // This returns the square of the scale (in GeV^2) relevant to tunneling
  // between the given minima for this potential.
  double PotentialFromPolynomialAndMasses::ScaleSquaredRelevantToTunneling(
                                           PotentialMinimum const& falseVacuum,
                                           PotentialMinimum const& trueVacuum )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "PotentialFromPolynomialAndMasses::ScaleSquaredRelevantToTunneling("
    << " ... )";
    std::cout << std::endl;
    return 0.0;/**/
  }

  // This evaluates the target system and places the values in
  // destinationVector.
  void PotentialFromPolynomialAndMasses::HomotopyContinuationSystemValues(
                             std::vector< double > fieldConfigurationWithScale,
                                     std::vector< double >& destinationVector )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "PotentialFromPolynomialAndMasses::HomotopyContinuationSystemValues("
    << " ... )";
    std::cout << std::endl;/**/
  }

  // This evaluates the derivatives of the target system and places the
  // values in destinationMatrix.
  void PotentialFromPolynomialAndMasses::HomotopyContinuationSystemGradients(
                             std::vector< double > fieldConfigurationWithScale,
                      std::vector< std::vector< double > >& destinationMatrix )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "PotentialFromPolynomialAndMasses::HomotopyContinuationSystemGradients("
    << " ... )";
    std::cout << std::endl;/**/
  }


  PotentialFromPolynomialAndMasses::PotentialFromPolynomialAndMasses() :
    HomotopyContinuationReadyPotential(),
    runningParameters(),
    treeLevelPotential(),
    massSquaredMatrices(),
    vectorMassCorrectionConstant( NAN ),
    polynomialGradient(),
    polynomialHessian(),
    scaleSlopeOfGradient()
  {
    // This protected constructor is just an initialization list only used by
    // derived classes which are going to fill up the data members in their own
    // constructors.
  }


  // This sets the renormalization scale and broadcasts it to the running
  // parameters.
  void PotentialFromPolynomialAndMasses::UpdateRenormalizationScale(
                                      std::vector< double > fieldConfiguration,
                                           double const evaluationTemperature )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "PotentialFromPolynomialAndMasses::UpdateRenormalizationScale( ... )";
    std::cout << std::endl;/**/
  }

  // This evaluates the sum of corrections for a set of real scalar degrees
  // of freedom with masses-squared given by massesSquared at a temperature
  // given by evaluationTemperature.
  double PotentialFromPolynomialAndMasses::ScalarBosonCorrection(
                                    std::vector< double > const& massesSquared,
                                           double const evaluationTemperature )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "PotentialFromPolynomialAndMasses::ScalarBosonCorrection( ... )";
    std::cout << std::endl;
    return 0.0;/**/
  }

  // This evaluates the sum of corrections for a set of Weyl fermion degrees
  // of freedom with masses-squared given by massesSquared at a temperature
  // given by evaluationTemperature.
  double PotentialFromPolynomialAndMasses::WeylFermionCorrection(
                                    std::vector< double > const& massesSquared,
                                           double const evaluationTemperature )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "PotentialFromPolynomialAndMasses::WeylFermionCorrection( ... )";
    std::cout << std::endl;
    return 0.0;/**/
  }

  // This evaluates the sum of corrections for a set of vector gauge boson
  // degrees of freedom (transverse modes) with masses-squared given by
  // massesSquared at a temperature given by evaluationTemperature.
  double PotentialFromPolynomialAndMasses::GaugeBosonCorrection(
                                    std::vector< double > const& massesSquared,
                                           double const evaluationTemperature )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "PotentialFromPolynomialAndMasses::GaugeBosonCorrection( ... )";
    std::cout << std::endl;
    return 0.0;/**/
  }

  // This prepares system of polynomials for the homotopy continuation as a
  // set of polynomials in the field variables with coefficients from the
  // polynomials of the tree-level potential and the polynomial loop
  // corrections. The coefficient of each polynomial term is fitted to a
  // polynomial of degree powerOfScale in the logarithm of the scale. After
  // this, the differentials of the system are derived from these polynomials
  // in the fields and the logarithm of the scale, and also for the
  // constraint relating the scale to the field variables.
  void
  PotentialFromPolynomialAndMasses::PrepareHomotopyContinuationPolynomials(
                                                       int const powerOfScale )
  {
    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "PotentialFromPolynomialAndMasses::"
    << "PrepareHomotopyContinuationPolynomials( ... )";
    std::cout << std::endl;/**/
  }

} /* namespace VevaciousPlusPlus */
