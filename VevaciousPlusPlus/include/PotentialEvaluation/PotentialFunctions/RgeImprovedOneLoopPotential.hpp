/*
 * RgeImprovedOneLoopPotential.hpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef RGEIMPROVEDONELOOPPOTENTIAL_HPP_
#define RGEIMPROVEDONELOOPPOTENTIAL_HPP_

#include "CommonIncludes.hpp"
#include "PotentialFromPolynomialAndMasses.hpp"
#include "SlhaManagement/RunningParameterManager.hpp"
#include "SlhaManagement/SlhaManager.hpp"
#include "PotentialEvaluation/MassesSquaredCalculator.hpp"
#include "PotentialMinimization/PotentialMinimum.hpp"
#include "PotentialMinimization/HomotopyContinuation/FieldPolynomialsWithScale.hpp"

namespace VevaciousPlusPlus
{

  class RgeImprovedOneLoopPotential : public PotentialFromPolynomialAndMasses
  {
  public:
    RgeImprovedOneLoopPotential( std::string const& modelFilename,
                                 double const scaleRangeMinimumFactor,
            bool const treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions,
                            RunningParameterManager& runningParameterManager );
    RgeImprovedOneLoopPotential(
    PotentialFromPolynomialAndMasses const& potentialFromPolynomialAndMasses );
    virtual ~RgeImprovedOneLoopPotential();


    // This returns the energy density in GeV^4 of the potential for a state
    // strongly peaked around expectation values (in GeV) for the fields given
    // by the values of fieldConfiguration and temperature in GeV given by
    // temperatureValue.
    virtual double operator()( std::vector< double > const& fieldConfiguration,
                               double const temperatureValue = 0.0 ) const;

    // This returns the tree-level potential energy density evaluated at the
    // correct scale.
    virtual double
    QuickApproximation( std::vector< double > const& fieldConfiguration,
                        double const temperatureValue = 0.0 )
    { return treeLevelPotential( fieldConfiguration,
                  ( 0.5 * log( RenormalizationScaleSquared( fieldConfiguration,
                                                    temperatureValue ) ) ) ); }

    // This returns the square of the Euclidean distance between the two vacua.
    virtual double
    ScaleSquaredRelevantToTunneling( PotentialMinimum const& falseVacuum,
                                    PotentialMinimum const& trueVacuum ) const;

    // This performs all relevant updates for the new SLHA data except for
    // propagating the push to the set of dependent SlhaUpdatePropagators.
    virtual void UpdateSelfForNewSlha( SlhaManager const& slhaManager );

    virtual PolynomialGradientTargetSystem* HomotopyContinuationTargetSystem()
    { return &homotopyContinuationTargetSystem; }


  protected:
    double logarithmOfMinimumRenormalizationScale;
    double logarithmOfMaximumRenormalizationScale;
    FieldPolynomialsWithScale homotopyContinuationTargetSystem;


    // This returns the square of an appropriate renormalization scale.
    double RenormalizationScaleSquared(
                               std::vector< double > const& fieldConfiguration,
                                   double const temperatureValue = 0.0 ) const;

    // This appends the masses-squared and multiplicity from each
    // MassesSquaredFromMatrix in massSquaredMatrices to
    // massSquaredMatrices, with all functionoids evaluated at the natural
    // exponent of logarithmOfScale.
    void AddMassesSquaredWithMultiplicity(
                               std::vector< double > const& fieldConfiguration,
            std::vector< MassesSquaredCalculator* > const& massSquaredMatrices,
               std::vector< DoubleVectorWithDouble >& massesSquaredWithFactors,
                                         double const logarithmOfScale ) const;

    // This should return a string that is valid Python indented by 4 spaces to
    // set the scale Q, given by Q^(-2) called invQSq, for an evaluation of the
    // potential assuming that the field configuration is given as an array
    // called "fv" and the temperature by T, given by T^(-2) called invTSq,
    // given the rest of the Python code written by WriteAsPython. By default Q
    // is left as the lowest scale given by the blocks of the SLHA file.
    virtual std::string SetScaleForPythonPotentialCall() const;
  };




  // This returns the square of an appropriate renormalization scale.
  inline double RgeImprovedOneLoopPotential::RenormalizationScaleSquared(
                               std::vector< double > const& fieldConfiguration,
                                          double const temperatureValue ) const
  {
    double renormalizationScaleSquared( squareOfMinimumRenormalizationScale
                                   + ( temperatureValue * temperatureValue ) );
    for( std::vector< double >::const_iterator
         whichField( fieldConfiguration.begin() );
         whichField < fieldConfiguration.end();
         ++whichField )
    {
      renormalizationScaleSquared += ( (*whichField) * (*whichField) );
    }
    return renormalizationScaleSquared;
  }

  // This appends the masses-squared and multiplicity from each
  // MassesSquaredFromMatrix in massSquaredMatrices to
  // massSquaredMatrices, with all functionoids evaluated at the last scale
  // which was used to update them.
  inline void
  RgeImprovedOneLoopPotential::AddMassesSquaredWithMultiplicity(
                               std::vector< double > const& fieldConfiguration,
            std::vector< MassesSquaredCalculator* > const& massSquaredMatrices,
               std::vector< DoubleVectorWithDouble >& massesSquaredWithFactors,
                                          double const logarithmOfScale ) const
  {
    for( std::vector< MassesSquaredCalculator* >::const_iterator
         whichMatrix( massSquaredMatrices.begin() );
         whichMatrix < massSquaredMatrices.end();
         ++whichMatrix )
    {
      massesSquaredWithFactors.push_back(
             std::make_pair( (*whichMatrix)->MassesSquared( fieldConfiguration,
                                                            logarithmOfScale ),
                                      (*whichMatrix)->MultiplicityFactor() ) );
    }
  }

  // This returns the square of the Euclidean distance between the two vacua.
  inline double RgeImprovedOneLoopPotential::ScaleSquaredRelevantToTunneling(
                                           PotentialMinimum const& falseVacuum,
                                     PotentialMinimum const& trueVacuum ) const
  {
    return falseVacuum.SquareDistanceTo( trueVacuum );
  }

  // This returns a string that is valid Python indented by 4 spaces to set the
  // scale Q, given by Q^(-2) called invQSq, for an evaluation of the potential
  // assuming that the field configuration is given as an array called "fv" and
  // the temperature by T, given by T^(-2) called invTSq, given the rest of the
  // Python code written by WriteAsPython. By default Q is left as the lowest
  // scale given by the blocks of the SLHA file.
  inline std::string
  RgeImprovedOneLoopPotential::SetScaleForPythonPotentialCall() const
  {
    std::stringstream stringBuilder;
    stringBuilder
    << "# We use the square of lowest scale of the SLHA blocks plus the sum\n"
    "# of the squares of the fields plus the square of the temperature as\n"
    "# the square of the scale.\n"
    "    QQ = (" << currentMinimumRenormalizationScale << ")**2\n"
    "    for f in fv:\n"
    "        QQ += f**2\n"
    "    if ( invTSq > 0.0 ):\n"
    "        QQ += ( 1.0 / invTSq )\n"
    "    global invQSq"
    "    invQSq = ( 1.0 / QQ )\n"
    "    UpdateRunningParameters( 0.5 * math.log( QQ ) )\n";
    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* RGEIMPROVEDONELOOPPOTENTIAL_HPP_ */
