/*
 * RgeImprovedOneLoopPotential.cpp
 *
 *  Created on: Nov 5, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/PotentialFunctions/RgeImprovedOneLoopPotential.hpp"

namespace VevaciousPlusPlus
{

  RgeImprovedOneLoopPotential::RgeImprovedOneLoopPotential(
                                              std::string const& modelFilename,
                               double const assumedPositiveOrNegativeTolerance,
                     LagrangianParameterManager& lagrangianParameterManager ) :
    PotentialFromPolynomialWithMasses( modelFilename,
                                       assumedPositiveOrNegativeTolerance,
                                       lagrangianParameterManager ),
    BOL::PushedToObserver< LagrangianParameterManager >(),
    minimumScaleSquared( -1.0 )
  {
    lagrangianParameterManager.registerObserver( this );
  }

  RgeImprovedOneLoopPotential::RgeImprovedOneLoopPotential(
                   PotentialFromPolynomialWithMasses const& potentialToCopy ) :
    PotentialFromPolynomialWithMasses( potentialToCopy ),
    minimumScaleSquared( -1.0 )
  {
    lagrangianParameterManager.registerObserver( this );
  }

  RgeImprovedOneLoopPotential::~RgeImprovedOneLoopPotential()
  {
    // This does nothing.
  }


  // This returns the energy density in GeV^4 of the potential for a state
  // strongly peaked around expectation values (in GeV) for the fields given
  // by the values of fieldConfiguration and temperature in GeV given by
  // temperatureValue.
  double RgeImprovedOneLoopPotential::operator()(
                               std::vector< double > const& fieldConfiguration,
                                          double const temperatureValue ) const
  {
    double scaleSquared( temperatureValue * temperatureValue );
    for( std::vector< double >::const_iterator
         fieldValue( fieldConfiguration.begin() );
         fieldValue < fieldConfiguration.end();
         ++fieldValue )
    {
      scaleSquared += ( (*fieldValue) * (*fieldValue) );
    }

    // Of course, if minimumScaleSquared is *greater than* scaleSquared then we
    // need to use minimumScaleSquared as scaleSquared is below the minimum
    // allowed, hence using std::max.
    scaleSquared = std::max( scaleSquared, minimumScaleSquared );
    // The logarithm of the scale is of course half the logarithm of the square
    // of the scale.
    std::vector< double > const
    parameterValues( lagrangianParameterManager.ParameterValues(
                                             ( 0.5 * log( scaleSquared ) ) ) );
    std::vector< DoubleVectorWithDouble > scalarMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( parameterValues,
                                      fieldConfiguration,
                                      scalarSquareMasses,
                                      scalarMassesSquaredWithFactors );
    std::vector< DoubleVectorWithDouble > fermionMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( parameterValues,
                                      fieldConfiguration,
                                      fermionSquareMasses,
                                      fermionMassesSquaredWithFactors );
    std::vector< DoubleVectorWithDouble > vectorMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( parameterValues,
                                      fieldConfiguration,
                                      vectorSquareMasses,
                                      vectorMassesSquaredWithFactors );
    return ( treeLevelPotential( parameterValues,
                                 fieldConfiguration )
             + polynomialLoopCorrections( parameterValues,
                                          fieldConfiguration )
             + LoopAndThermalCorrections( scalarMassesSquaredWithFactors,
                                          fermionMassesSquaredWithFactors,
                                          vectorMassesSquaredWithFactors,
                                          ( 1.0 / scaleSquared ),
                                          temperatureValue ) );
  }

  // This returns a string that is valid Python with no indentation to evaluate
  // the potential in three functions:
  // TreeLevelPotential( fv ), JustLoopCorrectedPotential( fv ), and
  // LoopAndThermallyCorrectedPotential( fv ).
   std::string RgeImprovedOneLoopPotential::WriteActualPythonFunction() const
   {
     std::stringstream stringBuilder;
     stringBuilder << std::setprecision( 12 );
     stringBuilder << "def AppropriateScaleSquared( fv ):\n"
     "  return max( " << minimumScaleSquared << ",\n"
     "              ( temperatureSquare + sum( f**2 for f in fv ) ) )\n"
     "\n"
     "def TreeLevelPotential( fv ):\n"
     "  scaleSquared = AppropriateScaleSquared( fv )\n"
     "  lp = LagrangianParameters( 0.5 * math.log( scaleSquared ) )\n"
     "  return TreeLevelContribution( fv, lp )\n"
     "\n"
     "def JustLoopCorrectedPotential( fv ):\n"
     "  scaleSquared = AppropriateScaleSquared( fv )\n"
     "  lp = LagrangianParameters( 0.5 * math.log( scaleSquared ) )\n"
     "  return ( TreeLevelContribution( fv, lp )\n"
     "           + PolynomialLoopCorrections( fv, lp )\n"
     "           + JustLoopCorrections( fv, lp, ( 1.0 / scaleSquared ) ) )\n"
     "\n"
     "def LoopAndThermallyCorrectedPotential( fv ):\n"
     "  scaleSquared = AppropriateScaleSquared( fv )\n"
     "  lp = LagrangianParameters( 0.5 * math.log( scaleSquared ) )\n"
     "  return ( TreeLevelContribution( fv, lp )\n"
     "           + PolynomialLoopCorrections( fv, lp )\n"
     "           + LoopAndThermalCorrections( fv,\n"
     "                                        lp,\n"
     "                                        ( 1.0 / scaleSquared ),\n"
     "                                        temperatureInverseSquare ) )\n";
     return stringBuilder.str();
   }

} /* namespace VevaciousPlusPlus */
