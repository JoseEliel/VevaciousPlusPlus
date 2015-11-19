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
    BOL::BasicObserver(),
    minimumScaleSquared( -1.0 )
  {
    lagrangianParameterManager.registerObserver( this );
  }

  RgeImprovedOneLoopPotential::RgeImprovedOneLoopPotential(
                   PotentialFromPolynomialWithMasses const& potentialToCopy ) :
    PotentialFromPolynomialWithMasses( potentialToCopy ),
    BOL::BasicObserver(),
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

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "RgeImprovedOneLoopPotential::operator()( "
    << FieldConfigurationAsMathematica( fieldConfiguration )
    << ", T =" << temperatureValue
    << " ) called. Before comparing, scaleSquared = " << scaleSquared
    << ", minimumScaleSquared = " << minimumScaleSquared;
    std::cout << std::endl;*/

    // Of course, if minimumScaleSquared is *greater than* scaleSquared then we
    // need to use minimumScaleSquared as scaleSquared is below the minimum
    // allowed, hence using std::max.
    scaleSquared = std::max( scaleSquared, minimumScaleSquared );
    // The logarithm of the scale is of course half the logarithm of the square
    // of the scale.
    std::vector< double > const
    parameterValues( lagrangianParameterManager.ParameterValues(
                                             ( 0.5 * log( scaleSquared ) ) ) );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "parameterValues = { ";
    for( std::vector< double >::const_iterator
         parameterValue( parameterValues.begin() );
         parameterValue < parameterValues.end();
         ++parameterValue )
    {
      std::cout << *parameterValue << ", ";
    }
    std::cout << "}, treeLevelPotential = "
    << treeLevelPotential.AsDebuggingString();
    std::vector< double > fieldOrigin( FieldValuesOrigin() );
    std::vector< DoubleVectorWithDouble > originScalarMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( parameterValues,
                                      fieldOrigin,
                                      scalarSquareMasses,
                                      originScalarMassesSquaredWithFactors );
    std::vector< DoubleVectorWithDouble >
    originFermionMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( parameterValues,
                                      fieldOrigin,
                                      fermionSquareMasses,
                                      originFermionMassesSquaredWithFactors );
    std::vector< DoubleVectorWithDouble > originVectorMassesSquaredWithFactors;
    AddMassesSquaredWithMultiplicity( parameterValues,
                                      fieldOrigin,
                                      vectorSquareMasses,
                                      originVectorMassesSquaredWithFactors );
    std::cout << "For origin:"
    << std::endl << "originScalarMassesSquaredWithFactors = {" << std::endl;
    for( std::vector< DoubleVectorWithDouble >::const_iterator
         massesWithFactor( originScalarMassesSquaredWithFactors.begin() );
         massesWithFactor < originScalarMassesSquaredWithFactors.end();
         ++massesWithFactor )
    {
      std::cout
      << "[ factor: " << massesWithFactor->second
      << ", masses-squared: { ";
      for( std::vector< double >::const_iterator
           massSquared( massesWithFactor->first.begin() );
           massSquared < massesWithFactor->first.end();
           ++massSquared )
      {
        std::cout << *massSquared << ",";
      }
      std::cout << " } ]" << std::endl;
    }
    std::cout << "}"
    << std::endl << "originFermionMassesSquaredWithFactors = {" << std::endl;
    for( std::vector< DoubleVectorWithDouble >::const_iterator
         massesWithFactor( originFermionMassesSquaredWithFactors.begin() );
         massesWithFactor < originFermionMassesSquaredWithFactors.end();
         ++massesWithFactor )
    {
      std::cout
      << "[ factor: " << massesWithFactor->second
      << ", masses-squared: { ";
      for( std::vector< double >::const_iterator
           massSquared( massesWithFactor->first.begin() );
           massSquared < massesWithFactor->first.end();
           ++massSquared )
      {
        std::cout << *massSquared << ",";
      }
      std::cout << " } ]" << std::endl;
    }
    std::cout << "}"
    << std::endl << "originVectorMassesSquaredWithFactors = {" << std::endl;
    for( std::vector< DoubleVectorWithDouble >::const_iterator
         massesWithFactor( originVectorMassesSquaredWithFactors.begin() );
         massesWithFactor < originVectorMassesSquaredWithFactors.end();
         ++massesWithFactor )
    {
      std::cout
      << "[ factor: " << massesWithFactor->second
      << ", masses-squared: { ";
      for( std::vector< double >::const_iterator
           massSquared( massesWithFactor->first.begin() );
           massSquared < massesWithFactor->first.end();
           ++massSquared )
      {
        std::cout << *massSquared << ",";
      }
      std::cout << " } ]" << std::endl;
    }
    std::cout << "}"
    << std::endl
    << "Before fixed scale update,"
    << " vectorMassSquaredMatrices.front().DebugCurrentValues() with"
    << " parameters = " << std::endl
    << vectorMassSquaredMatrices.front().DebugCurrentValues( parameterValues,
                                                             fieldOrigin )
    << std::endl
    << " vectorMassSquaredMatrices.front().DebugCurrentValues() without"
    << " parameters = " << std::endl
    << vectorMassSquaredMatrices.front().DebugCurrentValues( fieldOrigin )
    << std::endl;
    std::cout
    << " vectorMassSquaredMatrices.front().AsString() = { "
    << vectorMassSquaredMatrices.front().AsString()
    << std::endl;
    std::vector< RealMassesSquaredMatrix >
    updatedVectorMassSquareds( vectorMassSquaredMatrices );
    std::vector< MassesSquaredCalculator* > updatedPointers;
    for( std::vector< RealMassesSquaredMatrix >::iterator
         massMatrix( updatedVectorMassSquareds.begin() );
         massMatrix < updatedVectorMassSquareds.end();
         ++massMatrix )
    {
      massMatrix->UpdateForFixedScale( parameterValues );
      updatedPointers.push_back( &(*massMatrix) );
    }
    std::cout << std::endl
    << "After fixed scale update,"
    << " updatedVectorMassSquareds.front().DebugCurrentValues() with"
    << " parameters = " << std::endl
    << updatedVectorMassSquareds.front().DebugCurrentValues( parameterValues,
                                                             fieldOrigin )
    << std::endl
    << " updatedVectorMassSquareds.front().DebugCurrentValues() without"
    << " parameters = " << std::endl
    << updatedVectorMassSquareds.front().DebugCurrentValues( fieldOrigin )
    << std::endl;
    std::cout
    << " updatedVectorMassSquareds.front().AsString() = { "
    << updatedVectorMassSquareds.front().AsString()
    << std::endl;
    originVectorMassesSquaredWithFactors.clear();
    AddMassesSquaredWithMultiplicity( parameterValues,
                                      fieldOrigin,
                                      updatedPointers,
                                      originVectorMassesSquaredWithFactors );
    std::cout
    << "After fixed scale update, originVectorMassesSquaredWithFactors with"
    << " parameterValues passed in = {"
    << std::endl;
    for( std::vector< DoubleVectorWithDouble >::const_iterator
         massesWithFactor( originVectorMassesSquaredWithFactors.begin() );
         massesWithFactor < originVectorMassesSquaredWithFactors.end();
         ++massesWithFactor )
    {
      std::cout
      << "[ factor: " << massesWithFactor->second
      << ", masses-squared: { ";
      for( std::vector< double >::const_iterator
           massSquared( massesWithFactor->first.begin() );
           massSquared < massesWithFactor->first.end();
           ++massSquared )
      {
        std::cout << *massSquared << ",";
      }
      std::cout << " } ]" << std::endl;
    }
    std::cout << "}" << std::endl;
    originVectorMassesSquaredWithFactors.clear();
    AddMassesSquaredWithMultiplicity( fieldOrigin,
                                      updatedPointers,
                                      originVectorMassesSquaredWithFactors );
    std::cout
    << "After fixed scale update, originVectorMassesSquaredWithFactors without"
    << " parameterValues passed in = {"
    << std::endl;
    for( std::vector< DoubleVectorWithDouble >::const_iterator
         massesWithFactor( originVectorMassesSquaredWithFactors.begin() );
         massesWithFactor < originVectorMassesSquaredWithFactors.end();
         ++massesWithFactor )
    {
      std::cout
      << "[ factor: " << massesWithFactor->second
      << ", masses-squared: { ";
      for( std::vector< double >::const_iterator
           massSquared( massesWithFactor->first.begin() );
           massSquared < massesWithFactor->first.end();
           ++massSquared )
      {
        std::cout << *massSquared << ",";
      }
      std::cout << " } ]" << std::endl;
    }
    std::cout << " }" << std::endl;
    std::cout << std::endl;*/

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
