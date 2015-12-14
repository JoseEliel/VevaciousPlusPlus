/*
 * PotentialFromPolynomialWithMasses.cpp
 *
 *  Created on: Oct 09, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/PotentialFunctions/PotentialFromPolynomialWithMasses.hpp"

namespace VevaciousPlusPlus
{
  std::string const PotentialFromPolynomialWithMasses::digitChars(
                                        LHPC::ParsingUtilities::DigitChars() );
  std::string const PotentialFromPolynomialWithMasses::dotAndDigits(
                         PotentialFromPolynomialWithMasses::digitChars + "." );
  std::string const PotentialFromPolynomialWithMasses::allowedVariableInitials(
                               LHPC::ParsingUtilities::UppercaseAlphabetChars()
                          + LHPC::ParsingUtilities::LowercaseAlphabetChars() );
  std::string const PotentialFromPolynomialWithMasses::allowedVariableChars(
                     PotentialFromPolynomialWithMasses::allowedVariableInitials
                                + PotentialFromPolynomialWithMasses::digitChars
                                                                      + "_~" );
  double const PotentialFromPolynomialWithMasses::piSquared(
       boost::math::double_constants::pi * boost::math::double_constants::pi );
  double const
  PotentialFromPolynomialWithMasses::loopFactor( 1.0
                   / ( 64.0 * PotentialFromPolynomialWithMasses::piSquared ) );
  double const PotentialFromPolynomialWithMasses::thermalFactor( 1.0
                    / ( 2.0 * PotentialFromPolynomialWithMasses::piSquared ) );
  std::string const PotentialFromPolynomialWithMasses::positiveByConvention(
                                                      "PositiveByConvention" );
  std::string const PotentialFromPolynomialWithMasses::negativeByConvention(
                                                      "NegativeByConvention" );

  PotentialFromPolynomialWithMasses::PotentialFromPolynomialWithMasses(
                                              std::string const& modelFilename,
                               double const assumedPositiveOrNegativeTolerance,
                     LagrangianParameterManager& lagrangianParameterManager ) :
    PotentialFunction( lagrangianParameterManager ),
    treeLevelPotential(),
    polynomialLoopCorrections(),
    scalarSquareMasses(),
    fermionSquareMasses(),
    vectorSquareMasses(),
    scalarMassSquaredMatrices(),
    fermionMassMatrices(),
    fermionMassSquaredMatrices(),
    vectorMassSquaredMatrices(),
    vectorMassCorrectionConstant( 5.0 / 6.0 ),
    fieldsAssumedPositive(),
    fieldsAssumedNegative(),
    assumedPositiveOrNegativeTolerance( assumedPositiveOrNegativeTolerance )
  {
    LHPC::RestrictedXmlParser xmlParser;

    std::string xmlFieldVariables( "" );
    std::string xmlDsbMinimum( "" );
    std::string xmlTreeLevelPotential( "" );
    std::string xmlLoopCorrections( "" );
    xmlParser.OpenRootElementOfFile( modelFilename );
    while( xmlParser.ReadNextElement() )
    {
      if( xmlParser.CurrentName() == "FieldVariables" )
      {
        xmlFieldVariables = xmlParser.CurrentBody();
      }
      else if( xmlParser.CurrentName() == "DsbMinimum" )
      {
        xmlDsbMinimum = xmlParser.CurrentBody();
      }
      else if( xmlParser.CurrentName() == "TreeLevelPotential" )
      {
        xmlTreeLevelPotential = xmlParser.CurrentBody();
      }
      else if( xmlParser.CurrentName() == "LoopCorrections" )
      {
        std::string
        renormalizationScheme( xmlParser.CurrentAttributes().find(
                                           "RenormalizationScheme" )->second );
        LHPC::ParsingUtilities::TransformToUppercase( renormalizationScheme );
        if( renormalizationScheme == "MSBAR" )
        {
          vectorMassCorrectionConstant = ( 5.0 / 6.0 );
        }
        else if( renormalizationScheme == "DRBAR" )
        {
          vectorMassCorrectionConstant = 1.5;
        }
        else
        {
          std::stringstream errorBuilder;
          errorBuilder << "RenormalizationScheme was not MSBAR or DRBAR"
          << " (nothing else is currently supported)!";
          throw std::runtime_error( errorBuilder.str() );
        }
        xmlLoopCorrections = xmlParser.CurrentBody();
      }
    }
    if( xmlFieldVariables.empty() )
    {
      throw std::runtime_error( "Could not parse <FieldVariables>." );
    }
    if( xmlDsbMinimum.empty() )
    {
      throw std::runtime_error( "Could not parse <DsbMinimum>." );
    }
    if( xmlTreeLevelPotential.empty() )
    {
      throw std::runtime_error( "Could not parse <TreeLevelPotential>." );
    }
    if( xmlLoopCorrections.empty() )
    {
      throw std::runtime_error( "Could not parse <LoopCorrections>." );
    }

    // Now we parse <FieldVariables>.
    std::vector< std::string >
    fieldLines( LHPC::ParsingUtilities::SplitBySubstrings( xmlFieldVariables,
                                                           "\n" ) );
    std::string readFieldName( "" );
    for( std::vector< std::string >::const_iterator
         fieldLine( fieldLines.begin() );
         fieldLine != fieldLines.end();
         ++fieldLine )
    {
      readFieldName.assign(
        LHPC::ParsingUtilities::TrimWhitespaceFromFrontAndBack( *fieldLine ) );
      if( !(readFieldName.empty()) )
      {
        if( ( readFieldName.size() > positiveByConvention.size() )
            &&
            ( readFieldName.compare( ( readFieldName.size()
                                       - positiveByConvention.size() ),
                                     positiveByConvention.size(),
                                     positiveByConvention ) == 0 ) )
        {
          fieldsAssumedPositive.push_back( fieldNames.size() );
          fieldNames.push_back( lagrangianParameterManager.FormatVariable(
                       LHPC::ParsingUtilities::TrimWhitespaceFromFrontAndBack(
                                                       readFieldName.substr( 0,
                ( readFieldName.size() - positiveByConvention.size() ) ) ) ) );
        }
        else if( ( readFieldName.size() > negativeByConvention.size() )
                 &&
                 ( readFieldName.compare( ( readFieldName.size()
                                            - negativeByConvention.size() ),
                                          negativeByConvention.size(),
                                          negativeByConvention ) == 0 ) )
        {
          fieldsAssumedNegative.push_back( fieldNames.size() );
          fieldNames.push_back( lagrangianParameterManager.FormatVariable(
                        LHPC::ParsingUtilities::TrimWhitespaceFromFrontAndBack(
                                                       readFieldName.substr( 0,
                ( readFieldName.size() - negativeByConvention.size() ) ) ) ) );
        }
        else
        {
          fieldNames.push_back( lagrangianParameterManager.FormatVariable(
                                                             readFieldName ) );
        }
      }
    }
    numberOfFields = fieldNames.size();
    dsbFieldValueInputs.resize( numberOfFields );
    dsbFieldInputStrings.resize( numberOfFields );
    // </FieldVariables>
    // Now we parse <DsbMinimum>
    std::vector< std::string >
    dsbLines( LHPC::ParsingUtilities::SplitBySubstrings( xmlDsbMinimum,
                                                         "\n" ) );
    std::string trimmedLine( "" );
    for( std::vector< std::string >::const_iterator
         dsbLine( dsbLines.begin() );
         dsbLine != dsbLines.end();
         ++dsbLine )
    {
      trimmedLine.assign(
          LHPC::ParsingUtilities::TrimWhitespaceFromFrontAndBack( *dsbLine ) );
      if( trimmedLine.empty() )
      {
        continue;
      }
      size_t const equalsPosition( trimmedLine.find( '=' ) );
      if( !( equalsPosition < ( trimmedLine.size() - 1 ) ) )
      {
        std::stringstream errorBuilder;
        errorBuilder
        << "Field given no value in <DsbMinimum>! (Offending line = \""
        << *dsbLine << "\")";
        throw std::runtime_error( errorBuilder.str() );
      }

      readFieldName.assign(
                        LHPC::ParsingUtilities::TrimWhitespaceFromFrontAndBack(
                                                         trimmedLine.substr( 0,
                                                          equalsPosition ) ) );
      size_t const
      fieldIndex( FieldIndex( lagrangianParameterManager.FormatVariable(
                                                           readFieldName ) ) );
      if( !( fieldIndex < fieldNames.size() ) )
      {
        std::stringstream errorBuilder;
        errorBuilder
        << "Unknown field (\"" << trimmedLine.substr( 0,
                                                      equalsPosition )
        << "\") given value in DsbMinimum!";
        throw std::runtime_error( errorBuilder.str() );
      }
      dsbFieldInputStrings[ fieldIndex ]
      = lagrangianParameterManager.FormatVariable(
                                         LHPC::ParsingUtilities::TrimFromFront(
                                      trimmedLine.substr( equalsPosition + 1 ),
                                 LHPC::ParsingUtilities::WhitespaceChars() ) );
    }
    // </DsbMinimum>
    // Now we parse <TreeLevelPotential>
    ParseSumOfPolynomialTerms( xmlTreeLevelPotential,
                               treeLevelPotential );
    // </TreeLevelPotential>
    // <LoopCorrections>
    xmlParser.LoadString( xmlLoopCorrections );
    std::vector< std::string > matrixLines;
    while( xmlParser.ReadNextElement() )
    {
      //   <ExtraPolynomialPart>
      if( xmlParser.CurrentName() == "ExtraPolynomialPart" )
      {
        ParseSumOfPolynomialTerms( xmlParser.CurrentBody(),
                                   polynomialLoopCorrections );
      }
      //   </ExtraPolynomialPart>
      //   <RealBosonMassSquaredMatrix>
      else if( xmlParser.CurrentName() == "RealBosonMassSquaredMatrix" )
      {
        size_t const numberOfRows( PrepareMatrixLines( xmlParser.CurrentBody(),
                                                       matrixLines,
                                              "RealBosonMassSquaredMatrix" ) );
        RealMassesSquaredMatrix massSquaredMatrix( numberOfRows,
                                               xmlParser.CurrentAttributes() );
        for( size_t lineIndex( 0 );
             lineIndex < matrixLines.size();
             ++lineIndex )
        {
          ParseSumOfPolynomialTerms(
                        LHPC::ParsingUtilities::TrimWhitespaceFromFrontAndBack(
                                                    matrixLines[ lineIndex ] ),
                                     massSquaredMatrix.ElementAt( lineIndex ),
                                     false );
        }
        if( massSquaredMatrix.GetSpinType()
            == MassesSquaredCalculator::gaugeBoson )
        {
          vectorMassSquaredMatrices.push_back( massSquaredMatrix );
        }
        else
        {
          scalarMassSquaredMatrices.push_back( massSquaredMatrix );
        }
      }
      //   </RealBosonMassSquaredMatrix>
      //   <WeylFermionMassMatrix>
      else if( xmlParser.CurrentName() == "WeylFermionMassMatrix" )
      {
        size_t const numberOfRows( PrepareMatrixLines( xmlParser.CurrentBody(),
                                                       matrixLines,
                                                   "WeylFermionMassMatrix" ) );
        SymmetricComplexMassMatrix fermionMassMatrix( numberOfRows,
                                               xmlParser.CurrentAttributes() );
        for( size_t lineIndex( 0 );
             lineIndex < matrixLines.size();
             ++lineIndex )
        {
          ParseSumOfPolynomialTerms(
                        LHPC::ParsingUtilities::TrimWhitespaceFromFrontAndBack(
                                                    matrixLines[ lineIndex ] ),
                                    fermionMassMatrix.ElementAt( lineIndex ) );
        }
        fermionMassMatrices.push_back( fermionMassMatrix );
      }
      //   </WeylFermionMassMatrix>
      //   <ComplexWeylFermionMassSquaredMatrix>
      else if( xmlParser.CurrentName()
               == "ComplexWeylFermionMassSquaredMatrix" )
      {
        size_t const numberOfRows( PrepareMatrixLines( xmlParser.CurrentBody(),
                                                       matrixLines,
                                     "ComplexWeylFermionMassSquaredMatrix" ) );
        ComplexMassSquaredMatrix fermionMassSquaredMatrix( numberOfRows,
                                               xmlParser.CurrentAttributes() );
        for( size_t lineIndex( 0 );
             lineIndex < matrixLines.size();
             ++lineIndex )
        {
          ParseSumOfPolynomialTerms(
                        LHPC::ParsingUtilities::TrimWhitespaceFromFrontAndBack(
                                                    matrixLines[ lineIndex ] ),
                             fermionMassSquaredMatrix.ElementAt( lineIndex ) );
        }
        fermionMassSquaredMatrices.push_back( fermionMassSquaredMatrix );
      }
      //   </WeylFermionMassMatrix>
    }
    // </LoopCorrections>

    // Now we can fill the MassesSquaredCalculator* vectors, as their pointers
    // should remain valid as the other vectors do not change size any more
    // after the constructor.
    for( size_t pointerIndex( 0 );
         pointerIndex < scalarMassSquaredMatrices.size();
         ++pointerIndex )
    {
      scalarSquareMasses.push_back(
                                &(scalarMassSquaredMatrices[ pointerIndex ]) );
    }
    for( size_t pointerIndex( 0 );
         pointerIndex < fermionMassMatrices.size();
         ++pointerIndex )
    {
      fermionSquareMasses.push_back( &(fermionMassMatrices[ pointerIndex ]) );
    }
    for( size_t pointerIndex( 0 );
         pointerIndex < fermionMassSquaredMatrices.size();
         ++pointerIndex )
    {
      fermionSquareMasses.push_back(
                               &(fermionMassSquaredMatrices[ pointerIndex ]) );
    }
    for( size_t pointerIndex( 0 );
         pointerIndex < vectorMassSquaredMatrices.size();
         ++pointerIndex )
    {
      vectorSquareMasses.push_back(
                                &(vectorMassSquaredMatrices[ pointerIndex ]) );
    }
  }

  PotentialFromPolynomialWithMasses::~PotentialFromPolynomialWithMasses()
  {
    // This does nothing.
  }

  // This writes the potential as
  // def PotentialFunction( fv ): return ...
  // in pythonFilename for fv being an array of floating-point numbers in the
  // same order as they are for the field configurations as internal to this
  // C++ code. It uses the virtual function WriteActualPythonFunction.
  void PotentialFromPolynomialWithMasses::WriteAsPython(
                                      std::string const& pythonFilename ) const
  {
    std::ofstream pythonFile( pythonFilename.c_str() );
    pythonFile << std::setprecision( 12 );
    pythonFile << "# Automatically generated file! Modify at your own PERIL!\n"
    "# This file was created by Vevacious version "
    << VersionInformation::CurrentVersion() << "\n"
    "from __future__ import division\n"
    "import math\n"
    "import numpy\n"
    "\n"
    "\n"
    "# The temperature has to be held as a global variable because\n"
    "# CosmoTransitions.fullTunneling does not pass in the temperature.\n"
    "temperatureLinear = 0.0\n"
    "temperatureSquare = 0.0\n"
    "temperatureFourth = 0.0\n"
    "temperatureInverseSquare = -1.0\n"
    "def SetGlobalTemperature( globalTemperature ):\n"
    "  temperatureLinear = globalTemperature\n"
    "  temperatureSquare = ( temperatureLinear * temperatureLinear )\n"
    "  temperatureFourth = ( temperatureSquare * temperatureSquare )\n"
    "  if ( globalTemperature > 0.0 ):\n"
    "    temperatureInverseSquare = ( 1.0 / temperatureSquare )\n"
    "  else:\n"
    "    temperatureInverseSquare = -1.0\n"
    "\n"
    "\n"
    << ThermalFunctions::JFunctionsAsPython() << "\n"
    "\n"
    << lagrangianParameterManager.ParametersAsPython() << "\n"
    "# The Lagrangian parameters evaluated at the appropriate scale from the\n"
    "# C++ parameter manager are held in the following global array.\n"
    "fixedScaleLagrangianParameters = LagrangianParameters( math.log( "
    << lagrangianParameterManager.AppropriateSingleFixedScale()
    << " ) )\n"
    "\n"
    "# In the following functions, the field configuration is given by fv\n"
    "# and the Lagrangian parameters are given by lp.\n"
    "def TreeLevelContribution( fv,\n"
    "                           lp ):\n"
    "    return ( " << treeLevelPotential.AsPython() << " )\n"
    "\n"
    "# The field configuration is given by fv and the Lagrangian parameters\n"
    "# are given by lp.\n"
    "def PolynomialLoopCorrections( fv,\n"
    "                               lp ):\n"
    "    return ( " << polynomialLoopCorrections.AsPython() << " )\n"
    "\n"
    "def JustLoopCorrection( massesSquaredWithFactors,\n"
    "                        subtractionConstant,\n"
    "                        invQSq ):\n"
    "    cumulativeQuantumCorrection = 0.0\n"
    "    for massSquaredWithFactor in massesSquaredWithFactors:\n"
    "        currentQuantumCorrection = 0.0\n"
    "        for MSq in massSquaredWithFactor[ 0 ]:\n"
    "            MM = abs( MSq )\n"
    "            if ( MM > 0.0 ):\n"
    "                currentQuantumCorrection += ( MM**2\n"
    "                                       * ( math.log( MM * invQSq )\n"
    "                                           - subtractionConstant ) )\n"
    "        cumulativeQuantumCorrection += ( massSquaredWithFactor[ 1 ]\n"
    "                                         * currentQuantumCorrection )\n"
    "    return cumulativeQuantumCorrection\n"
    "\n"
    "def LoopAndThermalCorrection( massesSquaredWithFactors,\n"
    "                              subtractionConstant,\n"
    "                              invQSq,\n"
    "                              invTSq,\n"
    "                              JFunction ):\n"
    "    cumulativeQuantumCorrection = 0.0\n"
    "    cumulativeThermalCorrection = 0.0\n"
    "    for massSquaredWithFactor in massesSquaredWithFactors:\n"
    "        currentQuantumCorrection = 0.0\n"
    "        currentThermalCorrection = 0.0\n"
    "        for MSq in massSquaredWithFactor[ 0 ]:\n"
    "            MM = abs( MSq )\n"
    "            if ( MM > 1.0 ):\n"
    "                currentQuantumCorrection += ( MM**2\n"
    "                                       * ( math.log( MM * invQSq )\n"
    "                                           - subtractionConstant ) )\n"
    "            currentThermalCorrection += JFunction( MM * invTSq )\n"
    "        cumulativeQuantumCorrection += ( massSquaredWithFactor[ 1 ]\n"
    "                                         * currentQuantumCorrection )\n"
    "        cumulativeThermalCorrection += ( massSquaredWithFactor[ 1 ]\n"
    "                                         * currentThermalCorrection )\n"
    "    return [ cumulativeQuantumCorrection, cumulativeThermalCorrection ]\n"
    "\n"
    "def ScalarMassesSquaredWithFactors( fv,\n"
    "                                    lp ):\n"
    "    massesSquaredWithFactors = []\n";
    for( std::vector< RealMassesSquaredMatrix >::const_iterator
         whichMatrix( scalarMassSquaredMatrices.begin() );
         whichMatrix < scalarMassSquaredMatrices.end();
         ++whichMatrix )
    {
      pythonFile << "    massSquaredMatrix = numpy.array( [ ";
      for( std::vector< ParametersAndFieldsProductSum >::const_iterator
           matrixElement( whichMatrix->MatrixElements().begin() );
           matrixElement < whichMatrix->MatrixElements().end();
           ++matrixElement )
      {
        if( matrixElement != whichMatrix->MatrixElements().begin() )
        {
          pythonFile << ",\n";
        }
        pythonFile << matrixElement->AsPython();
      }
      pythonFile << " ] )\n"
      "    M = numpy.reshape( massSquaredMatrix, ( "
                                 << whichMatrix->NumberOfRows() << ", -1 ) )\n"
      "    massesSquaredWithFactors.append( [ numpy.linalg.eigvalsh( M ),\n"
      "                     " << whichMatrix->MultiplicityFactor() << " ] )\n";
    }
    pythonFile << "    return massesSquaredWithFactors\n"
    "\n"
    "def FermionMassesSquaredWithFactors( fv,\n"
    "                                     lp ):\n"
    "    massesSquaredWithFactors = []\n";
    for( std::vector< ComplexMassSquaredMatrix >::const_iterator
         whichMatrix( fermionMassSquaredMatrices.begin() );
         whichMatrix < fermionMassSquaredMatrices.end();
         ++whichMatrix )
    {
      pythonFile << "    massSquaredMatrix = numpy.array( [ ";
      for( std::vector< ComplexParametersAndFieldsProductSum >::const_iterator
           matrixElement( whichMatrix->MatrixElements().begin() );
           matrixElement < whichMatrix->MatrixElements().end();
           ++matrixElement )
      {
        if( matrixElement != whichMatrix->MatrixElements().begin() )
        {
          pythonFile << ",\n";
        }
        pythonFile << "( " << matrixElement->first.AsPython()
        << " + ( 1j * " << matrixElement->second.AsPython() << " ) )";
      }
      pythonFile << " ] )\n"
      "    M = numpy.reshape( massSquaredMatrix, ( "
                                 << whichMatrix->NumberOfRows() << ", -1 ) )\n"
      "    massesSquaredWithFactors.append( [ numpy.linalg.eigvalsh( M ),\n"
      "                     " << whichMatrix->MultiplicityFactor() << " ] )\n";
    }
    for( std::vector< SymmetricComplexMassMatrix >::const_iterator
         whichMatrix( fermionMassMatrices.begin() );
         whichMatrix < fermionMassMatrices.end();
         ++whichMatrix )
    {
      pythonFile << "    massMatrix = numpy.array( [ ";
      for( std::vector< ComplexParametersAndFieldsProductSum >::const_iterator
           matrixElement( whichMatrix->MatrixElements().begin() );
           matrixElement < whichMatrix->MatrixElements().end();
           ++matrixElement )
      {
        if( matrixElement != whichMatrix->MatrixElements().begin() )
        {
          pythonFile << ",\n";
        }
        pythonFile << "( " << matrixElement->first.AsPython()
        << " + ( 1j * " << matrixElement->second.AsPython() << " ) )";
      }
      pythonFile << " ] )\n"
      "    M = numpy.reshape( massMatrix, ( " << whichMatrix->NumberOfRows()
                                                                << ", -1 ) )\n"
      "    MM = M.dot( numpy.transpose( numpy.conjugate( M ) ) )\n"
      "    massesSquaredWithFactors.append( [ numpy.linalg.eigvalsh( MM ),\n"
                              << whichMatrix->MultiplicityFactor() << " ] )\n";
    }
    pythonFile << "    return massesSquaredWithFactors\n"
    "\n"
    "def VectorMassesSquaredWithFactors( fv,\n"
    "                                    lp ):\n"
    "    massesSquaredWithFactors = []\n";
    for( std::vector< RealMassesSquaredMatrix >::const_iterator
         whichMatrix( vectorMassSquaredMatrices.begin() );
         whichMatrix < vectorMassSquaredMatrices.end();
         ++whichMatrix )
    {
      pythonFile << "    massSquaredMatrix = numpy.array( [ ";
      for( std::vector< ParametersAndFieldsProductSum >::const_iterator
           matrixElement( whichMatrix->MatrixElements().begin() );
           matrixElement < whichMatrix->MatrixElements().end();
           ++matrixElement )
      {
        if( matrixElement != whichMatrix->MatrixElements().begin() )
        {
          pythonFile << ",\n";
        }
        pythonFile << matrixElement->AsPython();
      }
      pythonFile << " ] )\n"
      "    M = numpy.reshape( massSquaredMatrix, ( "
                                 << whichMatrix->NumberOfRows() << ", -1 ) )\n"
      "    massesSquaredWithFactors.append( [ numpy.linalg.eigvalsh( M ),\n"
      "                     " << whichMatrix->MultiplicityFactor() << " ] )\n";
    }
    pythonFile << "    return massesSquaredWithFactors\n"
    "\n"
    "\n"
    "loopFactor = ( 1.0 / ( 64.0 * math.pi * math.pi ) )\n"
    "thermalFactor = ( 1.0 / ( 2.0 * math.pi * math.pi ) )\n"
    "\n"
    "def JustLoopCorrections( fv,\n"
    "                         lp,\n"
    "                         invQSq ):\n"
    "    totalQuantumCorrections = 0.0\n"
    "    currentCorrection = JustLoopCorrection(\n"
    "                                    ScalarMassesSquaredWithFactors( fv,\n"
    "                                                                  lp ),\n"
    "                                            1.5,\n"
    "                                            invQSq )\n"
    "    totalQuantumCorrections += currentCorrection\n"
    "    currentCorrection = JustLoopCorrection(\n"
    "                                   FermionMassesSquaredWithFactors( fv,\n"
    "                                                                  lp ),\n"
    "                                            1.5,\n"
    "                                            invQSq )\n"
    "    totalQuantumCorrections -= ( 2.0 * currentCorrection )\n"
    "    currentCorrection = JustLoopCorrection(\n"
    "                                    VectorMassesSquaredWithFactors( fv,\n"
    "                                                                  lp ),\n"
    "                                " << vectorMassCorrectionConstant << ",\n"
    "                                            invQSq )\n"
    "    totalQuantumCorrections += ( 3.0 * currentCorrection )\n"
    "    return ( totalQuantumCorrections * loopFactor )\n"
    "\n"
    "def LoopAndThermalCorrections( fv,\n"
    "                               lp,\n"
    "                               invQSq,\n"
    "                               invTSq ):\n"
    "    totalQuantumCorrections = 0.0\n"
    "    totalThermalCorrections = 0.0\n"
    "    currentCorrections = LoopAndThermalCorrection(\n"
    "                                    ScalarMassesSquaredWithFactors( fv,\n"
    "                                                                  lp ),\n"
    "                                                   1.5,\n"
    "                                                   invQSq,\n"
    "                                                   invTSq,\n"
    "                                                   BosonicJ )\n"
    "# Real scalar degrees of freedom add to both quantum and thermal\n"
    "# corrections without any factors.\n"
    "    totalQuantumCorrections += currentCorrections[ 0 ]\n"
    "    totalThermalCorrections += currentCorrections[ 1 ]\n"
    "    currentCorrections = LoopAndThermalCorrection(\n"
    "                                   FermionMassesSquaredWithFactors( fv,\n"
    "                                                                  lp ),\n"
    "                                                   1.5,\n"
    "                                                   invQSq,\n"
    "                                                   invTSq,\n"
    "                                                   FermionicJ )\n"
    "# Weyl fermion degrees of freedom add to both quantum and thermal\n"
    "# corrections with a factor of 2, though there is an additional minus\n"
    "# sign for the quantum corrections. (The conventions used in Vevacious\n"
    "# are that the thermal correction functions already have any minus\n"
    "# signs for fermions.)\n"
    "    totalQuantumCorrections -= ( 2.0 * currentCorrections[ 0 ] )\n"
    "    totalThermalCorrections += ( 2.0 * currentCorrections[ 1 ] )\n"
    "    currentCorrections = LoopAndThermalCorrection(\n"
    "                                    VectorMassesSquaredWithFactors( fv,\n"
    "                                                                  lp ),\n"
    "                                " << vectorMassCorrectionConstant << ",\n"
    "                                                   invQSq,\n"
    "                                                   invTSq,\n"
    "                                                   BosonicJ )\n"
    "# Vector boson degrees of freedom add to quantum corrections with a\n"
    "# factor of 3 in dimensional regularization schemes, and to thermal\n"
    "# corrections with a factor of 2.\n"
    "    totalQuantumCorrections += ( 3.0 * currentCorrections[ 0 ] )\n"
    "    totalThermalCorrections += ( 2.0 * currentCorrections[ 1 ] )\n"
    "    return ( ( totalQuantumCorrections * loopFactor )\n"
    "             + ( ( totalThermalCorrections * thermalFactor )\n"
    "                 / ( invTSq**2 ) ) )\n"
    "\n"
    << WriteActualPythonFunction() << "\n"
    "\n"
    "# This can be changed to use just the tree-level potential, or the\n"
    "# loop-corrected potential without thermal corrections.\n"
    "UnderlyingPotential = LoopAndThermallyCorrectedPotential\n"
    "\n"
    "stepSize = 1.0\n"
    "\n"
    "def UnderlyingGradient( fv ):\n"
    "    potentialAtPoint = UnderlyingPotential( fv )\n"
    "    gradientArray = numpy.zeros( len( fv ) )\n"
    "    for whichField in range( len( fv ) ):\n"
    "        displacedPoint = fv.copy()\n"
    "        displacedPoint[ whichField ] += stepSize\n"
    "        gradientArray[ whichField ] = ( ( UnderlyingPotential(\n"
    "                                                       displacedPoint )\n"
    "                                          - potentialAtPoint )\n"
    "                                        / stepSize )\n"
    "    return gradientArray\n"
    "\n"
    "def PotentialForCosmotransitions( arrayOfArrays ):\n"
    "    if ( arrayOfArrays.shape == ( " << numberOfFields << ", ) ):\n"
    "        return UnderlyingPotential( arrayOfArrays )\n"
    "    elif ( arrayOfArrays.shape == ( len( arrayOfArrays ), "
                                                 << numberOfFields << " ) ):\n"
    "        returnArray = numpy.zeros( len( arrayOfArrays ) )\n"
    "        for whichIndex in range( len( arrayOfArrays ) ):\n"
    "            returnArray[ whichIndex ] = UnderlyingPotential(\n"
    "                                          arrayOfArrays[ whichIndex ] )\n"
    "        return returnArray\n"
    "    else:\n"
    "        return None\n"
    "\n"
    "def GradientForCosmotransitions( arrayOfArrays ):\n"
    "    if ( arrayOfArrays.shape == ( " << numberOfFields << ", ) ):\n"
    "        return UnderlyingGradient( arrayOfArrays )\n"
    "    elif ( arrayOfArrays.shape == ( len( arrayOfArrays ), "
                                                 << numberOfFields << " ) ):\n"
    "        returnMatrix = arrayOfArrays.copy()\n"
    "        for whichIndex in range( len( arrayOfArrays ) ):\n"
    "            returnMatrix[ whichIndex ] = UnderlyingGradient(\n"
    "                                          arrayOfArrays[ whichIndex ] )\n"
    "        return returnMatrix\n"
    "    else:\n"
    "        return None\n"
    "\n"
    "# End of automatically generated file!\n";
    pythonFile.close();
  }

  // This is just for derived classes.
  PotentialFromPolynomialWithMasses::PotentialFromPolynomialWithMasses(
                     LagrangianParameterManager& lagrangianParameterManager ) :
    PotentialFunction( lagrangianParameterManager ),
    treeLevelPotential(),
    polynomialLoopCorrections(),
    scalarSquareMasses(),
    fermionSquareMasses(),
    vectorSquareMasses(),
    scalarMassSquaredMatrices(),
    fermionMassMatrices(),
    fermionMassSquaredMatrices(),
    vectorMassSquaredMatrices(),
    vectorMassCorrectionConstant( 5.0 / 6.0 ),
    fieldsAssumedPositive(),
    fieldsAssumedNegative(),
    assumedPositiveOrNegativeTolerance( -1.0 )
  {
    // This protected constructor is just an initialization list only used by
    // derived classes which are going to fill up the data members in their own
    // constructors.
  }

  // This is just for derived classes.
  PotentialFromPolynomialWithMasses::PotentialFromPolynomialWithMasses(
                        PotentialFromPolynomialWithMasses const& copySource ) :
    PotentialFunction( copySource ),
    treeLevelPotential( copySource.treeLevelPotential ),
    polynomialLoopCorrections( copySource.polynomialLoopCorrections ),
    scalarSquareMasses(),
    fermionSquareMasses(),
    vectorSquareMasses(),
    scalarMassSquaredMatrices( copySource.scalarMassSquaredMatrices ),
    fermionMassMatrices( copySource.fermionMassMatrices ),
    fermionMassSquaredMatrices( copySource.fermionMassSquaredMatrices ),
    vectorMassSquaredMatrices( copySource.vectorMassSquaredMatrices ),
    vectorMassCorrectionConstant( copySource.vectorMassCorrectionConstant ),
    fieldsAssumedPositive( copySource.fieldsAssumedPositive ),
    fieldsAssumedNegative( copySource.fieldsAssumedNegative ),
    assumedPositiveOrNegativeTolerance(
                                copySource.assumedPositiveOrNegativeTolerance )
  {
    // Now we can fill the MassesSquaredCalculator* vectors, as their pointers
    // should remain valid as the other vectors do not change size any more
    // after the constructor.
    for( size_t pointerIndex( 0 );
         pointerIndex < scalarMassSquaredMatrices.size();
         ++pointerIndex )
    {
      scalarSquareMasses.push_back(
                                &(scalarMassSquaredMatrices[ pointerIndex ]) );
    }
    for( size_t pointerIndex( 0 );
         pointerIndex < fermionMassMatrices.size();
         ++pointerIndex )
    {
      fermionSquareMasses.push_back( &(fermionMassMatrices[ pointerIndex ]) );
    }
    for( size_t pointerIndex( 0 );
         pointerIndex < fermionMassSquaredMatrices.size();
         ++pointerIndex )
    {
      fermionSquareMasses.push_back(
                               &(fermionMassSquaredMatrices[ pointerIndex ]) );
    }
    for( size_t pointerIndex( 0 );
         pointerIndex < vectorMassSquaredMatrices.size();
         ++pointerIndex )
    {
      vectorSquareMasses.push_back(
                                &(vectorMassSquaredMatrices[ pointerIndex ]) );
    }
  }


  // This evaluates the one-loop potential with thermal corrections assuming
  // that the squared masses were evaluated at the given scale correctly.
  double
  PotentialFromPolynomialWithMasses::LoopAndThermalCorrections(
          std::vector< DoubleVectorWithDouble > scalarMassesSquaredWithFactors,
         std::vector< DoubleVectorWithDouble > fermionMassesSquaredWithFactors,
          std::vector< DoubleVectorWithDouble > vectorMassesSquaredWithFactors,
                                              double const inverseScaleSquared,
                                          double const temperatureValue ) const
  {
    bool const temperatureGreaterThanZero( temperatureValue > 0.0 );
    double const inverseTemperatureSquared( temperatureGreaterThanZero ?
                            ( 1.0 / ( temperatureValue * temperatureValue ) ) :
                                            -1.0 );
    double totalQuantumCorrections( 0.0 );
    double totalThermalCorrections( 0.0 );

    double scalarQuantumCorrections( 0.0 );
    double scalarThermalCorrections( 0.0 );
    AddToCorrections( scalarMassesSquaredWithFactors,
                      inverseScaleSquared,
                      temperatureGreaterThanZero,
                      inverseTemperatureSquared,
                      1.5,
                      &(ThermalFunctions::BosonicJ),
                      scalarQuantumCorrections,
                      scalarThermalCorrections );

    // Real scalar degrees of freedom add to both quantum and thermal
    // corrections without any factors.
    totalQuantumCorrections += scalarQuantumCorrections;
    totalThermalCorrections += scalarThermalCorrections;

    double fermionQuantumCorrections( 0.0 );
    double fermionThermalCorrections( 0.0 );
    AddToCorrections( fermionMassesSquaredWithFactors,
                      inverseScaleSquared,
                      temperatureGreaterThanZero,
                      inverseTemperatureSquared,
                      1.5,
                      &(ThermalFunctions::FermionicJ),
                      fermionQuantumCorrections,
                      fermionThermalCorrections );

    // Weyl fermion degrees of freedom add to both quantum and thermal
    // corrections with a factor of 2, though there is an additional minus sign
    // for the quantum corrections. (The conventions used in Vevacious are that
    // the thermal correction functions already have any minus signs for
    // fermions.)
    totalQuantumCorrections -= ( 2.0 * fermionQuantumCorrections );
    totalThermalCorrections += ( 2.0 * fermionThermalCorrections );

    double vectorQuantumCorrections( 0.0 );
    double vectorThermalCorrections( 0.0 );
    AddToCorrections( vectorMassesSquaredWithFactors,
                      inverseScaleSquared,
                      temperatureGreaterThanZero,
                      inverseTemperatureSquared,
                      vectorMassCorrectionConstant,
                      &(ThermalFunctions::BosonicJ),
                      vectorQuantumCorrections,
                      vectorThermalCorrections );

    // Vector boson degrees of freedom add to quantum corrections with a factor
    // of 3 in dimensional regularization schemes, and to thermal corrections
    // with a factor of 2.
    totalQuantumCorrections += ( 3.0 * vectorQuantumCorrections );
    totalThermalCorrections += ( 2.0 * vectorThermalCorrections );

    return ( ( totalQuantumCorrections * loopFactor )
             + ( totalThermalCorrections * thermalFactor
                 * temperatureValue * temperatureValue
                 * temperatureValue * temperatureValue ) );
  }

  // This interprets stringToParse as a sum of polynomial terms and sets
  // polynomialSum accordingly.
  void PotentialFromPolynomialWithMasses::ParseSumOfPolynomialTerms(
                                              std::string const& stringToParse,
                          ComplexParametersAndFieldsProductSum& polynomialSum )
  {
    polynomialSum.first.ParametersAndFieldsProducts().clear();
    polynomialSum.second.ParametersAndFieldsProducts().clear();
    if( stringToParse.empty() )
    {
      return;
    }

    bool positiveTerm( true );

    // We know that stringToParse[ 0 ] is not whitespace.
    size_t wordStart( 0 );
    if( stringToParse[ 0 ] == '-' )
    {
      positiveTerm = false;
      wordStart = 1;
    }
    else if( stringToParse[ 0 ] == '+' )
    {
      wordStart = 1;
    }

    // Now we have skipped any initial '+' or '-', so we start the 1st term.
    ParametersAndFieldsProductTerm polynomialTerm;
    if( !positiveTerm )
    {
      polynomialTerm.MultiplyByConstant( -1.0 );
    }

    bool imaginaryTerm( false );
    wordStart = PutNextNumberOrVariableIntoPolynomial( stringToParse,
                                                       wordStart,
                                                       polynomialTerm,
                                                       imaginaryTerm );

    // After parsing the 1st word, we keep parsing words until the end of
    // stringToParse is reached.
    while( wordStart != std::string::npos )
    {
      // If the last word ended with '+' or '-' (internal '+'/'-' within
      // floating point numbers and within [...] are considered part of a word
      // and do not mark the end of the word), then we are about to move on to
      // a new polynomial term.
      if( ( stringToParse[ wordStart ] == '+' )
          ||
          ( stringToParse[ wordStart ] == '-' ) )
      {
        if( polynomialTerm.IsValid() )
        {
          if( imaginaryTerm )
          {
            polynomialSum.second.ParametersAndFieldsProducts().push_back(
                                                              polynomialTerm );
          }
          else
          {
            polynomialSum.first.ParametersAndFieldsProducts().push_back(
                                                              polynomialTerm );
          }
        }
        polynomialTerm.ResetValues();
        if( stringToParse[ wordStart ] == '-' )
        {
          polynomialTerm.MultiplyByConstant( -1.0 );
        }
      }

      // Now we parse the next word.
      wordStart = PutNextNumberOrVariableIntoPolynomial( stringToParse,
                                                         (++wordStart),
                                                         polynomialTerm,
                                                         imaginaryTerm );
    }

    if( polynomialTerm.IsValid() )
    {
      if( imaginaryTerm )
      {
        polynomialSum.second.ParametersAndFieldsProducts().push_back(
                                                              polynomialTerm );
      }
      else
      {
        polynomialSum.first.ParametersAndFieldsProducts().push_back(
                                                              polynomialTerm );
      }
    }
  }

  // This reads in a whole number or variable (including possible raising to
  // a power), applies the correct operation to polynomialTerm, and then
  // returns the position of the character just after the interpreted word.
  // If there was a factor of "i", "I", "j", or "J", imaginaryTerm is set to
  // true.
  size_t
  PotentialFromPolynomialWithMasses::PutNextNumberOrVariableIntoPolynomial(
                                              std::string const& stringToParse,
                                                              size_t wordStart,
                                ParametersAndFieldsProductTerm& polynomialTerm,
                                                          bool& imaginaryTerm )
  {
    size_t wordEnd( 0 );
    while( wordStart < stringToParse.size() )
    {
      // First we check to see if it's a number:
      if( dotAndDigits.find( stringToParse[ wordStart ] )
          != std::string::npos )
      {
        wordEnd = stringToParse.find_first_not_of( dotAndDigits,
                                                   wordStart );
        if( ( wordEnd != std::string::npos )
            &&
            ( ( stringToParse[ wordEnd ] == 'e' )
              ||
              ( stringToParse[ wordEnd ] == 'E' ) ) )
        {
          // If the number is in "e notation", we have to check that it is not
          // malformed:
          if( ( wordEnd < ( stringToParse.size() - 2 ) )
              &&
              ( digitChars.find( stringToParse[ wordEnd + 1 ] )
                != std::string::npos ) )
          {
            wordEnd = stringToParse.find_first_not_of( digitChars,
                                                       ( wordEnd + 1 ) );
          }
          else if( ( wordEnd < ( stringToParse.size() - 3 ) )
                   &&
                   ( ( stringToParse[ wordEnd + 1 ] == '+' )
                     ||
                     ( stringToParse[ wordEnd + 1 ] == '-' ) )
                   &&
                   ( digitChars.find( stringToParse[ wordEnd + 2 ] )
                     != std::string::npos ) )
          {
            wordEnd = stringToParse.find_first_not_of( digitChars,
                                                       ( wordEnd + 2 ) );
          }
          else
          {
            throw std::runtime_error(
                           "Model file had malformed scientific E notation." );
          }
        }
        // If it is a number, we multiply the polynomial term and return the
        // position of the char just after the chars that we have just parsed.
        polynomialTerm.MultiplyByConstant(
                                        LHPC::ParsingUtilities::StringToDouble(
                                               stringToParse.substr( wordStart,
                                                 ( wordEnd - wordStart ) ) ) );
        return wordEnd;
      }
      else if( allowedVariableInitials.find( stringToParse[ wordStart ] )
               != std::string::npos )
      {
        // Otherwise we have to check to see if it is a string that could be a
        // field variable name or a string mapping to a Lagrangian parameter.
        wordEnd = stringToParse.find_first_not_of( allowedVariableChars,
                                                   wordStart );
        if( ( wordEnd < ( stringToParse.size() - 1 ) )
            &&
            ( stringToParse[ wordEnd ] == '[' ) )
        {
          wordEnd = stringToParse.find( ']',
                                        wordEnd );
          if( wordEnd == std::string::npos )
          {
            throw std::runtime_error( "Model file had unclosed [...]." );
          }
          // Since stringToParse[ wordEnd ] is ']', it is OK to increment
          // wordEnd by 1 so that the ']' is included in the substring taken
          // for parsing later.
          wordEnd += 1;
        }

        // For comparison, we need the string to be in the proper format.
        std::string variableString( lagrangianParameterManager.FormatVariable(
                                               stringToParse.substr( wordStart,
                                                 ( wordEnd - wordStart ) ) ) );
        unsigned int powerInt( 1 );
        if( ( wordEnd < stringToParse.size() )
            &&
            ( stringToParse[ wordEnd ] == '^' ) )
        {
          if( ( wordEnd == ( stringToParse.size() - 1 ) )
              ||
              ( digitChars.find( stringToParse[ wordEnd + 1 ] )
                == std::string::npos ) )
          {
            // Only positive integers are allowed as exponents. It's easier to
            // also demand that no '+' or whitespace characters are allowed.
            throw std::runtime_error(
                              "Model file had invalid exponent after \'^\'." );
          }
          wordStart = ( wordEnd + 1 );
          wordEnd = stringToParse.find_first_not_of( digitChars,
                                                     wordStart );
          powerInt = static_cast< unsigned int >(
                                    LHPC::ParsingUtilities::BaseTenStringToInt(
                                               stringToParse.substr( wordStart,
                                                 ( wordEnd - wordStart ) ) ) );
        }

        // First we check for the imaginary unit.
        if( ( variableString.size() == 1 )
            &&
            ( ( variableString[ 0 ] == 'i' )
              ||
              ( variableString[ 0 ] == 'I' )
              ||
              ( variableString[ 0 ] == 'j' )
              ||
              ( variableString[ 0 ] == 'J' ) ) )
        {
          // If the term was already imaginary, we take its pre-existing factor
          // of the imaginary unit into the power just read.
          if( imaginaryTerm )
          {
            ++powerInt;
          }
          imaginaryTerm = ( ( powerInt % 2 ) == 1 );
          if( ( powerInt % 4 ) > 1 )
          {
            polynomialTerm.MultiplyByConstant( -1.0 );
          }
          return wordEnd;
        }

        // Next we check for a field name.
        for( size_t fieldIndex( 0 );
             fieldIndex < fieldNames.size();
             ++fieldIndex )
        {
          if( fieldNames[ fieldIndex ] == variableString )
          {
            // If it is a field name, we raise its power in the polynomial term
            // and return position of the char just after the parsed
            // characters.
            polynomialTerm.RaiseFieldPower( fieldIndex,
                                            powerInt );
            return wordEnd;
          }
        }

        // If we get to here, it wasn't a field name, so we check whether it's
        // a Lagrangian parameter known to lagrangianParameterManager.
        std::pair< bool, size_t > parameterValidityAndIndex(
              lagrangianParameterManager.RegisterParameter( variableString ) );
        if( parameterValidityAndIndex.first )
        {
          polynomialTerm.MultiplyByParameter( parameterValidityAndIndex.second,
                                              powerInt );
        }
        else
        {
          polynomialTerm.MarkAsInvalid();
        }
        return wordEnd;
      }
      else if( ( stringToParse[ wordStart ] == '+' )
               ||
               ( stringToParse[ wordStart ] == '-' ) )
      {
        return wordStart;
      }
      else
      {
        wordStart++;
      }
    }
    return std::string::npos;
  }

  // This evaluates the sum of corrections for the degrees of freedom with
  // masses-squared given by massesSquaredWithFactors with
  // subtractFromLogarithm as the constant to subtract from the logarithm of
  // the ratio of mass-squared to square of renormalization scale, at a
  // temperature given by inverseTemperatureSquared^(-1/2) using a thermal
  // correction function given by ThermalFunction, and adds them to
  // cumulativeQuantumCorrection and cumulativeThermalCorrection.
  void PotentialFromPolynomialWithMasses::AddToCorrections(
         std::vector< DoubleVectorWithDouble > const& massesSquaredWithFactors,
                                              double const inverseScaleSquared,
                                         bool const temperatureGreaterThanZero,
                                        double const inverseTemperatureSquared,
                                            double const subtractFromLogarithm,
                                     double (*ThermalFunction)( double const ),
                                           double& cumulativeQuantumCorrection,
                                    double& cumulativeThermalCorrection ) const
  {
    double currentQuantumCorrection( 0.0 );
    double currentThermalCorrection( 0.0 );
    double massSquared( 0.0 );
    for( std::vector< DoubleVectorWithDouble >::const_iterator
         massesSquared( massesSquaredWithFactors.begin() );
         massesSquared < massesSquaredWithFactors.end();
         ++massesSquared )
    {
      currentQuantumCorrection = 0.0;
      currentThermalCorrection = 0.0;
      for( std::vector< double >::const_iterator
           massSquaredIterator( massesSquared->first.begin() );
           massSquaredIterator < massesSquared->first.end();
           ++massSquaredIterator )
      {
        massSquared = *massSquaredIterator;
        if( massSquared < 0.0 )
        {
          massSquared = -massSquared;
        }
        if( massSquared > 0.0 )
        {
          currentQuantumCorrection += ( massSquared * massSquared
                                        * ( log( massSquared
                                                 * inverseScaleSquared )
                                            - subtractFromLogarithm ) );
        }
        if( temperatureGreaterThanZero )
        {
          currentThermalCorrection += (*ThermalFunction)( massSquared
                                                 * inverseTemperatureSquared );
        }
      }
      cumulativeQuantumCorrection
      += ( massesSquared->second * currentQuantumCorrection );
      cumulativeThermalCorrection
      += ( massesSquared->second * currentThermalCorrection );
    }
  }

  // This is for debugging.
  std::string PotentialFromPolynomialWithMasses::AsDebuggingString() const
  {
    std::stringstream returnStream;
    returnStream << "fieldNames, dsbFieldInputStrings:" << std::endl;
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      returnStream << fieldNames[ fieldIndex ] << ": "
      << dsbFieldInputStrings[ fieldIndex ] << std::endl;
    }
    returnStream
    << std::endl
    << "treeLevelPotential = " << treeLevelPotential.AsDebuggingString()
    << std::endl
    << "polynomialLoopCorrections = "
    << polynomialLoopCorrections.AsDebuggingString()
    << std::endl
    << "scalarSquareMasses = " << std::endl;
    for( std::vector< RealMassesSquaredMatrix >::const_iterator
         whichScalar( scalarMassSquaredMatrices.begin() );
         whichScalar < scalarMassSquaredMatrices.end();
         ++whichScalar )
    {
      returnStream << whichScalar->AsString();
    }
    returnStream << std::endl;
    returnStream << std::endl
    << "fermionMasses = " << std::endl;
    for( std::vector< SymmetricComplexMassMatrix >::const_iterator
         whichFermion( fermionMassMatrices.begin() );
         whichFermion < fermionMassMatrices.end();
         ++whichFermion )
    {
      returnStream << whichFermion->AsString();
    }
    returnStream << std::endl;
    returnStream << std::endl
    << "fermionMassSquareds = " << std::endl;
    for( std::vector< ComplexMassSquaredMatrix >::const_iterator
         whichFermion( fermionMassSquaredMatrices.begin() );
         whichFermion < fermionMassSquaredMatrices.end();
         ++whichFermion )
    {
      returnStream << whichFermion->AsString();
    }
    returnStream << std::endl;
    returnStream << std::endl
    << "vectorSquareMasses = " << std::endl;
    for( std::vector< RealMassesSquaredMatrix >::const_iterator
         whichVector( vectorMassSquaredMatrices.begin() );
         whichVector < vectorMassSquaredMatrices.end();
         ++whichVector )
    {
      returnStream << whichVector->AsString();
    }
    returnStream << std::endl
    << "vectorMassCorrectionConstant = " << vectorMassCorrectionConstant
    << std::endl;
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
