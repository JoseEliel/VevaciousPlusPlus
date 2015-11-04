/*
 * PotentialFromPolynomialAndMasses.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/PotentialFunctions/OldPotentialFromPolynomialAndMasses.hpp"

#include "VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{
  std::string const
  OldPotentialFromPolynomialAndMasses::digitChars( "0123456789" );
  std::string const
  OldPotentialFromPolynomialAndMasses::dotAndDigits( "."
                           + OldPotentialFromPolynomialAndMasses::digitChars );
  std::string const
  OldPotentialFromPolynomialAndMasses::allowedVariableInitials(
                                                   "qwertyuiopasdfghjklzxcvbnm"
                                                "QWERTYUIOPASDFGHJKLZXCVBNM" );
  std::string const OldPotentialFromPolynomialAndMasses::allowedVariableChars(
                   OldPotentialFromPolynomialAndMasses::allowedVariableInitials
                              + OldPotentialFromPolynomialAndMasses::digitChars
                                                                      + "_~" );
  double const OldPotentialFromPolynomialAndMasses::piSquared(
       boost::math::double_constants::pi * boost::math::double_constants::pi );
  double const
  OldPotentialFromPolynomialAndMasses::loopFactor( 1.0
                    / ( 64.0 * OldPotentialFromPolynomialAndMasses::piSquared ) );
  double const OldPotentialFromPolynomialAndMasses::thermalFactor( 1.0
                     / ( 2.0 * OldPotentialFromPolynomialAndMasses::piSquared ) );
  std::string const OldPotentialFromPolynomialAndMasses::positiveByConvention(
                                                      "PositiveByConvention" );
  std::string const OldPotentialFromPolynomialAndMasses::negativeByConvention(
                                                      "NegativeByConvention" );

  OldPotentialFromPolynomialAndMasses::OldPotentialFromPolynomialAndMasses(
                                              std::string const& modelFilename,
                                          double const scaleRangeMinimumFactor,
            bool const treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions,
                               double const assumedPositiveOrNegativeTolerance,
                           RunningParameterManager& runningParameterManager ) :
    PotentialFunction(),
    IWritesPythonPotential(),
    ParameterUpdatePropagator( runningParameterManager ),
    runningParameters( runningParameterManager ),
    dsbFieldValuePolynomials(),
    currentMinimumRenormalizationScale( -1.0 ),
    squareOfMinimumRenormalizationScale( -1.0 ),
    currentMaximumRenormalizationScale( -1.0 ),
    squareOfMaximumRenormalizationScale( -1.0 ),
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
    treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions(
                     treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions ),
    scaleRangeMinimumFactor( scaleRangeMinimumFactor ),
    fieldsAssumedPositive(),
    fieldsAssumedNegative(),
    assumedPositiveOrNegativeTolerance( assumedPositiveOrNegativeTolerance )
  {
    BOL::AsciiXmlParser fileParser( false );
    BOL::AsciiXmlParser elementParser( false );
    BOL::VectorlikeArray< std::string > elementLines;
    bool successfullyReadElement(
                           fileParser.openRootElementOfFile( modelFilename ) );
    if( !successfullyReadElement )
    {
      throw std::runtime_error( "Could not parse XML of " + modelFilename );
    }
    // The model file should always have the XML elements in the correct order!
    // <ModelFileDetails>
    successfullyReadElement = fileParser.readNextElement();
    if( !successfullyReadElement )
    {
      throw std::runtime_error( "Could not parse <ModelFileDetails>." );
    }
    // </ModelFileDetails>
    // <AllowedNonZeroVariables>
    successfullyReadElement = fileParser.readNextElement();
    if( !successfullyReadElement )
    {
      throw std::runtime_error( "Could not parse <AllowedNonZeroVariables>." );
    }
    successfullyReadElement
    = elementParser.loadString( fileParser.getCurrentElementContent() );
    if( !successfullyReadElement )
    {
      throw std::runtime_error( "Could not parse <AllowedNonZeroVariables>." );
    }
    //   <FieldVariables>
    successfullyReadElement = elementParser.readNextElement();
    if( !successfullyReadElement )
    {
      throw std::runtime_error( "Could not parse <FieldVariables>." );
    }
    BOL::StringParser::parseByChar(
                               elementParser.getTrimmedCurrentElementContent(),
                                    elementLines,
                                    '\n');
    std::string readFieldName( "" );
    for( int lineIndex( 0 );
         lineIndex < elementLines.getSize();
         ++lineIndex )
    {
      readFieldName.assign( BOL::StringParser::trimFromFrontAndBack(
                                                 elementLines[ lineIndex ] ) );
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
          fieldNames.push_back( FormatVariable(
                                       BOL::StringParser::trimFromFrontAndBack(
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
          fieldNames.push_back( FormatVariable(
                                       BOL::StringParser::trimFromFrontAndBack(
                                                       readFieldName.substr( 0,
                ( readFieldName.size() - negativeByConvention.size() ) ) ) ) );
        }
        else
        {
          fieldNames.push_back( FormatVariable( readFieldName ) );
        }
      }
    }
    numberOfFields = fieldNames.size();
    dsbFieldValuePolynomials.resize( numberOfFields );
    dsbFieldValueInputs.resize( numberOfFields );
    //   </FieldVariables>
    //   <SlhaBlocks>
    successfullyReadElement = elementParser.readNextElement();
    if( !successfullyReadElement )
    {
      throw std::runtime_error( "Could not parse <SlhaBlocks>." );
    }
    std::string slhaString;
    std::string aliasString;
    elementLines.clearEntries();
    BOL::StringParser::parseByChar(
                               elementParser.getTrimmedCurrentElementContent(),
                                    elementLines,
                                    '\n');
    for( int lineIndex( 0 );
         lineIndex < elementLines.getSize();
         ++lineIndex )
    {
      slhaString = BOL::StringParser::trimFromFrontAndBack(
                     BOL::StringParser::firstWordOf( elementLines[ lineIndex ],
                                                     &aliasString,
                                                     "=" ) );
      if( !(slhaString.empty()) )
      {
        runningParameters.AddValidSlhaBlock( slhaString,
                      BOL::StringParser::trimFromFrontAndBack( aliasString ) );
      }
    }
    //   </SlhaBlocks>
    //   <DerivedParameters>
    successfullyReadElement = elementParser.readNextElement();
    if( !successfullyReadElement )
    {
      throw std::runtime_error( "Could not parse <DerivedParameters>." );
    }
    std::string readableName;
    elementLines.clearEntries();
    BOL::StringParser::parseByChar(
                               elementParser.getTrimmedCurrentElementContent(),
                                    elementLines,
                                    '\n');
    for( int lineIndex( 0 );
         lineIndex < elementLines.getSize();
         ++lineIndex )
    {
      aliasString
      = BOL::StringParser::trimFromFrontAndBack( elementLines[ lineIndex ] );
      if( !(aliasString.empty()) )
      {
        runningParameters.CreateDerivedParameter( aliasString );
      }
    }
    //   </DerivedParameters>
    //   <DsbMinimum>
    successfullyReadElement = elementParser.readNextElement();
    if( !successfullyReadElement )
    {
      throw std::runtime_error( "Could not parse <DsbMinimum>." );
    }
    elementLines.clearEntries();
    BOL::StringParser::parseByChar(
                               elementParser.getTrimmedCurrentElementContent(),
                                    elementLines,
                                    '\n');
    for( int lineIndex( 0 );
         lineIndex < elementLines.getSize();
         ++lineIndex )
    {
      size_t equalsPosition( elementLines[ lineIndex ].find( '=' ) );
      if( !( equalsPosition < ( elementLines[ lineIndex ].size() - 1 ) ) )
      {
        std::string errorMessage(
                 "Field given no value in DsbMinimum! (Offending line = \"" );
        errorMessage.append( elementLines[ lineIndex ] );
        errorMessage.append( "\")" );
        throw std::runtime_error( errorMessage );
      }
      size_t fieldIndex( FieldIndex( BOL::StringParser::trimFromFrontAndBack(
                                           elementLines[ lineIndex ].substr( 0,
                                                        equalsPosition ) ) ) );
      if( !( fieldIndex < fieldNames.size() ) )
      {
        std::string errorMessage( "Unknown field (\"" );
        errorMessage.append( elementLines[ lineIndex ].substr( 0,
                                                            equalsPosition ) );
        errorMessage.append( "\") given value in DsbMinimum!" );
        throw std::runtime_error( errorMessage );
      }
      ParseSumOfPolynomialTerms( BOL::StringParser::trimFromFrontAndBack(
                      elementLines[ lineIndex ].substr( equalsPosition + 1 ) ),
                                 dsbFieldValuePolynomials[ fieldIndex ] );
    }
    //   </DsbMinimum>
    // </AllowedNonZeroVariables>
    // <TreeLevelPotential>
    successfullyReadElement = fileParser.readNextElement();
    if( !successfullyReadElement )
    {
      throw std::runtime_error( "Could not parse <TreeLevelPotential>." );
    }
    ParseSumOfPolynomialTerms( fileParser.getTrimmedCurrentElementContent(),
                               treeLevelPotential );
    // </TreeLevelPotential>
    // <LoopCorrections>
    successfullyReadElement = fileParser.readNextElement();
    if( !successfullyReadElement )
    {
      throw std::runtime_error( "Could not parse <LoopCorrections>." );
    }
    std::string
    renormalizationScheme( fileParser.getCurrentElementAttributes().find(
                                           "RenormalizationScheme" )->second );
    if( renormalizationScheme.compare( "MSBAR" ) == 0 )
    {
      vectorMassCorrectionConstant = ( 5.0 / 6.0 );
    }
    else if( renormalizationScheme.compare( "DRBAR" ) == 0 )
    {
      vectorMassCorrectionConstant = 1.5;
    }
    else
    {
      throw std::runtime_error( "RenormalizationScheme was not MSBAR or DRBAR"
                                " (nothing else is currently supported)!" );
    }
    elementParser.loadString( fileParser.getCurrentElementContent() );
    while( elementParser.readNextElement() )
    {
      //   <ExtraPolynomialPart>
      if( elementParser.currentElementNameMatches( "ExtraPolynomialPart" ) )
      {
        ParseSumOfPolynomialTerms(
                              elementParser.getTrimmedCurrentElementContent(),
                                   polynomialLoopCorrections );
      }
      //   </ExtraPolynomialPart>
      //   <RealBosonMassSquaredMatrix>
      else if( elementParser.currentElementNameMatches(
                                               "RealBosonMassSquaredMatrix" ) )
      {
        elementLines.clearEntries();
        BOL::StringParser::parseByChar(
                               elementParser.getTrimmedCurrentElementContent(),
                                        elementLines,
                                        '\n');
        int numberOfRows(
                     sqrt( static_cast< double >( elementLines.getSize() ) ) );
        if( ( numberOfRows * numberOfRows ) != elementLines.getSize() )
        {
          throw std::runtime_error( "Number of elements for"
                     " RealBosonMassSquaredMatrix was not a square integer!" );
        }
        OldRealMassesSquaredMatrix massSquaredMatrix( numberOfRows,
                                 elementParser.getCurrentElementAttributes() );
        for( int lineIndex( 0 );
             lineIndex < elementLines.getSize();
             ++lineIndex )
        {
          ParseSumOfPolynomialTerms( BOL::StringParser::trimFromFrontAndBack(
                                                   elementLines[ lineIndex ] ),
                                    massSquaredMatrix.ElementAt( lineIndex ) );
        }
        if( massSquaredMatrix.GetSpinType()
            == OldMassesSquaredCalculator::gaugeBoson )
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
      else if( elementParser.currentElementNameMatches(
                                                    "WeylFermionMassMatrix" ) )
      {
        elementLines.clearEntries();
        BOL::StringParser::parseByChar(
                               elementParser.getTrimmedCurrentElementContent(),
                                        elementLines,
                                        '\n');
        int numberOfRows(
                     sqrt( static_cast< double >( elementLines.getSize() ) ) );
        if( ( numberOfRows * numberOfRows ) != elementLines.getSize() )
        {
          throw std::runtime_error( "Number of elements for"
                          " WeylFermionMassMatrix was not a square integer!" );
        }
        OldSymmetricComplexMassMatrix fermionMassMatrix( numberOfRows,
                                 elementParser.getCurrentElementAttributes() );
        for( int lineIndex( 0 );
             lineIndex < elementLines.getSize();
             ++lineIndex )
        {
          ParseSumOfPolynomialTerms( BOL::StringParser::trimFromFrontAndBack(
                                                   elementLines[ lineIndex ] ),
                                    fermionMassMatrix.ElementAt( lineIndex ) );
        }
        fermionMassMatrices.push_back( fermionMassMatrix );
      }
      //   </WeylFermionMassMatrix>
      //   <ComplexWeylFermionMassSquaredMatrix>
      else if( elementParser.currentElementNameMatches(
                                      "ComplexWeylFermionMassSquaredMatrix" ) )
      {
        elementLines.clearEntries();
        BOL::StringParser::parseByChar(
                               elementParser.getTrimmedCurrentElementContent(),
                                        elementLines,
                                        '\n');
        int numberOfRows(
                     sqrt( static_cast< double >( elementLines.getSize() ) ) );
        if( ( numberOfRows * numberOfRows ) != elementLines.getSize() )
        {
          throw std::runtime_error( "Number of elements for"
            " ComplexWeylFermionMassSquaredMatrix was not a square integer!" );
        }
        OldComplexMassSquaredMatrix fermionMassSquaredMatrix( numberOfRows,
                                 elementParser.getCurrentElementAttributes() );
        for( int lineIndex( 0 );
             lineIndex < elementLines.getSize();
             ++lineIndex )
        {
          ParseSumOfPolynomialTerms( BOL::StringParser::trimFromFrontAndBack(
                                                   elementLines[ lineIndex ] ),
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

  OldPotentialFromPolynomialAndMasses::~OldPotentialFromPolynomialAndMasses()
  {
    // This does nothing.
  }


  // This writes the potential as
  // def PotentialFunction( fv ): return ...
  // in pythonFilename for fv being an array of floating-point numbers in the
  // same order as they are for the field configurations as internal to this
  // C++ code. It uses the virtual function SetScaleForPythonPotentialCall.
  void OldPotentialFromPolynomialAndMasses::WriteAsPython(
                                       std::string const pythonFilename ) const
  {
    std::ofstream pythonFile( pythonFilename.c_str() );
    pythonFile << std::setprecision( 12 );
    pythonFile << "# Automatically generated file! Modify at your own PERIL!\n"
    "# This file was created by Vevacious version "
    << VersionInformation::currentVersion << "\n"
    "from __future__ import division\n"
    "import math\n"
    "import numpy\n"
    "\n"
    "# Global variables:\n"
    "# Temperature T given as T^(-2):\n"
    "invTSq = -1.0E+20\n"
    "# A very large negative number should be used for T=0, so that the code\n"
    "# to set the scale knows that T=0, and any thermal corrections that may\n"
    "# accidentally get called will be shot out of the interpolation region\n"
    "# of the thermal functions so will be add zero.\n"
    "# Renormalization scale Q given as Q^(-2):\n"
    "invQSq = -1.0\n"
    "# The scale is set later (implicitly through invQSq), and may be set\n"
    "# for each evaluation of the potential.\n"
    "\n" << ThermalFunctions::JFunctionsAsPython()
    << "\n"
    "\n" << runningParameters.RunningParametersAsPython()
    << "\n"
    "def TreeLevelPotential( fv ):\n"
    "    return ( " << treeLevelPotential.AsPython() << " )\n"
    "\n"
    "def PolynomialLoopCorrections( fv ):\n"
    "    return ( " << polynomialLoopCorrections.AsPython() << " )\n"
    "\n"
    "def JustLoopCorrection( massesSquaredWithFactors,\n"
    "                        subtractionConstant ):\n"
    "    cumulativeQuantumCorrection = 0.0\n"
    "    for massSquaredWithFactor in massesSquaredWithFactors:\n"
    "        currentQuantumCorrection = 0.0\n"
    "        for MSq in massSquaredWithFactor[ 0 ]:\n"
    "            MM = abs( MSq )\n"
    "            if ( MM > 1.0 ):\n"
    "                currentQuantumCorrection += ( MM**2\n"
    "                                       * ( math.log( MM * invQSq )\n"
    "                                           - subtractionConstant ) )\n"
    "        cumulativeQuantumCorrection += ( massSquaredWithFactor[ 1 ]\n"
    "                                         * currentQuantumCorrection )\n"
    "    return cumulativeQuantumCorrection\n"
    "\n"
    "def LoopAndThermalCorrection( massesSquaredWithFactors,\n"
    "                              subtractionConstant,\n"
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
    "def ScalarMassesSquaredWithFactors( fv ):\n"
    "    massesSquaredWithFactors = []\n";
    for( std::vector< OldRealMassesSquaredMatrix >::const_iterator
         whichMatrix( scalarMassSquaredMatrices.begin() );
         whichMatrix < scalarMassSquaredMatrices.end();
         ++whichMatrix )
    {
      pythonFile << "    massSquaredMatrix = numpy.array( [ ";
      for( std::vector< PolynomialSum >::const_iterator
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
    "def FermionMassesSquaredWithFactors( fv ):\n"
    "    massesSquaredWithFactors = []\n";
    for( std::vector< OldComplexMassSquaredMatrix >::const_iterator
         whichMatrix( fermionMassSquaredMatrices.begin() );
         whichMatrix < fermionMassSquaredMatrices.end();
         ++whichMatrix )
    {
      pythonFile << "    massSquaredMatrix = numpy.array( [ ";
      for( std::vector< std::pair< PolynomialSum,
                                   PolynomialSum > >::const_iterator
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
    for( std::vector< OldSymmetricComplexMassMatrix >::const_iterator
         whichMatrix( fermionMassMatrices.begin() );
         whichMatrix < fermionMassMatrices.end();
         ++whichMatrix )
    {
      pythonFile << "    massMatrix = numpy.array( [ ";
      for( std::vector< std::pair< PolynomialSum,
                                   PolynomialSum > >::const_iterator
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
    "def VectorMassesSquaredWithFactors( fv ):\n"
    "    massesSquaredWithFactors = []\n";
    for( std::vector< OldRealMassesSquaredMatrix >::const_iterator
         whichMatrix( vectorMassSquaredMatrices.begin() );
         whichMatrix < vectorMassSquaredMatrices.end();
         ++whichMatrix )
    {
      pythonFile << "    massSquaredMatrix = numpy.array( [ ";
      for( std::vector< PolynomialSum >::const_iterator
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
    "def JustLoopCorrections( fv ):\n"
    "    totalQuantumCorrections = 0.0\n"
    "    currentCorrection = JustLoopCorrection(\n"
    "                                  ScalarMassesSquaredWithFactors( fv ),\n"
    "                                           1.5 )\n"
    "    totalQuantumCorrections += currentCorrection\n"
    "    currentCorrection = JustLoopCorrection(\n"
    "                                 FermionMassesSquaredWithFactors( fv ),\n"
    "                                           1.5 )\n"
    "    totalQuantumCorrections -= ( 2.0 * currentCorrection )\n"
    "    currentCorrection = JustLoopCorrection(\n"
    "                                  VectorMassesSquaredWithFactors( fv ),\n"
    "                               " << vectorMassCorrectionConstant << " )\n"
    "    totalQuantumCorrections += ( 3.0 * currentCorrection )\n"
    "    return ( totalQuantumCorrections * loopFactor )\n"
    "\n"
    "# The full one-loop potential without thermal corrections:\n"
    "def JustLoopCorrectedPotential( fv ):\n"
    << SetScaleForPythonPotentialCall()
    << "    return ( TreeLevelPotential( fv )\n"
    "             + PolynomialLoopCorrections( fv )\n"
    "             + JustLoopCorrections( fv ) )\n"
    "\n"
    "def LoopAndThermalCorrections( fv ):\n"
    "    totalQuantumCorrections = 0.0\n"
    "    totalThermalCorrections = 0.0\n"
    "    currentCorrections = LoopAndThermalCorrection(\n"
    "                                  ScalarMassesSquaredWithFactors( fv ),\n"
    "                                                  1.5,\n"
    "                                                  BosonicJ )\n"
    "# Real scalar degrees of freedom add to both quantum and thermal\n"
    "# corrections without any factors.\n"
    "    totalQuantumCorrections += currentCorrections[ 0 ]\n"
    "    totalThermalCorrections += currentCorrections[ 1 ]\n"
    "    currentCorrections = LoopAndThermalCorrection(\n"
    "                                 FermionMassesSquaredWithFactors( fv ),\n"
    "                                           1.5,\n"
    "                                           FermionicJ )\n"
    "# Weyl fermion degrees of freedom add to both quantum and thermal\n"
    "# corrections with a factor of 2, though there is an additional minus\n"
    "# sign for the quantum corrections. (The conventions used in Vevacious\n"
    "# are that the thermal correction functions already have any minus\n"
    "# signs for fermions.)\n"
    "    totalQuantumCorrections -= ( 2.0 * currentCorrections[ 0 ] )\n"
    "    totalThermalCorrections += ( 2.0 * currentCorrections[ 1 ] )\n"
    "    currentCorrections= LoopAndThermalCorrection(\n"
    "                                  VectorMassesSquaredWithFactors( fv ),\n"
    "                                " << vectorMassCorrectionConstant << ",\n"
    "                                          BosonicJ )\n"
    "# Vector boson degrees of freedom add to quantum corrections with a\n"
    "# factor of 3 in dimensional regularization schemes, and to thermal\n"
    "# corrections with a factor of 2.\n"
    "    totalQuantumCorrections += ( 3.0 * currentCorrections[ 0 ] )\n"
    "    totalThermalCorrections += ( 2.0 * currentCorrections[ 1 ] )\n"
    "    return ( ( totalQuantumCorrections * loopFactor )\n"
    "             + ( ( totalThermalCorrections * thermalFactor )\n"
    "                 / ( invTSq**2 ) ) )\n"
    "\n"
    "# The full one-loop, thermally-corrected potential:\n"
    "def LoopAndThermallyCorrectedPotential( fv ):\n"
    << SetScaleForPythonPotentialCall()
    << "    return ( TreeLevelPotential( fv )\n"
    "             + PolynomialLoopCorrections( fv )\n"
    "             + LoopAndThermalCorrections( fv ) )\n"
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
  OldPotentialFromPolynomialAndMasses::OldPotentialFromPolynomialAndMasses(
                           RunningParameterManager& runningParameterManager ) :
    PotentialFunction( runningParameterManager ),
    IWritesPythonPotential(),
    runningParameters( runningParameterManager ),
    dsbFieldValuePolynomials(),
    currentMinimumRenormalizationScale( -1.0 ),
    squareOfMinimumRenormalizationScale( -1.0 ),
    currentMaximumRenormalizationScale( -1.0 ),
    squareOfMaximumRenormalizationScale( -1.0 ),
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
    treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions( false ),
    scaleRangeMinimumFactor( -1.0 ),
    fieldsAssumedPositive(),
    fieldsAssumedNegative(),
    assumedPositiveOrNegativeTolerance( -1.0 )
  {
    // This protected constructor is just an initialization list only used by
    // derived classes which are going to fill up the data members in their own
    // constructors.
  }

  // This is just for derived classes.
  OldPotentialFromPolynomialAndMasses::OldPotentialFromPolynomialAndMasses(
                         OldPotentialFromPolynomialAndMasses const& copySource ) :
    PotentialFunction( copySource ),
    IWritesPythonPotential(),
    runningParameters( copySource.runningParameters ),
    dsbFieldValuePolynomials( copySource.dsbFieldValuePolynomials ),
    currentMinimumRenormalizationScale(
                               copySource.currentMinimumRenormalizationScale ),
    squareOfMinimumRenormalizationScale(
                              copySource.squareOfMinimumRenormalizationScale ),
    currentMaximumRenormalizationScale(
                               copySource.currentMaximumRenormalizationScale ),
    squareOfMaximumRenormalizationScale(
                               copySource.squareOfMaximumRenormalizationScale ),
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
    treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions(
          copySource.treeLevelMinimaOnlyAsValidHomotopyContinuationSolutions ),
    scaleRangeMinimumFactor( copySource.scaleRangeMinimumFactor ),
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
  // that the scale has been set correctly.
  double
  OldPotentialFromPolynomialAndMasses::LoopAndThermalCorrections(
                               std::vector< double > const& fieldConfiguration,
          std::vector< DoubleVectorWithDouble > scalarMassesSquaredWithFactors,
         std::vector< DoubleVectorWithDouble > fermionMassesSquaredWithFactors,
          std::vector< DoubleVectorWithDouble > vectorMassesSquaredWithFactors,
                                              double const inverseScaleSquared,
                                          double const temperatureValue ) const
  {
    double inverseTemperatureSquared( 1.0 );
    if( temperatureValue > 0.0 )
    {
      inverseTemperatureSquared
      = ( 1.0 / ( temperatureValue * temperatureValue ) );
    }
    double totalQuantumCorrections( 0.0 );
    double totalThermalCorrections( 0.0 );

    double scalarQuantumCorrections( 0.0 );
    double scalarThermalCorrections( 0.0 );
    AddToCorrections( scalarMassesSquaredWithFactors,
                      inverseScaleSquared,
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
  void OldPotentialFromPolynomialAndMasses::ParseSumOfPolynomialTerms(
                                              std::string const& stringToParse,
                    std::pair< PolynomialSum, PolynomialSum >& polynomialSums )
  {
    polynomialSums.first.PolynomialTerms().clear();
    polynomialSums.second.PolynomialTerms().clear();
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
    if( stringToParse[ 0 ] == '+' )
    {
      wordStart = 1;
    }
    // Now we have skipped any initial '+' or '-', so we start the 1st term.
    PolynomialTerm polynomialTerm;
    if( !positiveTerm )
    {
      polynomialTerm.MultiplyBy( -1.0 );
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
            polynomialSums.second.PolynomialTerms().push_back(
                                                              polynomialTerm );
          }
          else
          {
            polynomialSums.first.PolynomialTerms().push_back( polynomialTerm );
          }
        }
        polynomialTerm.ResetValues();
        if( stringToParse[ wordStart ] == '-' )
        {
          polynomialTerm.MultiplyBy( -1.0 );
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
        polynomialSums.second.PolynomialTerms().push_back( polynomialTerm );
      }
      else
      {
        polynomialSums.first.PolynomialTerms().push_back( polynomialTerm );
      }
    }
  }

  // This reads in a whole number or variable (including possible raising to
  // a power), applies the correct operation to polynomialTerm, and then
  // returns the position of the character just after the interpreted word.
  // If there was a factor of "i", "I", "j", or "J", imaginaryTerm is set to
  // true.
  size_t
  OldPotentialFromPolynomialAndMasses::PutNextNumberOrVariableIntoPolynomial(
                                              std::string const& stringToParse,
                                                              size_t wordStart,
                                               PolynomialTerm& polynomialTerm,
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
        polynomialTerm.MultiplyBy( BOL::StringParser::stringToDouble(
                                               stringToParse.substr( wordStart,
                                                 ( wordEnd - wordStart ) ) ) );
        return wordEnd;
      }
      else if( allowedVariableInitials.find( stringToParse[ wordStart ] )
               != std::string::npos )
      {
        // Otherwise we have to check to see if it is a string that could be a
        // field variable name or a string mapping to a functionoid.
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
        // For comparison, we need the string to be in the proper format:
        std::string
        variableString( FormatVariable( stringToParse.substr( wordStart,
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
            // also demand that no '+' or whitespace chars are allowed.
            throw std::runtime_error(
                              "Model file had invalid exponent after \'^\'." );
          }
          wordStart = ( wordEnd + 1 );
          wordEnd = stringToParse.find_first_not_of( digitChars,
                                                     wordStart );
          powerInt
          = BOL::StringParser::stringToInt( stringToParse.substr( wordStart,
                                                   ( wordEnd - wordStart ) ) );
        }
        // First we check for the imaginary unit:
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
          if( ( powerInt % 2 ) == 1 )
          {
            imaginaryTerm = true;
          }
          if( ( powerInt % 4 ) > 1 )
          {
            polynomialTerm.MultiplyBy( -1.0 );
          }
          return wordEnd;
        }
        // Next we check for a field name:
        for( unsigned int fieldIndex( 0 );
             fieldIndex < fieldNames.size();
             ++fieldIndex )
        {
          if( fieldNames[ fieldIndex ].compare( variableString ) == 0 )
          {
            // If it is a field name, we raise its power in the polynomial term
            // and return position of the char just after the parsed chars.
            polynomialTerm.RaiseFieldPower( fieldIndex,
                                            powerInt );
            return wordEnd;
          }
        }
        // If we get to here, it wasn't a field name, so we try to get a
        // functionoid from runningParameters, which will return NULL if the
        // string doesn't map to any functionoid, and polynomialTerm.MultiplyBy
        // will set its internal validity flag to false on being given NULL.
        polynomialTerm.MultiplyBy( runningParameters.GetFunctionoid(
                                                              variableString ),
                                   powerInt );
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
  void OldPotentialFromPolynomialAndMasses::AddToCorrections(
         std::vector< DoubleVectorWithDouble > const& massesSquaredWithFactors,
                                              double const inverseScaleSquared,
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
        if( massSquared > 1.0 )
        {
          currentQuantumCorrection += ( massSquared * massSquared
                                        * ( log( massSquared
                                                 * inverseScaleSquared )
                                            - subtractFromLogarithm ) );
        }
        if( inverseTemperatureSquared > 0.0 )
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
  std::string OldPotentialFromPolynomialAndMasses::AsDebuggingString() const
  {
    std::stringstream returnStream;
    returnStream << "fieldNames, dsbFieldValuePolynomials:" << std::endl;
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      returnStream << fieldNames[ fieldIndex ] << ": "
      << dsbFieldValuePolynomials[ fieldIndex ].AsDebuggingString()
      << std::endl;
    }
    returnStream
    << std::endl
    << "currentMinimumRenormalizationScale = "
    << currentMinimumRenormalizationScale
    << std::endl
    << "squareOfMinimumRenormalizationScale = "
    << squareOfMinimumRenormalizationScale
    << std::endl
    << "currentMaximumRenormalizationScale = "
    << currentMaximumRenormalizationScale
    << std::endl
    << "squareOfMaximumRenormalizationScale = "
    << squareOfMaximumRenormalizationScale
    << std::endl
    << std::endl
    << "treeLevelPotential = " << treeLevelPotential.AsDebuggingString()
    << std::endl
    << "polynomialLoopCorrections = "
    << polynomialLoopCorrections.AsDebuggingString()
    << std::endl
    << "scalarSquareMasses = " << std::endl;
    for( std::vector< OldRealMassesSquaredMatrix >::const_iterator
         whichScalar( scalarMassSquaredMatrices.begin() );
         whichScalar < scalarMassSquaredMatrices.end();
         ++whichScalar )
    {
      returnStream << whichScalar->AsString();
    }
    returnStream << std::endl;
    returnStream << std::endl
    << "fermionMasses = " << std::endl;
    for( std::vector< OldSymmetricComplexMassMatrix >::const_iterator
         whichFermion( fermionMassMatrices.begin() );
         whichFermion < fermionMassMatrices.end();
         ++whichFermion )
    {
      returnStream << whichFermion->AsString();
    }
    returnStream << std::endl;
    returnStream << std::endl
    << "fermionMassSquareds = " << std::endl;
    for( std::vector< OldComplexMassSquaredMatrix >::const_iterator
         whichFermion( fermionMassSquaredMatrices.begin() );
         whichFermion < fermionMassSquaredMatrices.end();
         ++whichFermion )
    {
      returnStream << whichFermion->AsString();
    }
    returnStream << std::endl;
    returnStream << std::endl
    << "vectorSquareMasses = " << std::endl;
    for( std::vector< OldRealMassesSquaredMatrix >::const_iterator
         whichVector( vectorMassSquaredMatrices.begin() );
         whichVector < vectorMassSquaredMatrices.end();
         ++whichVector )
    {
      returnStream << whichVector->AsString();
    }
    returnStream << std::endl
    << "vectorMassCorrectionConstant = " << vectorMassCorrectionConstant
    << std::endl << "scaleRangeMinimumFactor = " << scaleRangeMinimumFactor
    << std::endl;
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
