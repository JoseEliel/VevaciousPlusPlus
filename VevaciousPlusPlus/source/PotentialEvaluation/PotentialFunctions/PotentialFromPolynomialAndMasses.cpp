/*
 * PotentialFromPolynomialAndMasses.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{
  std::string const
  PotentialFromPolynomialAndMasses::digitChars( "0123456789" );
  std::string const
  PotentialFromPolynomialAndMasses::dotAndDigits( "."
                              + PotentialFromPolynomialAndMasses::digitChars );
  std::string const PotentialFromPolynomialAndMasses::allowedVariableInitials(
                                                   "qwertyuiopasdfghjklzxcvbnm"
                                                "QWERTYUIOPASDFGHJKLZXCVBNM" );
  std::string const PotentialFromPolynomialAndMasses::allowedVariableChars(
                      PotentialFromPolynomialAndMasses::allowedVariableInitials
                                 + PotentialFromPolynomialAndMasses::digitChars
                                                                      + "_~" );
  double const
  PotentialFromPolynomialAndMasses::loopFactor( 1.0 / ( 64.0 * M_PI * M_PI ) );
  double const PotentialFromPolynomialAndMasses::thermalFactor( 1.0
                                                     / ( 2.0 * M_PI * M_PI ) );

  PotentialFromPolynomialAndMasses::PotentialFromPolynomialAndMasses(
                                              std::string const& modelFilename,
                           RunningParameterManager& runningParameterManager ) :
    PotentialFunction( runningParameterManager ),
    IWritesPythonPotential(),
    runningParameters( runningParameterManager ),
    dsbFieldValuePolynomials(),
    currentMinimumRenormalizationScale( NAN ),
    squareOfMinimumRenormalizationScale( NAN ),
    currentMaximumRenormalizationScale( NAN ),
    treeLevelPotential(),
    polynomialLoopCorrections(),
    scalarSquareMasses(),
    fermionSquareMasses(),
    vectorSquareMasses(),
    scalarMassSquaredMatrices(),
    fermionMassMatrices(),
    fermionMassSquaredMatrices(),
    vectorMassSquaredMatrices(),
    vectorMassCorrectionConstant( NAN ),
    needToUpdateHomotopyContinuation( false )
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
                               FormatVariable( elementLines[ lineIndex ] ) ) );
      if( !(readFieldName.empty()) )
      {
        fieldNames.push_back( readFieldName );
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
      unsigned int fieldIndex( FieldIndex(
                                       BOL::StringParser::trimFromFrontAndBack(
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
        int numberOfRows( sqrt( (double)(elementLines.getSize()) ) );
        if( ( numberOfRows * numberOfRows ) != elementLines.getSize() )
        {
          throw std::runtime_error( "Number of elements for"
                     " RealBosonMassSquaredMatrix was not a square integer!" );
        }
        RealMassesSquaredMatrix massSquaredMatrix( numberOfRows,
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
      else if( elementParser.currentElementNameMatches(
                                                    "WeylFermionMassMatrix" ) )
      {
        elementLines.clearEntries();
        BOL::StringParser::parseByChar(
                               elementParser.getTrimmedCurrentElementContent(),
                                        elementLines,
                                        '\n');
        int numberOfRows( sqrt( (double)(elementLines.getSize()) ) );
        if( ( numberOfRows * numberOfRows ) != elementLines.getSize() )
        {
          throw std::runtime_error( "Number of elements for"
                          " WeylFermionMassMatrix was not a square integer!" );
        }
        SymmetricComplexMassMatrix fermionMassMatrix( numberOfRows,
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
        int numberOfRows( sqrt( (double)(elementLines.getSize()) ) );
        if( ( numberOfRows * numberOfRows ) != elementLines.getSize() )
        {
          throw std::runtime_error( "Number of elements for"
            " ComplexWeylFermionMassSquaredMatrix was not a square integer!" );
        }
        ComplexMassSquaredMatrix fermionMassSquaredMatrix( numberOfRows,
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
    for( unsigned int pointerIndex( 0 );
         pointerIndex < scalarMassSquaredMatrices.size();
         ++pointerIndex )
    {
      scalarSquareMasses.push_back(
                                &(scalarMassSquaredMatrices[ pointerIndex ]) );
    }
    for( unsigned int pointerIndex( 0 );
         pointerIndex < fermionMassMatrices.size();
         ++pointerIndex )
    {
      fermionSquareMasses.push_back( &(fermionMassMatrices[ pointerIndex ]) );
    }
    for( unsigned int pointerIndex( 0 );
         pointerIndex < fermionMassSquaredMatrices.size();
         ++pointerIndex )
    {
      fermionSquareMasses.push_back(
                               &(fermionMassSquaredMatrices[ pointerIndex ]) );
    }
    for( unsigned int pointerIndex( 0 );
         pointerIndex < vectorMassSquaredMatrices.size();
         ++pointerIndex )
    {
      vectorSquareMasses.push_back(
                                &(vectorMassSquaredMatrices[ pointerIndex ]) );
    }

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "end of PotentialFromPolynomialAndMasses::"
    << "PotentialFromPolynomialAndMasses( \"" << modelFilename << "\" )"
    << std::endl << "fieldNames, dsbFieldValuePolynomials:" << std::endl;
    for( unsigned int fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      std::cout << fieldNames[ fieldIndex ] << ": "
      << dsbFieldValuePolynomials[ fieldIndex ].AsString() << std::endl;
    }
    std::cout
    << std::endl
    << "renormalizationScaleSquared = " << renormalizationScaleSquared
    << std::endl
    << "minimumRenormalizationScaleSquared = "
    << minimumRenormalizationScaleSquared
    << std::endl
    << "treeLevelPotential = " << treeLevelPotential.AsString()
    << std::endl
    << "polynomialLoopCorrections = " << polynomialLoopCorrections.AsString()
    << std::endl
    << "scalarSquareMasses = " << std::endl;
    for( std::vector< RealMassesSquaredMatrix >::iterator
         whichScalar( scalarSquareMasses.begin() );
         whichScalar < scalarSquareMasses.end();
         ++whichScalar )
    {
      std::cout << whichScalar->AsString();
    }
    std::cout << std::endl;
    std::cout << std::endl
    << "fermionMasses = " << std::endl;
    for( std::vector< SymmetricComplexMassMatrix >::iterator
         whichFermion( fermionMasses.begin() );
         whichFermion < fermionMasses.end();
         ++whichFermion )
    {
      std::cout << whichFermion->AsString();
    }
    std::cout << std::endl;
    std::cout << std::endl
    << "fermionMassSquareds = " << std::endl;
    for( std::vector< ComplexMassSquaredMatrix >::iterator
         whichFermion( fermionMassSquareds.begin() );
         whichFermion < fermionMassSquareds.end();
         ++whichFermion )
    {
      std::cout << whichFermion->AsString();
    }
    std::cout << std::endl;
    std::cout << std::endl
    << "vectorSquareMasses = " << std::endl;
    for( std::vector< RealMassesSquaredMatrix >::iterator
         whichVector( vectorSquareMasses.begin() );
         whichVector < vectorSquareMasses.end();
         ++whichVector )
    {
      std::cout << whichVector->AsString();
    }
    std::cout << std::endl;*/
  }

  PotentialFromPolynomialAndMasses::~PotentialFromPolynomialAndMasses()
  {
    // This does nothing.
  }


  // This writes the potential as
  // def PotentialFunction( fv ): return ...
  // in pythonFilename for fv being an array of floating-point numbers in the
  // same order as they are for the field configurations as internal to this
  // C++ code. It uses the purely virtual function SetScaleInPythonFunction.
  void PotentialFromPolynomialAndMasses::WriteAsPython(
                                       std::string const pythonFilename ) const
  {
    std::ofstream pythonFile( pythonFilename.c_str() );
    pythonFile << std::setprecision( 12 );
    pythonFile << "# Automatically generated file! Modify at your own PERIL!\n"
    "# This file was created by Vevacious version "
    << VevaciousPlusPlus::versionString << "\n"
    "from __future__ import division\n"
    "import math\n"
    "import numpy\n"
    "\n"
    "# Global variables:\n"
    "# Temperature, T:\n"
    "T = 0.0\n"
    "# T^(-2):\n"
    "invTSq = -1.0\n"
    "# Renormalization scale, Q:\n"
    "Q = -1.0\n"
    "# Q^(-2):\n"
    "invQSq = -1.0\n"
    "# Q is set later, and may be set for each evaluation of the potential.\n"
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
    for( std::vector< RealMassesSquaredMatrix >::const_iterator
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
    for( std::vector< ComplexMassSquaredMatrix >::const_iterator
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
    for( std::vector< SymmetricComplexMassMatrix >::const_iterator
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
    for( std::vector< RealMassesSquaredMatrix >::const_iterator
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
    "             + ( totalThermalCorrections * thermalFactor * T**4 ) )\n";


    // placeholder:
    /**/std::cout << std::endl
    << "Placeholder: "
    << "need to write CT functions in .py";
    std::cout << std::endl;/**/

    pythonFile.close();
  }

  // This is just for derived classes.
  PotentialFromPolynomialAndMasses::PotentialFromPolynomialAndMasses(
                           RunningParameterManager& runningParameterManager ) :
    PotentialFunction( runningParameterManager ),
    IWritesPythonPotential(),
    runningParameters( runningParameterManager ),
    dsbFieldValuePolynomials(),
    currentMinimumRenormalizationScale( NAN ),
    squareOfMinimumRenormalizationScale( NAN ),
    currentMaximumRenormalizationScale( NAN ),
    treeLevelPotential(),
    polynomialLoopCorrections(),
    scalarSquareMasses(),
    fermionSquareMasses(),
    vectorSquareMasses(),
    scalarMassSquaredMatrices(),
    fermionMassMatrices(),
    fermionMassSquaredMatrices(),
    vectorMassSquaredMatrices(),
    vectorMassCorrectionConstant( NAN ),
    needToUpdateHomotopyContinuation( false )
  {
    // This protected constructor is just an initialization list only used by
    // derived classes which are going to fill up the data members in their own
    // constructors.
  }

  // This is just for derived classes.
  PotentialFromPolynomialAndMasses::PotentialFromPolynomialAndMasses(
                         PotentialFromPolynomialAndMasses const& copySource ) :
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
    needToUpdateHomotopyContinuation(
                                  copySource.needToUpdateHomotopyContinuation )
  {
    // Now we can fill the MassesSquaredCalculator* vectors, as their pointers
    // should remain valid as the other vectors do not change size any more
    // after the constructor.
    for( unsigned int pointerIndex( 0 );
         pointerIndex < scalarMassSquaredMatrices.size();
         ++pointerIndex )
    {
      scalarSquareMasses.push_back(
                                &(scalarMassSquaredMatrices[ pointerIndex ]) );
    }
    for( unsigned int pointerIndex( 0 );
         pointerIndex < fermionMassMatrices.size();
         ++pointerIndex )
    {
      fermionSquareMasses.push_back( &(fermionMassMatrices[ pointerIndex ]) );
    }
    for( unsigned int pointerIndex( 0 );
         pointerIndex < fermionMassSquaredMatrices.size();
         ++pointerIndex )
    {
      fermionSquareMasses.push_back(
                               &(fermionMassSquaredMatrices[ pointerIndex ]) );
    }
    for( unsigned int pointerIndex( 0 );
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
  PotentialFromPolynomialAndMasses::LoopAndThermalCorrections(
                             std::vector< double > const& fieldConfiguration,
          std::vector< DoubleVectorWithDouble > scalarMassesSquaredWithFactors,
         std::vector< DoubleVectorWithDouble > fermionMassesSquaredWithFactors,
          std::vector< DoubleVectorWithDouble > vectorMassesSquaredWithFactors,
                                              double const inverseScaleSquared,
                                          double const temperatureValue ) const
  {
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "PotentialFromPolynomialAndMasses::LoopAndThermalCorrections( {";
    for( std::vector< double >::const_iterator
         fieldValue( fieldConfiguration.begin() );
         fieldValue < fieldConfiguration.end();
         ++fieldValue )
    {
      std::cout << "  " << *fieldValue;
    }
    std::cout
    << "  }, [scalarMassesSquaredWithFactors, size = "
    << scalarMassesSquaredWithFactors.size()
    << "], [fermionMassesSquaredWithFactors, size = "
    << fermionMassesSquaredWithFactors.size()
    << "], [vectorMassesSquaredWithFactors, size = "
    << vectorMassesSquaredWithFactors.size() << "], inverseScaleSquared = "
    << inverseScaleSquared << ", temperatureValue = " << temperatureValue
    << " ) called.";
    std::cout << std::endl;*/

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

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "scalarQuantumCorrections = " << scalarQuantumCorrections;
    std::cout << std::endl;*/


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

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "-2.0 * fermionQuantumCorrections = "
    << ( -2.0 * fermionQuantumCorrections );
    std::cout << std::endl;*/

    // Weyl fermion degrees of freedom add to both quantum and thermal
    // corrections with a factor of 2, though there is an additional minus sign
    // for the quantum corrections. (The conventions used in Vevacious are that
    // the thermal correction functions already have any minus signs for
    // fermions.)
    totalQuantumCorrections -= fermionQuantumCorrections;
    totalQuantumCorrections -= fermionQuantumCorrections;
    totalThermalCorrections += fermionThermalCorrections;
    totalThermalCorrections += fermionThermalCorrections;

    double vectorQuantumCorrections( 0.0 );
    double vectorThermalCorrections( 0.0 );
    AddToCorrections( vectorMassesSquaredWithFactors,
                      inverseScaleSquared,
                      inverseTemperatureSquared,
                      vectorMassCorrectionConstant,
                      &(ThermalFunctions::BosonicJ),
                      vectorQuantumCorrections,
                      vectorThermalCorrections );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "3.0 * vectorQuantumCorrections = "
    << ( 3.0 * vectorQuantumCorrections );
    std::cout << std::endl;*/

    // Vector boson degrees of freedom add to quantum corrections with a factor
    // of 3 in dimensional regularization schemes, and to thermal corrections
    // with a factor of 2.
    totalQuantumCorrections += vectorQuantumCorrections;
    totalQuantumCorrections += vectorQuantumCorrections;
    totalQuantumCorrections += vectorQuantumCorrections;
    totalThermalCorrections += vectorThermalCorrections;
    totalThermalCorrections += vectorThermalCorrections;

    return ( ( totalQuantumCorrections * loopFactor )
             + ( totalThermalCorrections * thermalFactor
                 * temperatureValue * temperatureValue
                 * temperatureValue * temperatureValue ) );
  }

  // This interprets stringToParse as a sum of polynomial terms and sets
  // polynomialSum accordingly.
  void PotentialFromPolynomialAndMasses::ParseSumOfPolynomialTerms(
                                              std::string const& stringToParse,
                    std::pair< PolynomialSum, PolynomialSum >& polynomialSums )
  {
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "PotentialFromPolynomialAndMasses::ParseSumOfPolynomialTerms( \""
    << stringToParse << "\", ... ) called.";
    std::cout << std::endl;*/

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

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "Just made initial polynomialTerm, =" << std::endl
    << polynomialTerm.AsString();
    std::cout << std::endl;*/


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
        // debugging:
        /*std::cout << std::endl << "debugging:"
        << std::endl
        << "imaginaryTerm = " << imaginaryTerm << std::endl
        << "polynomialTerm:"
        << std::endl << polynomialTerm.AsString();
        std::cout << std::endl;*/

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

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "imaginaryTerm = " << imaginaryTerm << std::endl
    << "polynomialTerm:"
    << std::endl << polynomialTerm.AsString();
    std::cout << std::endl;*/

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
  PotentialFromPolynomialAndMasses::PutNextNumberOrVariableIntoPolynomial(
                                              std::string const& stringToParse,
                                                              size_t wordStart,
                                               PolynomialTerm& polynomialTerm,
                                                          bool& imaginaryTerm )
  {
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "PotentialFromPolynomialAndMasses::"
    << "PutNextNumberOrVariableIntoPolynomial( \"" << stringToParse << "\", "
    << wordStart << ", ..., " << imaginaryTerm << " ) called. stringToParse[ "
    << wordStart << " ] = \'" << stringToParse[ wordStart ] << "\'";
    std::cout << std::endl;*/

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

        // debugging:
        /*std::cout << std::endl << "debugging:"
        << std::endl
        << "multiplied by number: " << stringToParse.substr( wordStart,
                                                      ( wordEnd - wordStart ) )
        << ", returning " << wordEnd;
        std::cout << std::endl;*/

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
          // debugging:
          /*std::cout << std::endl << "debugging:"
          << std::endl
          << "power of imaginary unit = " << powerInt
          << ", returning " << wordEnd;
          std::cout << std::endl;*/

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

            // debugging:
            /*std::cout << std::endl << "debugging:"
            << std::endl
            << "multiplied by field[ " << fieldIndex << " ] = \""
            << fieldNames[ fieldIndex ] << "\" to power of " << powerInt
            << ", returning " << wordEnd;
            std::cout << std::endl;*/

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

        // debugging:
        /*std::cout << std::endl << "debugging:"
        << std::endl
        << "multiplied by functionoid ";
        ParameterFunctionoid* multiplyingFunctionoid(
                          runningParameters.GetFunctionoid( variableString ) );
        if( multiplyingFunctionoid == NULL )
        {
          std::cout << "NULL";
        }
        else
        {
          std::cout << multiplyingFunctionoid->AsString();
        }
        std::cout << " to power of " << powerInt
        << ", returning " << wordEnd;
        std::cout << std::endl;*/

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
  void PotentialFromPolynomialAndMasses::AddToCorrections(
         std::vector< DoubleVectorWithDouble > const& massesSquaredWithFactors,
                                              double const inverseScaleSquared,
                                        double const inverseTemperatureSquared,
                                           double const subtractFromLogarithm,
                                    double (*ThermalFunction)( double const ),
                                           double& cumulativeQuantumCorrection,
                                    double& cumulativeThermalCorrection ) const
  {
    // debugging:
    /*std::cout << std::setprecision( 12 ) << std::endl << "debugging:"
    << std::endl
    << "inverseScaleSquared = " << inverseScaleSquared;
    std::cout << std::endl;*/

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

          // debugging:
          /*std::cout << std::endl << "debugging:"
          << std::endl
          << "massSquared = " << massSquared
          << std::endl
          << "massSquared * inverseScaleSquared = "
          << ( massSquared * inverseScaleSquared )
          << std::endl
          << "log(massSquared * inverseScaleSquared ) = "
          << log( massSquared * inverseScaleSquared )
          << std::endl
          << "log(massSquared * inverseScaleSquared )"
          <<  " - subtractFromLogarithm = "
          << ( log( massSquared * inverseScaleSquared )
               - subtractFromLogarithm )
          << std::endl
          << "total correction = "
          << ( massSquared * massSquared
              * ( log( massSquared
                       * inverseScaleSquared )
                  - subtractFromLogarithm ) );
          std::cout << std::endl;*/
        }
        if( inverseTemperatureSquared > 0.0 )
        {
          currentThermalCorrection += (*ThermalFunction)( massSquared
                                                 * inverseTemperatureSquared );
        }

        // debugging:
        /*std::cout << std::endl << "debugging:"
        << std::endl
        << "*massSquaredIterator = " << *massSquaredIterator
        << ", currentQuantumCorrection = " << currentQuantumCorrection;
        std::cout << std::endl;*/
      }
      cumulativeQuantumCorrection
      += ( massesSquared->second * currentQuantumCorrection );
      cumulativeThermalCorrection
      += ( massesSquared->second * currentThermalCorrection );

      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "cumulativeQuantumCorrection = " << cumulativeQuantumCorrection;
      std::cout << std::endl;*/
    }
  }

} /* namespace VevaciousPlusPlus */
