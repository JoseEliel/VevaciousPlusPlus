/*
 * PotentialFromPolynomialWithMasses.cpp
 *
 *  Created on: Oct 09, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "PotentialEvaluation/PotentialFunctions/PotentialFromPolynomialWithMasses.hpp"

#include "VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{
  std::string const
  PotentialFromPolynomialWithMasses::digitChars( "0123456789" );
  std::string const
  PotentialFromPolynomialWithMasses::dotAndDigits( "."
                             + PotentialFromPolynomialWithMasses::digitChars );
  std::string const PotentialFromPolynomialWithMasses::allowedVariableInitials(
                                                   "qwertyuiopasdfghjklzxcvbnm"
                                                "QWERTYUIOPASDFGHJKLZXCVBNM" );
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
    PotentialFunction(),
    IWritesPythonPotential(),
    lagrangianParameterManager( lagrangianParameterManager ),
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
    // <FieldVariables>
    successfullyReadElement = fileParser.readNextElement();
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
          fieldNames.push_back( lagrangianParameterManager.FormatVariable(
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
          fieldNames.push_back( lagrangianParameterManager.FormatVariable(
                                       BOL::StringParser::trimFromFrontAndBack(
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
    // <DsbMinimum>
    successfullyReadElement = fileParser.readNextElement();
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

      readFieldName.assign( BOL::StringParser::trimFromFrontAndBack(
                                           elementLines[ lineIndex ].substr( 0,
                                                          equalsPosition ) ) );
      size_t fieldIndex( FieldIndex( lagrangianParameterManager.FormatVariable(
                                                           readFieldName ) ) );
      if( !( fieldIndex < fieldNames.size() ) )
      {
        std::string errorMessage( "Unknown field (\"" );
        errorMessage.append( elementLines[ lineIndex ].substr( 0,
                                                            equalsPosition ) );
        errorMessage.append( "\") given value in DsbMinimum!" );
        throw std::runtime_error( errorMessage );
      }
      dsbFieldInputStrings[ fieldIndex ]
      = lagrangianParameterManager.FormatVariable(
                                       BOL::StringParser::trimFromFrontAndBack(
                   elementLines[ lineIndex ].substr( equalsPosition + 1 ) ) ) ;
    }
    // </DsbMinimum>
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
        int numberOfRows(
                     sqrt( static_cast< double >( elementLines.getSize() ) ) );
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
        int numberOfRows(
                     sqrt( static_cast< double >( elementLines.getSize() ) ) );
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


  // This is just for derived classes.
  PotentialFromPolynomialWithMasses::PotentialFromPolynomialWithMasses(
                     LagrangianParameterManager& lagrangianParameterManager ) :
    PotentialFunction(),
    IWritesPythonPotential(),
    lagrangianParameterManager( lagrangianParameterManager ),
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
    IWritesPythonPotential(),
    lagrangianParameterManager( copySource.lagrangianParameterManager ),
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
    ParametersAndFieldsProduct polynomialTerm;
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
                                    ParametersAndFieldsProduct& polynomialTerm,
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
        polynomialTerm.MultiplyByConstant( BOL::StringParser::stringToDouble(
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
        // For comparison, we need the string to be in the proper format:
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
            polynomialTerm.MultiplyByConstant( -1.0 );
          }
          return wordEnd;
        }
        // Next we check for a field name:
        for( size_t fieldIndex( 0 );
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
        // If we get to here, it wasn't a field name, so we check whether it's
        // a Lagrangian parameter known to lagrangianParameterManager.
        std::pair< bool, size_t >
        parameterValidityAndIndex(
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
