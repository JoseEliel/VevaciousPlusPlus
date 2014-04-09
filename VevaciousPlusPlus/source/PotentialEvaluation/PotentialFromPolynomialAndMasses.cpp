/*
 * PotentialFromPolynomialAndMasses.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

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
    HomotopyContinuationReadyPolynomial( runningParameterManager ),
    runningParameters( runningParameterManager ),
    dsbFieldValuePolynomials(),
    renormalizationScaleSquared( NAN ),
    minimumRenormalizationScaleSquared( NAN ),
    treeLevelPotential(),
    polynomialLoopCorrections(),
    scalarSquareMasses(),
    fermionMasses(),
    fermionMassSquareds(),
    vectorSquareMasses(),
    vectorMassCorrectionConstant( NAN )
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
    fieldOrigin = std::vector< double >( numberOfFields,
                                         0.0 );
    //   </FieldVariables>
    //   <MinimumRenormalizationScale>
    successfullyReadElement = elementParser.readNextElement();
    if( !successfullyReadElement )
    {
      throw
      std::runtime_error( "Could not parse <MinimumRenormalizationScale>." );
    }
    minimumRenormalizationScaleSquared = BOL::StringParser::stringToDouble(
                             elementParser.getTrimmedCurrentElementContent() );
    renormalizationScaleSquared = ( minimumRenormalizationScaleSquared
                                    * minimumRenormalizationScaleSquared );
    minimumRenormalizationScaleSquared = renormalizationScaleSquared;
    //   </MinimumRenormalizationScale>
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
            == MassesSquaredFromPolynomials::gaugeBoson )
        {
          vectorSquareMasses.push_back( massSquaredMatrix );
        }
        else
        {
          scalarSquareMasses.push_back( massSquaredMatrix );
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
        fermionMasses.push_back( SymmetricComplexMassMatrix( numberOfRows,
                               elementParser.getCurrentElementAttributes() ) );
        for( int lineIndex( 0 );
             lineIndex < elementLines.getSize();
             ++lineIndex )
        {
          ParseSumOfPolynomialTerms( BOL::StringParser::trimFromFrontAndBack(
                                                   elementLines[ lineIndex ] ),
                                 fermionMasses.back().ElementAt( lineIndex ) );
        }
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
        fermionMassSquareds.push_back( ComplexMassSquaredMatrix( numberOfRows,
                               elementParser.getCurrentElementAttributes() ) );
        for( int lineIndex( 0 );
             lineIndex < elementLines.getSize();
             ++lineIndex )
        {
          ParseSumOfPolynomialTerms( BOL::StringParser::trimFromFrontAndBack(
                                                   elementLines[ lineIndex ] ),
                           fermionMassSquareds.back().ElementAt( lineIndex ) );
        }
      }
      //   </WeylFermionMassMatrix>
    }
    // </LoopCorrections>

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


  // This is just for derived classes.
  PotentialFromPolynomialAndMasses::PotentialFromPolynomialAndMasses(
                           RunningParameterManager& runningParameterManager ) :
    HomotopyContinuationReadyPolynomial( runningParameterManager ),
    runningParameters( runningParameterManager ),
    dsbFieldValuePolynomials(),
    renormalizationScaleSquared( NAN ),
    minimumRenormalizationScaleSquared( NAN ),
    treeLevelPotential(),
    polynomialLoopCorrections(),
    scalarSquareMasses(),
    fermionMasses(),
    fermionMassSquareds(),
    vectorSquareMasses(),
    vectorMassCorrectionConstant( NAN )
  {
    // This protected constructor is just an initialization list only used by
    // derived classes which are going to fill up the data members in their own
    // constructors.
  }

  // This is just for derived classes.
  PotentialFromPolynomialAndMasses::PotentialFromPolynomialAndMasses(
                               PotentialFromPolynomialAndMasses& copySource ) :
    HomotopyContinuationReadyPolynomial( copySource ),
    runningParameters( copySource.runningParameters ),
    dsbFieldValuePolynomials( copySource.dsbFieldValuePolynomials ),
    renormalizationScaleSquared( copySource.renormalizationScaleSquared ),
    minimumRenormalizationScaleSquared(
                               copySource.minimumRenormalizationScaleSquared ),
    treeLevelPotential( copySource.treeLevelPotential ),
    polynomialLoopCorrections( copySource.polynomialLoopCorrections ),
    scalarSquareMasses( copySource.scalarSquareMasses ),
    fermionMasses( copySource.fermionMasses ),
    fermionMassSquareds( copySource.fermionMassSquareds ),
    vectorSquareMasses( copySource.vectorSquareMasses ),
    vectorMassCorrectionConstant( copySource.vectorMassCorrectionConstant )
  {
    // This is just an initialization-list copy-constructor.
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

  // This evaluates the sum of corrections for the real scalar degrees
  // of freedom with masses-squared given by scalarSquareMasses evaluated for
  // the field configuration given by fieldConfiguration at a temperature
  // given by evaluationTemperature.
  double PotentialFromPolynomialAndMasses::ScalarBosonCorrections(
                               std::vector< double > const& fieldConfiguration,
                                        double const inverseTemperatureSquared,
                                             double const temperatureFourthed )
  {
    double totalQuantumCorrection( 0.0 );
    double currentQuantumCorrection( 0.0 );
    double totalThermalCorrection( 0.0 );
    double currentThermalCorrection( 0.0 );
    double absoluteMassSquared( 0.0 );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "PotentialFromPolynomialAndMasses::ScalarBosonCorrections(...)";
    std::cout << std::endl;*/


    for( std::vector< RealMassesSquaredMatrix >::iterator
         scalarSet( scalarSquareMasses.begin() );
         scalarSet < scalarSquareMasses.end();
         ++scalarSet )
    {
      currentQuantumCorrection = 0.0;
      currentThermalCorrection = 0.0;
      std::vector< double > const&
      massesSquared( scalarSet->MassesSquared( fieldConfiguration ) );

      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "massesSquared.size() =" << massesSquared.size();
      std::cout << std::endl;*/

      for( unsigned int whichMass( 0 );
           whichMass < massesSquared.size();
           ++whichMass )
      {
        absoluteMassSquared = massesSquared[ whichMass ];
        if( absoluteMassSquared < 0.0 )
        {
          absoluteMassSquared = -absoluteMassSquared;
        }
        // debugging:
        /*std::cout << std::endl << "debugging:"
        << std::endl
        << "massesSquared[ " << whichMass << " ] = "
        << massesSquared[ whichMass ];
        std::cout << std::endl;*/

        if( absoluteMassSquared > 1.0 )
        {
          currentQuantumCorrection += ( absoluteMassSquared
                                        * absoluteMassSquared
                                        * ( log( absoluteMassSquared
                                                / renormalizationScaleSquared )
                                            - 1.5 ) );

          // debugging:
          /*std::cout << std::endl << "debugging:"
          << std::endl
          << "currentQuantumCorrection = " << currentQuantumCorrection;
          std::cout << std::endl;*/
        }
        if( inverseTemperatureSquared > 1.0 )
        {
          currentThermalCorrection += bosonThermalFunction( absoluteMassSquared
                                                 * inverseTemperatureSquared );
        }
      }
      totalQuantumCorrection
      += ( scalarSet->MultiplicityFactor() * currentQuantumCorrection );
      totalThermalCorrection
      += ( scalarSet->MultiplicityFactor() * currentThermalCorrection );
    }
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "PotentialFromPolynomialAndMasses::ScalarBosonCorrections( {";
    for( std::vector< double >::const_iterator
         whichField( fieldConfiguration.begin() );
         whichField < fieldConfiguration.end();
         ++whichField )
    {
      std::cout << " " << *whichField;
    }
    std::cout << " } ) finishing." << std::endl
    << "totalQuantumCorrection = " << totalQuantumCorrection << std::endl
    << "( loopFactor * totalQuantumCorrection ) = "
    << ( loopFactor * totalQuantumCorrection ) << std::endl
    << "totalThermalCorrection = " << totalThermalCorrection << std::endl
    << "( thermalFactor * totalThermalCorrection * evaluationTemperature"
    << " * evaluationTemperature * evaluationTemperature"
    << " * evaluationTemperature ) = "
    << ( thermalFactor * totalThermalCorrection
         * temperatureFourthed ) << std::endl
    << "returning = " << ( ( loopFactor * totalQuantumCorrection )
                           + ( thermalFactor * totalThermalCorrection
                               * temperatureFourthed ) ) << std::endl;
    std::cout << std::endl;*/

    return ( ( loopFactor * totalQuantumCorrection )
             + ( thermalFactor * totalThermalCorrection
                 * temperatureFourthed ) );
  }

  // This evaluates the sum of corrections for a set of Weyl fermion degrees
  // of freedom with masses-squared given by fermionMasses evaluated for
  // the field configuration given by fieldConfiguration at a temperature
  // given by evaluationTemperature.
  double PotentialFromPolynomialAndMasses::WeylFermionCorrections(
                               std::vector< double > const& fieldConfiguration,
                                        double const inverseTemperatureSquared,
                                             double const temperatureFourthed )
  {
    double totalQuantumCorrection( 0.0 );
    double currentQuantumCorrection( 0.0 );
    double totalThermalCorrection( 0.0 );
    double currentThermalCorrection( 0.0 );
    double absoluteMassSquared( 0.0 );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "PotentialFromPolynomialAndMasses::WeylFermionCorrections(...)";
    std::cout << std::endl;*/

    for( std::vector< SymmetricComplexMassMatrix >::iterator
         fermionSet( fermionMasses.begin() );
         fermionSet < fermionMasses.end();
         ++fermionSet )
    {
      currentQuantumCorrection = 0.0;
      currentThermalCorrection = 0.0;
      std::vector< double > const&
      massesSquared( fermionSet->MassesSquared( fieldConfiguration ) );

      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "(SymmetricComplexMassMatrix) massesSquared.size() =" << massesSquared.size();
      std::cout << std::endl;*/

      for( unsigned int whichMass( 0 );
           whichMass < massesSquared.size();
           ++whichMass )
      {
        absoluteMassSquared = massesSquared[ whichMass ];
        if( absoluteMassSquared < 0.0 )
        {
          absoluteMassSquared = -absoluteMassSquared;
        }
        // debugging:
        /*std::cout << std::endl << "debugging:"
        << std::endl
        << "massesSquared[ " << whichMass << " ] = "
        << massesSquared[ whichMass ];
        std::cout << std::endl;*/

        if( absoluteMassSquared > 1.0 )
        {
          currentQuantumCorrection += ( absoluteMassSquared
                                        * absoluteMassSquared
                                        * ( log( absoluteMassSquared
                                                / renormalizationScaleSquared )
                                          - 1.5 ) );

          // debugging:
          /*std::cout << std::endl << "debugging:"
          << std::endl
          << "currentQuantumCorrection = " << currentQuantumCorrection;
          std::cout << std::endl;*/
        }
        if( inverseTemperatureSquared > 1.0 )
        {
          // debugging:
          /*std::cout << std::endl << "debugging:"
          << std::endl
          << "*massSquared = " << *massSquared << ", temperatureFourthed = "
          << temperatureFourthed
          << ", ( (*massSquared) * inverseTemperatureSquared ) = "
          << ( (*massSquared) * inverseTemperatureSquared );
          std::cout << std::endl;*/

          currentThermalCorrection
          += fermionThermalFunction( absoluteMassSquared
                                     * inverseTemperatureSquared );
        }
      }
      totalQuantumCorrection += ( fermionSet->MultiplicityFactor()
                                  * currentQuantumCorrection );
      totalThermalCorrection += ( fermionSet->MultiplicityFactor()
                                  * currentThermalCorrection );
    }

    for( std::vector< ComplexMassSquaredMatrix >::iterator
         fermionSet( fermionMassSquareds.begin() );
         fermionSet < fermionMassSquareds.end();
         ++fermionSet )
    {
      currentQuantumCorrection = 0.0;
      currentThermalCorrection = 0.0;
      std::vector< double > const&
      massesSquared( fermionSet->MassesSquared( fieldConfiguration ) );

      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "(ComplexMassSquaredMatrix) massesSquared.size() ="
      << massesSquared.size();
      std::cout << std::endl;*/

      for( unsigned int whichMass( 0 );
           whichMass < massesSquared.size();
           ++whichMass )
      {
        absoluteMassSquared = massesSquared[ whichMass ];
        if( absoluteMassSquared < 0.0 )
        {
          absoluteMassSquared = -absoluteMassSquared;
        }
        // debugging:
        /*std::cout << std::endl << "debugging:"
        << std::endl
        << "massesSquared[ " << whichMass << " ] = "
        << massesSquared[ whichMass ];
        std::cout << std::endl;*/

        if( absoluteMassSquared > 1.0 )
        {
          currentQuantumCorrection += ( absoluteMassSquared
                                        * absoluteMassSquared
                                        * ( log( absoluteMassSquared
                                                / renormalizationScaleSquared )
                                            - 1.5 ) );

          // debugging:
          /*std::cout << std::endl << "debugging:"
          << std::endl
          << "currentQuantumCorrection = " << currentQuantumCorrection;
          std::cout << std::endl;*/
        }
        if( inverseTemperatureSquared > 1.0 )
        {
          // debugging:
          /*std::cout << std::endl << "debugging:"
          << std::endl
          << "*massSquared = " << *massSquared << ", temperatureFourthed = "
          << temperatureFourthed
          << ", ( (*massSquared) * inverseTemperatureSquared ) = "
          << ( (*massSquared) * inverseTemperatureSquared );
          std::cout << std::endl;*/

          currentThermalCorrection
          += fermionThermalFunction( absoluteMassSquared
                                     * inverseTemperatureSquared );
        }
      }
      totalQuantumCorrection += ( fermionSet->MultiplicityFactor()
                                  * currentQuantumCorrection );
      totalThermalCorrection += ( fermionSet->MultiplicityFactor()
                                  * currentThermalCorrection );
    }
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "PotentialFromPolynomialAndMasses::WeylFermionCorrections( {";
    for( std::vector< double >::const_iterator
         whichField( fieldConfiguration.begin() );
         whichField < fieldConfiguration.end();
         ++whichField )
    {
      std::cout << " " << *whichField;
    }
    std::cout << " } ) finishing." << std::endl
    << "totalQuantumCorrection = " << totalQuantumCorrection << std::endl
    << "( loopFactor * totalQuantumCorrection ) = "
    << ( loopFactor * totalQuantumCorrection ) << std::endl
    << "totalThermalCorrection = " << totalThermalCorrection << std::endl
    << "( thermalFactor * totalThermalCorrection * temperatureFourthed ) = "
    << ( thermalFactor * totalThermalCorrection
         * temperatureFourthed ) << std::endl
    << "returning = " << ( 2.0 * ( ( thermalFactor * totalThermalCorrection
        * temperatureFourthed )
      - ( loopFactor * totalQuantumCorrection ) ) ) << std::endl;
    std::cout << std::endl;*/

    // The thermal corrections already account for the sign of fermionic
    // contributions, while the zero-temperature corrections need to be
    // subtracted.
    return ( 2.0 * ( ( thermalFactor * totalThermalCorrection
                       * temperatureFourthed )
                     - ( loopFactor * totalQuantumCorrection ) ) );
  }

  // This evaluates the sum of corrections for the real scalar degrees
  // of freedom with masses-squared given by vectorSquareMasses evaluated for
  // the field configuration given by fieldConfiguration at a temperature
  // given by evaluationTemperature.
  double PotentialFromPolynomialAndMasses::GaugeBosonCorrections(
                               std::vector< double > const& fieldConfiguration,
                                        double const inverseTemperatureSquared,
                                             double const temperatureFourthed )
  {
    double totalQuantumCorrection( 0.0 );
    double currentQuantumCorrection( 0.0 );
    double totalThermalCorrection( 0.0 );
    double currentThermalCorrection( 0.0 );
    double absoluteMassSquared( 0.0 );

    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "PotentialFromPolynomialAndMasses::GaugeBosonCorrections(...)";
    std::cout << std::endl;*/

    for( std::vector< RealMassesSquaredMatrix >::iterator
         vectorSet( vectorSquareMasses.begin() );
         vectorSet < vectorSquareMasses.end();
         ++vectorSet )
    {
      currentQuantumCorrection = 0.0;
      currentThermalCorrection = 0.0;
      std::vector< double > const&
      massesSquared( vectorSet->MassesSquared( fieldConfiguration ) );

      // debugging:
      /*std::cout << std::endl << "debugging:"
      << std::endl
      << "massesSquared.size() ="
      << massesSquared.size();
      std::cout << std::endl;*/

      for( unsigned int whichMass( 0 );
           whichMass < massesSquared.size();
           ++whichMass )
      {
        absoluteMassSquared = massesSquared[ whichMass ];
        if( absoluteMassSquared < 0.0 )
        {
          absoluteMassSquared = -absoluteMassSquared;
        }
        // debugging:
        /*std::cout << std::endl << "debugging:"
        << std::endl
        << "massesSquared[ " << whichMass << " ] = "
        << massesSquared[ whichMass ];
        std::cout << std::endl;*/

        if( absoluteMassSquared > 1.0 )
        {
          currentQuantumCorrection += ( absoluteMassSquared
                                        * absoluteMassSquared
                                        * ( log( absoluteMassSquared
                                               / renormalizationScaleSquared )
                                          - vectorMassCorrectionConstant ) );

          // debugging:
          /*std::cout << std::endl << "debugging:"
          << std::endl
          << "currentQuantumCorrection = " << currentQuantumCorrection;
          std::cout << std::endl;*/
        }
        if( inverseTemperatureSquared > 0.0 )
        {
          currentThermalCorrection += bosonThermalFunction( absoluteMassSquared
                                                 * inverseTemperatureSquared );
        }
      }
      totalQuantumCorrection
      += ( vectorSet->MultiplicityFactor() * currentQuantumCorrection );
      totalThermalCorrection
      += ( vectorSet->MultiplicityFactor() * currentThermalCorrection );
    }
    // debugging:
    /*std::cout << std::endl << "debugging:"
    << std::endl
    << "PotentialFromPolynomialAndMasses::GaugeBosonCorrections( {";
    for( std::vector< double >::const_iterator
         whichField( fieldConfiguration.begin() );
         whichField < fieldConfiguration.end();
         ++whichField )
    {
      std::cout << " " << *whichField;
    }
    std::cout << " } ) finishing." << std::endl
    << "totalQuantumCorrection = " << totalQuantumCorrection << std::endl
    << "( loopFactor * totalQuantumCorrection ) = "
    << ( loopFactor * totalQuantumCorrection ) << std::endl
    << "totalThermalCorrection = " << totalThermalCorrection << std::endl
    << "( thermalFactor * totalThermalCorrection * temperatureFourthed ) = "
    << ( thermalFactor * totalThermalCorrection
         * temperatureFourthed ) << std::endl
    << "returning = " << ( ( 3.0 * loopFactor * totalQuantumCorrection )
        + ( 2.0 * thermalFactor * totalThermalCorrection
            * temperatureFourthed ) ) << std::endl;
    std::cout << std::endl;*/

    // Note that vectors have different degrees of freedom for the zero-
    // and non-zero-temperature corrections!
    return ( ( 3.0 * loopFactor * totalQuantumCorrection )
             + ( 2.0 * thermalFactor * totalThermalCorrection
                 * temperatureFourthed ) );
  }

} /* namespace VevaciousPlusPlus */
