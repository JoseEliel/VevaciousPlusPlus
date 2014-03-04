/*
 * PotentialFromPolynomialAndMasses.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/VevaciousPlusPlus.hpp"

namespace VevaciousPlusPlus
{
  std::string const PotentialFromPolynomialAndMasses::digitChars(
                                               BOL::StringParser::digitChars );
  std::string const
  PotentialFromPolynomialAndMasses::dotAndDigits( "."
                              + PotentialFromPolynomialAndMasses::digitChars );
  std::string const PotentialFromPolynomialAndMasses::allowedVariableInitials(
                                      BOL::StringParser::lowercaseAlphabetChars
                                 + BOL::StringParser::uppercaseAlphabetChars );
  std::string const PotentialFromPolynomialAndMasses::allowedVariableChars(
                      PotentialFromPolynomialAndMasses::allowedVariableInitials
                                 + PotentialFromPolynomialAndMasses::digitChars
                                                                      + "_~" );

  PotentialFromPolynomialAndMasses::PotentialFromPolynomialAndMasses(
                                           std::string const& modelFilename ) :
    HomotopyContinuationReadyPotential(),
    runningParameters(),
    renormalizationScaleSquared( NAN ),
    minimumRenormalizationScaleSquared( NAN ),
    treeLevelPotential(),
    polynomialLoopCorrections(),
    massSquaredMatrices(),
    vectorMassCorrectionConstant( NAN ),
    polynomialGradient(),
    polynomialHessian(),
    scaleSlopeOfGradient()
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
    for( int lineIndex( 0 );
         lineIndex < elementLines.getSize();
         ++lineIndex )
    {
      fieldNames.push_back( FormatVariable( elementLines[ lineIndex ] ) );
    }
    numberOfFields = fieldNames.size();
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
    minimumRenormalizationScaleSquared *= minimumRenormalizationScaleSquared;
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
      runningParameters.AddValidSlhaBlock( slhaString,
                      BOL::StringParser::trimFromFrontAndBack( aliasString ) );
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
      runningParameters.CreateDerivedParameter( elementLines[ lineIndex ] );
    }
    //   </DerivedParameters>
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
    if( renormalizationScheme.compare( "\"MSBAR\"" ) == 0 )
    {
      vectorMassCorrectionConstant = ( 5.0 / 6.0 );
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
      //   <MassSquaredMatrix>
      else if( elementParser.currentElementNameMatches( "MassSquaredMatrix" ) )
      {
        massSquaredMatrices.push_back( MassSquaredMatrix(
                               elementParser.getCurrentElementAttributes() ) );
        elementLines.clearEntries();
        BOL::StringParser::parseByChar(
                               elementParser.getTrimmedCurrentElementContent(),
                                        elementLines,
                                        '\n');
        for( int lineIndex( 0 );
             lineIndex < elementLines.getSize();
             ++lineIndex )
        {
          ParseSumOfPolynomialTerms( elementLines[ lineIndex ],
                                  massSquaredMatrices.back().AddNewElement() );
        }
      }
      //   </MassSquaredMatrix>
    }
    // </LoopCorrections>
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
    << "PotentialFromPolynomialAndMasses::operator()(...)";
    std::cout << std::endl;/**/

    UpdateRenormalizationScale( fieldConfiguration,
                                temperatureValue );
    double returnValue( treeLevelPotential( fieldConfiguration )
                        + polynomialLoopCorrections( fieldConfiguration ) );
    double correctionValue( 0.0 );
    for( std::vector< MassSquaredMatrix >::iterator
         massSquaredMatrix( massSquaredMatrices.begin() );
         massSquaredMatrix < massSquaredMatrices.end();
         ++massSquaredMatrix )
    {
      if( massSquaredMatrix->GetSpinType() == MassSquaredMatrix::scalarBoson )
      {
        correctionValue
        = ScalarBosonCorrection( massSquaredMatrix->MassesSquared(
                                                          fieldConfiguration ),
                                 temperatureValue );
      }
      else if( massSquaredMatrix->GetSpinType()
               == MassSquaredMatrix::weylFermion )
      {
        correctionValue
        = WeylFermionCorrection( massSquaredMatrix->MassesSquared(
                                                          fieldConfiguration ),
                                 temperatureValue );
      }
      else if( massSquaredMatrix->GetSpinType()
               == MassSquaredMatrix::gaugeBoson )
      {
        correctionValue
        = GaugeBosonCorrection( massSquaredMatrix->MassesSquared(
                                                          fieldConfiguration ),
                                temperatureValue );
      }
      else
      {
        correctionValue = 0.0;
      }
      returnValue += ( correctionValue
                       * massSquaredMatrix->MultiplicityFactor() );
    }
    return returnValue;
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
    << "PotentialFromPolynomialAndMasses::UpdateParameters( \"" << slhaFilename
    << "\" )";
    std::cout << std::endl;/**/

    runningParameters.UpdateSlhaParameters( slhaFilename );
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
    renormalizationScaleSquared( NAN ),
    minimumRenormalizationScaleSquared( NAN ),
    treeLevelPotential(),
    polynomialLoopCorrections(),
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


  // This interprets stringToParse as a sum of polynomial terms and sets
  // polynomialSum accordingly.
  void PotentialFromPolynomialAndMasses::ParseSumOfPolynomialTerms(
                                              std::string const& stringToParse,
                                                 PolynomialSum& polynomialSum )
  {
    std::vector< PolynomialTerm >&
    polynomialTerms( polynomialSum.PolynomialTerms() );
    polynomialTerms.clear();
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
    polynomialTerms.push_back( PolynomialTerm() );
    if( !positiveTerm )
    {
      polynomialTerms.back().MultiplyBy( -1.0 );
    }
    wordStart = PutNextNumberOrVariableIntoPolynomial( stringToParse,
                                                       wordStart,
                                                      polynomialTerms.back() );
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
        if( !(polynomialTerms.back().IsValid()) )
        {
          // If the term which we just built was invalid, we remove it from
          // the polynomial sum.
          polynomialTerms.pop_back();
        }
        // Now we make the new term based on whether it is added or subtracted
        // in the sum.
        polynomialTerms.push_back( PolynomialTerm() );
        if( stringToParse[ wordStart ] == '-' )
        {
          polynomialTerms.back().MultiplyBy( -1.0 );
        }
      }
      // Now we parse the next word.
      wordStart = PutNextNumberOrVariableIntoPolynomial( stringToParse,
                                                         (++wordStart),
                                                      polynomialTerms.back() );
    }
    if( !(polynomialTerms.back().IsValid()) )
    {
      // If last term, which we just built in the last iteration of the loop,
      // was invalid, we remove it from the polynomial sum.
      polynomialTerms.pop_back();
    }
  }


  // This reads in a whole number or variable (including possible raising to
  // a power), applies the correct operation to polynomialTerm, and then
  // returns the position of the character just after the interpreted word.
  size_t
  PotentialFromPolynomialAndMasses::PutNextNumberOrVariableIntoPolynomial(
                                              std::string const& stringToParse,
                                                              size_t wordStart,
                                               PolynomialTerm& polynomialTerm )
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
        // First we check for a field name:
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

  // This sets the renormalization scale and broadcasts it to the running
  // parameters.
  void PotentialFromPolynomialAndMasses::UpdateRenormalizationScale(
                                      std::vector< double > fieldConfiguration,
                                           double const evaluationTemperature )
  {
    renormalizationScaleSquared
    = ( minimumRenormalizationScaleSquared
        + ( evaluationTemperature * evaluationTemperature ) );
    for( std::vector< double >::iterator
        whichField( fieldConfiguration.begin() );
        whichField < fieldConfiguration.end();
        ++whichField )
    {
      renormalizationScaleSquared += ( (*whichField) * (*whichField) );
    }
    runningParameters.UpdateRunningParameters(
                                         sqrt( renormalizationScaleSquared ) );

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
