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
                       "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
  std::string const PotentialFromPolynomialAndMasses::allowedVariableChars(
                      PotentialFromPolynomialAndMasses::allowedVariableInitials
                                 + PotentialFromPolynomialAndMasses::digitChars
                                                                      + "_~" );

  PotentialFromPolynomialAndMasses::PotentialFromPolynomialAndMasses(
                                           std::string const& modelFilename ) :
    HomotopyContinuationReadyPotential(),
    runningParameters(),
    treeLevelPotential(),
    polynomialLoopCorrections(),
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


    BOL::AsciiXmlParser fileParser( false );
    BOL::AsciiXmlParser elementParser( false );
    BOL::VectorlikeArray< std::string > elementLines;
    fileParser.openRootElementOfFile( modelFilename );
    // The model file should always have the XML elements in the correct order!
    // <AllowedNonZeroVariables>
    fileParser.readNextElement();
    elementParser.loadString( fileParser.getCurrentElementContent() );
    //   <FieldVariables>
    elementParser.readNextElement();
    BOL::StringParser::parseByChar(
                               elementParser.getTrimmedCurrentElementContent(),
                                    elementLines,
                                    '\n');
    for( unsigned int lineIndex( 0 );
         lineIndex < elementLines.getSize();
         ++lineIndex )
    {
      fieldNames.push_back( FormatVariable( elementLines[ lineIndex ] ) );
    }
    numberOfFields = fieldNames.size();
    //   </FieldVariables>
    //   <SlhaBlocks>
    elementParser.readNextElement();
    std::string slhaString;
    std::string aliasString;
    BOL::StringParser::parseByChar(
                               elementParser.getTrimmedCurrentElementContent(),
                                    elementLines,
                                    '\n');
    for( unsigned int lineIndex( 0 );
         lineIndex < elementLines.getSize();
         ++lineIndex )
    {
      slhaString = BOL::StringParser::trimFromFrontAndBack(
                     BOL::StringParser::firstWordOf( elementLines[ lineIndex ],
                                                     &aliasString,
                                                     "=" ) );
      runningParameters.AddValidSlhaBlock( slhaString,
                                           aliasString );
    }
    //   </SlhaBlocks>
    //   <DerivedParameters>
    elementParser.readNextElement();
    std::string readableName;
    BOL::StringParser::parseByChar(
                               elementParser.getTrimmedCurrentElementContent(),
                                    elementLines,
                                    '\n');
    for( unsigned int lineIndex( 0 );
         lineIndex < elementLines.getSize();
         ++lineIndex )
    {
      runningParameters.CreateDerivedParameter( elementLines[ lineIndex ] );
    }
    //   </DerivedParameters>
    // </AllowedNonZeroVariables>
    // <TreeLevelPotential>
    fileParser.readNextElement();
    ParseSumOfPolynomialTerms( fileParser.getTrimmedCurrentElementContent(),
                               treeLevelPotential );
    // </TreeLevelPotential>
    // <LoopCorrections>
    fileParser.readNextElement();
    std::string renormalizationScheme(
         fileParser.getCurrentElementAttributes()[ "RenormalizationScheme" ] );
    if( renormalizationScheme.compare( "\"MSBAR\"" ) == 0 )
    {
      vectorMassCorrectionConstant = ( 5.0 / 6.0 );
    }
    elementParser.loadString( fileParser.getCurrentElementContent() );
    //   <ExtraPolynomialPart>
    elementParser.readNextElement();
    ParseSumOfPolynomialTerms( elementParser.getTrimmedCurrentElementContent(),
                               polynomialLoopCorrections );
    //   </ExtraPolynomialPart>
    //   <MassSquaredMatrix> (start of first MassSquaredMatrix)
    while( elementParser.readNextElement() )
    {
      if( !(elementParser.currentElementNameMatches( "MassSquaredMatrix" )) )
      {
        break;
      }
      massSquaredMatrices.push_back( MassSquaredMatrix(
                               elementParser.getCurrentElementAttributes() ) );
      BOL::StringParser::parseByChar(
                               elementParser.getTrimmedCurrentElementContent(),
                                      elementLines,
                                      '\n');
      for( unsigned int lineIndex( 0 );
           lineIndex < elementLines.getSize();
           ++lineIndex )
      {
        ParseSumOfPolynomialTerms( elementLines[ lineIndex ],
                                  massSquaredMatrices.back().AddNewElement() );
      }
    }
    //   </MassSquaredMatrix> (end of last MassSquaredMatrix)
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


  // This puts all index brackets into a consistent form.
  std::string PotentialFromPolynomialAndMasses::FormatVariable(
                                 std::string const& unformattedVariable ) const
  {
    size_t openBracket( unformattedVariable.find( '[' ) );
    if( openBracket == std::string::npos )
    {
      return std::string( unformattedVariable );
    }
    if( unformattedVariable[ unformattedVariable.size() - 1 ] != ']' )
    {
      throw std::runtime_error(
                         "In parsing model file, [...] not closed properly." );
    }
    std::vector< int > indicesVector( BOL::StringParser::stringToIntVector(
                              unformattedVariable.substr(  ( openBracket + 1 ),
                        ( unformattedVariable.size() - openBracket - 2 ) ) ) );
    std::stringstream indicesStream;
    indicesStream << unformattedVariable.substr( 0,
                                                 openBracket );
    indicesStream << '[';
    for( std::vector< int >::iterator
         whichIndex( indicesVector.begin() );
         whichIndex < indicesVector.end();
         ++whichIndex )
    {
      if( whichIndex != indicesVector.begin() )
      {
        indicesStream << ',';
      }
      indicesStream << *whichIndex;
    }
    indicesStream  << ']';
    return std::string( indicesStream.str() );
  }

  // This interprets stringToParse as a sum of polynomial terms and sets
  // polynomialSum accordingly.
  void PotentialFromPolynomialAndMasses::ParseSumOfPolynomialTerms(
                                              std::string const& stringToParse,
                                                 PolynomialSum& polynomialSum )
  {
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "PotentialFromPolynomialAndMasses::ParseSumOfPolynomialTerms( \""
    << stringToParse << "\", ... ) called.";
    std::cout << std::endl;/**/

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
    polynomialTerms.push_back( PolynomialTerm( positiveTerm ) );
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
        polynomialTerms.push_back( PolynomialTerm( ( stringToParse[ wordStart ]
                                                     == '+' ) ) );
      }
      // Now we parse the next word.
      wordStart = PutNextNumberOrVariableIntoPolynomial( stringToParse,
                                                         wordStart,
                                                      polynomialTerms.back() );
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
                == std::string::npos ) )
          {
            wordEnd = stringToParse.find_first_not_of( digitChars,
                                                       ( wordEnd + 1 ) );
          }
          else if( ( wordEnd < ( stringToParse.size() - 3 ) )
                   &&
                   ( ( stringToParse[ wordEnd ] == '+' )
                     ||
                     ( stringToParse[ wordEnd ] == '-' ) )
                   &&
                   ( digitChars.find( stringToParse[ wordEnd + 2 ] )
                   == std::string::npos ) )
          {
            wordEnd = stringToParse.find_first_not_of( digitChars,
                                                       ( wordEnd + 2 ) );
          }
          else
          {
            throw std::runtime_error(
                           "Model file had malformed scientific E notation." );
          }
          polynomialTerm.MultiplyBy( BOL::StringParser::stringToDouble(
                                               stringToParse.substr( wordStart,
                                                 ( wordEnd - wordStart ) ) ) );
          return wordEnd;
        }
      }
      else if( allowedVariableInitials.find( stringToParse[ wordStart ] )
               != std::string::npos )
      {
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
        }
        int powerInt( 1 );
        if( ( wordEnd < stringToParse.size() )
            &&
            ( stringToParse[ wordEnd ] == '^' ) )
        {
          if( ( wordEnd == ( stringToParse.size() - 1 ) )
              ||
              ( digitChars.find( stringToParse[ wordEnd + 1 ] )
                == std::string::npos ) )
          {
            throw std::runtime_error(
                              "Model file had invalid exponent after \'^\'." );
          }
          powerInt = BOL::StringParser::stringToInt(
              stringToParse.substr( ( wordEnd + 1 ),
                                   stringToParse.find_first_not_of( digitChars,
                                                         ( wordEnd + 1 ) ) ) );
        }
        polynomialTerm.MultiplyBy( runningParameters.GetFunctionoid(
                               FormatVariable( stringToParse.substr( wordStart,
                                                 ( wordEnd - wordStart ) ) ) ),
                                   powerInt );
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
