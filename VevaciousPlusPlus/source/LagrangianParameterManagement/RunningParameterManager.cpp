/*
 * RunningParameterManager.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../../include/LagrangianParameterManagement/RunningParameterManager.hpp"

namespace VevaciousPlusPlus
{

  RunningParameterManager::RunningParameterManager() :
    SlhaManager(),
    parameterFunctionoidPointers(),
    parameterFunctionoidMap(),
    parameterFinder(),
    constantFunctionoidMap(),
    constantFinder(),
    slhaParser(),
    slhaBlockPointers(),
    slhaBlockMap(),
    slhaBlockFinder(),
    slhaAliasMap(),
    slhaAliasFinder(),
    lowestBlockScale( -1.0 ),
    highestBlockScale( -1.0 )
  {
    // This constructor is just an initialization list.
  }

  RunningParameterManager::~RunningParameterManager()
  {
    for( int deletionIndex( parameterFunctionoidPointers.size() - 1 );
         deletionIndex >= 0;
         --deletionIndex )
    {
      delete parameterFunctionoidPointers.at( deletionIndex );
    }
    for( int deletionIndex( slhaBlockPointers.size() - 1 );
         deletionIndex >= 0;
         --deletionIndex )
    {
      delete slhaBlockPointers.at( deletionIndex );
    }
  }


  // This adds to the list of SLHA blocks that will be read in when updating
  // running parameters. An alias for the block can be provided.
  void
  RunningParameterManager::AddValidSlhaBlock( std::string slhaBlock,
                                              std::string const& aliasString )
  {
    BOL::StringParser::transformToUppercase( slhaBlock );
    bool notAlreadyDeclared( slhaBlockMap.insert(
                     std::pair< std::string, RestrictedSlhaBlock* >( slhaBlock,
                                                             NULL ) ).second );
    if( !notAlreadyDeclared )
    {
      throw std::runtime_error(
              "Multiple declarations of the same SLHA block in model file!" );
    }
    if( !(aliasString.empty()) )
    {
      notAlreadyDeclared = slhaAliasMap.insert(
                            std::pair< std::string, std::string >( aliasString,
                                                        slhaBlock ) ).second;
      if( !notAlreadyDeclared )
      {
        throw std::runtime_error(
                        "Same alias for multiple SLHA blocks in model file!" );
      }
    }
  }

  // This creates a new ParameterFunctionoid and a string mapping to it.
  void RunningParameterManager::CreateDerivedParameter(
                                           std::string const& parameterString )
  {
    size_t wordEnd( parameterString.find_first_of( '=' ) );
    if( ( wordEnd == std::string::npos )
        ||
        ( parameterString[ parameterString.size() - 1 ] != ']' ) )
    {
      throw
      std::runtime_error( "Model file had malformed derived parameter!" );
    }
    std::string aliasString( BOL::StringParser::trimFromFrontAndBack(
                                                     parameterString.substr( 0,
                                                                 wordEnd ) ) );
    size_t wordStart( wordEnd );
    wordEnd = parameterString.find_first_of( '[',
                                             wordStart );
    // The extra (+/-) 1 and 2 in the substr arguments just remove the '[' and
    // ']' from being included.
    std::string typeString( BOL::StringParser::trimFromFrontAndBack(
                                     parameterString.substr( ( wordStart + 1 ),
                                             ( wordEnd - wordStart - 1 ) ) ) );
    std::string bracketedString( parameterString.substr( ( wordEnd + 1 ),
                                 ( parameterString.size() - wordEnd - 2 ) ) );

    ParameterFunctionoid* createdFunctionoid( NULL );
    bool notAlreadyDeclared( true );
    if( typeString.compare( "SLHAELEMENT" ) == 0 )
    {
      wordEnd = bracketedString.find_first_of( ',' );
      createdFunctionoid
      = GetSlhaFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                                                     bracketedString.substr( 0,
                                                                   wordEnd ) ),
              FormatIndexBracketContent( bracketedString.substr( wordEnd ) ),
                            parameterString );
      if( createdFunctionoid == NULL )
      {
        throw std::runtime_error( "Derived SLHAELEMENT parameter referenced"
                                  " undeclared SLHA block in model file!" );
      }
    }
    else if( typeString.compare( "PLUS" ) == 0 )
    {
      std::pair< ParameterFunctionoid*, ParameterFunctionoid* >
      pointerPair( FindFunctionoidPair( bracketedString ) );
      if( pointerPair.second == NULL )
      {
        throw std::runtime_error(
                           "Derived PLUS parameter incorrect in model file!" );
      }
      createdFunctionoid = new BinaryOperationFunctionoid(
                                   &(BinaryOperationFunctionoid::PlusFunction),
                                                           pointerPair.first,
                                                           pointerPair.second,
                                                           parameterString,
                                                           FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "MINUS" ) == 0 )
    {
      std::pair< ParameterFunctionoid*, ParameterFunctionoid* >
      pointerPair( FindFunctionoidPair( bracketedString ) );
      if( pointerPair.second == NULL )
      {
        throw std::runtime_error(
                          "Derived MINUS parameter incorrect in model file!" );
      }
      createdFunctionoid = new BinaryOperationFunctionoid(
                                  &(BinaryOperationFunctionoid::MinusFunction),
                                                           pointerPair.first,
                                                           pointerPair.second,
                                                           parameterString,
                                                           FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "TIMES" ) == 0 )
    {
      std::pair< ParameterFunctionoid*, ParameterFunctionoid* >
      pointerPair( FindFunctionoidPair( bracketedString ) );
      if( pointerPair.second == NULL )
      {
        throw std::runtime_error(
                          "Derived TIMES parameter incorrect in model file!" );
      }
      createdFunctionoid = new BinaryOperationFunctionoid(
                               &(BinaryOperationFunctionoid::MultiplyFunction),
                                                           pointerPair.first,
                                                           pointerPair.second,
                                                           parameterString,
                                                           FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "DIVIDE" ) == 0 )
    {
      std::pair< ParameterFunctionoid*, ParameterFunctionoid* >
      pointerPair( FindFunctionoidPair( bracketedString ) );
      if( pointerPair.second == NULL )
      {
        throw std::runtime_error(
                         "Derived DIVIDE parameter incorrect in model file!" );
      }
      createdFunctionoid = new BinaryOperationFunctionoid(
                                 &(BinaryOperationFunctionoid::DivideFunction),
                                                           pointerPair.first,
                                                           pointerPair.second,
                                                           parameterString,
                                                           FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "POW" ) == 0 )
    {
      std::pair< ParameterFunctionoid*, ParameterFunctionoid* >
      pointerPair( FindFunctionoidPair( bracketedString ) );
      if( pointerPair.second == NULL )
      {
        throw std::runtime_error(
                            "Derived POW parameter incorrect in model file!" );
      }
      createdFunctionoid = new BinaryOperationFunctionoid( &pow,
                                                           pointerPair.first,
                                                           pointerPair.second,
                                                           parameterString,
                                                           FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "IFNONZERO" ) == 0 )
    {
      std::pair< ParameterFunctionoid*, ParameterFunctionoid* >
      pointerPair( FindFunctionoidPair( bracketedString ) );
      if( pointerPair.second == NULL )
      {
        throw std::runtime_error(
                      "Derived IFNONZERO parameter incorrect in model file!" );
      }
      createdFunctionoid = new BinaryOperationFunctionoid(
                              &(BinaryOperationFunctionoid::IfNonZeroFunction),
                                                           pointerPair.first,
                                                           pointerPair.second,
                                                           parameterString,
                                                           FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "SQRT" ) == 0 )
    {
      ParameterFunctionoid*
      foundPointer( GetFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                                                         bracketedString ) ) );
      if( foundPointer == NULL )
      {
        throw std::runtime_error(
                           "Derived SQRT parameter incorrect in model file!" );
      }
      createdFunctionoid = new UnaryOperationFunctionoid( &sqrt,
                                                          foundPointer,
                                                          parameterString,
                                                          FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "EXP" ) == 0 )
    {
      ParameterFunctionoid*
      foundPointer( GetFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                                                         bracketedString ) ) );
      if( foundPointer == NULL )
      {
        throw std::runtime_error(
                            "Derived EXP parameter incorrect in model file!" );
      }
      createdFunctionoid = new UnaryOperationFunctionoid( &exp,
                                                          foundPointer,
                                                          parameterString,
                                                          FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "LOG" ) == 0 )
    {
      ParameterFunctionoid*
      foundPointer( GetFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                                                         bracketedString ) ) );
      if( foundPointer == NULL )
      {
        throw std::runtime_error(
                            "Derived LOG parameter incorrect in model file!" );
      }
      createdFunctionoid = new UnaryOperationFunctionoid( &log,
                                                          foundPointer,
                                                          parameterString,
                                                          FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "SIN" ) == 0 )
    {
      ParameterFunctionoid*
      foundPointer( GetFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                                                         bracketedString ) ) );
      if( foundPointer == NULL )
      {
        throw std::runtime_error(
                            "Derived SIN parameter incorrect in model file!" );
      }
      createdFunctionoid = new UnaryOperationFunctionoid( &sin,
                                                          foundPointer,
                                                          parameterString,
                                                          FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "COS" ) == 0 )
    {
      ParameterFunctionoid*
      foundPointer( GetFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                                                         bracketedString ) ) );
      if( foundPointer == NULL )
      {
        throw std::runtime_error(
                            "Derived COS parameter incorrect in model file!" );
      }
      createdFunctionoid = new UnaryOperationFunctionoid( &cos,
                                                          foundPointer,
                                                          parameterString,
                                                          FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "TAN" ) == 0 )
    {
      ParameterFunctionoid*
      foundPointer( GetFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                                                         bracketedString ) ) );
      if( foundPointer == NULL )
      {
        throw std::runtime_error(
                            "Derived TAN parameter incorrect in model file!" );
      }
      createdFunctionoid = new UnaryOperationFunctionoid( &tan,
                                                          foundPointer,
                                                          parameterString,
                                                          FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "ASIN" ) == 0 )
    {
      ParameterFunctionoid*
      foundPointer( GetFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                                                         bracketedString ) ) );
      if( foundPointer == NULL )
      {
        throw std::runtime_error(
                           "Derived ASIN parameter incorrect in model file!" );
      }
      createdFunctionoid = new UnaryOperationFunctionoid( &asin,
                                                          foundPointer,
                                                          parameterString,
                                                          FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "ACOS" ) == 0 )
    {
      ParameterFunctionoid*
      foundPointer( GetFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                                                         bracketedString ) ) );
      if( foundPointer == NULL )
      {
        throw std::runtime_error(
                           "Derived ACOS parameter incorrect in model file!" );
      }
      createdFunctionoid = new UnaryOperationFunctionoid( &acos,
                                                          foundPointer,
                                                          parameterString,
                                                          FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "ATAN" ) == 0 )
    {
      ParameterFunctionoid*
      foundPointer( GetFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                                                         bracketedString ) ) );
      if( foundPointer == NULL )
      {
        throw std::runtime_error(
                           "Derived ATAN parameter incorrect in model file!" );
      }
      createdFunctionoid = new UnaryOperationFunctionoid( &atan,
                                                          foundPointer,
                                                          parameterString,
                                                          FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "SINH" ) == 0 )
    {
      ParameterFunctionoid*
      foundPointer( GetFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                                                         bracketedString ) ) );
      if( foundPointer == NULL )
      {
        throw std::runtime_error(
                           "Derived SINH parameter incorrect in model file!" );
      }
      createdFunctionoid = new UnaryOperationFunctionoid( &sinh,
                                                          foundPointer,
                                                          parameterString,
                                                          FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "COSH" ) == 0 )
    {
      ParameterFunctionoid*
      foundPointer( GetFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                                                         bracketedString ) ) );
      if( foundPointer == NULL )
      {
        throw std::runtime_error(
                           "Derived COSH parameter incorrect in model file!" );
      }
      createdFunctionoid = new UnaryOperationFunctionoid( &cosh,
                                                          foundPointer,
                                                          parameterString,
                                                          FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    else if( typeString.compare( "TANH" ) == 0 )
    {
      ParameterFunctionoid*
      foundPointer( GetFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                                                         bracketedString ) ) );
      if( foundPointer == NULL )
      {
        throw std::runtime_error(
                           "Derived TANH parameter incorrect in model file!" );
      }
      createdFunctionoid = new UnaryOperationFunctionoid( &tanh,
                                                          foundPointer,
                                                          parameterString,
                                                          FunctionoidName() );
      parameterFunctionoidPointers.push_back( createdFunctionoid );
    }
    notAlreadyDeclared = parameterFunctionoidMap.insert(
                  std::pair< std::string, ParameterFunctionoid* >( aliasString,
                                                 createdFunctionoid ) ).second;
    if( !notAlreadyDeclared )
    {
      throw std::runtime_error(
                 "Same alias for multiple derived parameters in model file!" );
    }
  }

  // This returns a pointer to the ParameterFunctionoid mapped to by
  // parameterString, returning NULL if there is none.
  ParameterFunctionoid*
  RunningParameterManager::GetFunctionoid( std::string const& parameterString )
  {
    parameterFinder = parameterFunctionoidMap.find( parameterString );
    if( parameterFinder != parameterFunctionoidMap.end() )
    {
      return parameterFinder->second;
    }
    double constantValue( 0.0 );
    if( BOL::StringParser::stringIsDouble( parameterString,
                                           constantValue ) )
    {
      constantFinder = constantFunctionoidMap.find( constantValue );
      if( constantFinder != constantFunctionoidMap.end() )
      {
        return constantFinder->second;
      }
      parameterFunctionoidPointers.push_back( new ConstantFunctionoid(
                                                                 constantValue,
                                                               parameterString,
                                                         FunctionoidName() ) );
      constantFunctionoidMap.insert(
                     std::pair< double, ParameterFunctionoid* >( constantValue,
                                       parameterFunctionoidPointers.back() ) );
      return parameterFunctionoidPointers.back();
    }
    size_t bracketPosition( parameterString.find_first_of( '[' ) );
    if( bracketPosition < ( parameterString.size() - 2 ) )
    {
      return GetSlhaFunctionoid( parameterString.substr( 0,
                                                         bracketPosition ),
                               parameterString.substr( ( bracketPosition + 1 ),
                            ( parameterString.size() - bracketPosition - 2 ) ),
                                 parameterString );
    }
    return NULL;
  }

  // This updates all SlhaFunctionoids with data from a new SLHA file.
  void RunningParameterManager::ReadFile( std::string const& slhaFilename )
  {
    slhaParser.readFile( slhaFilename );
    lowestBlockScale = -1.0;
    highestBlockScale = -1.0;
    for( std::vector< RestrictedSlhaBlock* >::iterator
         whichBlock( slhaBlockPointers.begin() );
         whichBlock < slhaBlockPointers.end();
         ++whichBlock )
    {
      if( (*(*whichBlock))[ 0 ].getScale() > 0.0 )
      {
        if( ( lowestBlockScale <= 0.0 )
            ||
            ( lowestBlockScale > (*(*whichBlock))[ 0 ].getScale() ) )
        {
          lowestBlockScale = (*(*whichBlock))[ 0 ].getScale();
        }
        if( ( highestBlockScale <= 0.0 )
            ||
            ( highestBlockScale < (*(*whichBlock))[ 0 ].getScale() ) )
        {
          highestBlockScale = (*(*whichBlock))[ 0 ].getScale();
        }
      }
    }
    for( std::vector< SlhaFunctionoid* >::iterator
         whichParameter( slhaFunctionoidPointers.begin() );
         whichParameter < slhaFunctionoidPointers.end();
         ++whichParameter )
    {
      (*whichParameter)->UpdateForNewSlhaParameters();
    }
    double const logarithmOfScale( log( lowestBlockScale ) );
    for( std::vector< ParameterFunctionoid* >::iterator
         whichParameter( parameterFunctionoidPointers.begin() );
         whichParameter < parameterFunctionoidPointers.end();
         ++whichParameter )
    {
      (*whichParameter)->UpdateForNewLogarithmOfScale( logarithmOfScale );
      (*whichParameter)->updateObservers();
    }
  }

  // This returns a pointer to the ParameterFunctionoid for the block named
  // blockName with indices given by indexString, making a new one if there
  // is none already. If the block name was invalid though, nothing is made
  // and NULL is returned.
  ParameterFunctionoid*
  RunningParameterManager::GetSlhaFunctionoid( std::string const& blockName,
                                               std::string const& indexString,
                                           std::string const& parameterString )
  {
    ParameterFunctionoid* returnPointer( NULL );

    // First we check to see if blockName was an alias:
    std::string const* properBlockName( &blockName );
    slhaAliasFinder = slhaAliasMap.find( blockName );
    if( slhaAliasFinder != slhaAliasMap.end() )
    {
      properBlockName = &(slhaAliasFinder->second);
    }

    slhaBlockFinder = slhaBlockMap.find( *properBlockName );
    if( slhaBlockFinder != slhaBlockMap.end() )
    {
      // If *properBlockName was a valid declared SLHA block, we form the
      // string which the functionoid should be mapped to.
      std::string lookupString( *properBlockName );
      lookupString.append( "[" );
      lookupString.append( indexString );
      lookupString.append( "]" );
      // We now check to see if the proper mapping string already exists
      // (which may happen if the original blockName was an alias).
      parameterFinder = parameterFunctionoidMap.find( lookupString );
      if( parameterFinder != parameterFunctionoidMap.end() )
      {
        // If the mapping did exist for the alias, we insert the alias into the
        // map before returning the pointer.
        lookupString = blockName;
        lookupString.append( "[" );
        lookupString.append( indexString );
        lookupString.append( "]" );
        parameterFunctionoidMap.insert(
                 std::pair< std::string, ParameterFunctionoid* >( lookupString,
                                                 parameterFinder->second ) );
        return parameterFinder->second;
      }

      // Otherwise we need to make a new functionoid. We make the functionoid
      // before its SLHA block because it sets up the number of indices which
      // is needed for the block constructor.
      SlhaFunctionoid* slhaFunctionoid = new SlhaFunctionoid( indexString,
                                                             parameterString,
                                                           FunctionoidName() );
      parameterFunctionoidPointers.push_back( slhaFunctionoid );
      parameterFunctionoidMap.insert(
                 std::pair< std::string, ParameterFunctionoid* >( lookupString,
                                                           slhaFunctionoid ) );
      slhaFunctionoidPointers.push_back( slhaFunctionoid );
      if( slhaBlockFinder->second == NULL )
      {
        slhaBlockFinder->second = new RestrictedSlhaBlock( *properBlockName,
                                            slhaFunctionoid->NumberOfIndices(),
                                                           0.0 );
        slhaParser.registerBlock( *(slhaBlockFinder->second) );
        slhaBlockPointers.push_back( slhaBlockFinder->second );
      }
      slhaFunctionoid->SetBlockPointer( slhaBlockFinder->second );
      returnPointer = slhaFunctionoid;
    }
    return returnPointer;
  }

  // This returns a pair of functionoid pointers which match the arguments
  // in commaSeparatedAliases.
  std::pair< ParameterFunctionoid*, ParameterFunctionoid* >
  RunningParameterManager::FindFunctionoidPair(
                                     std::string const& commaSeparatedAliases )
  {
    std::pair< ParameterFunctionoid*, ParameterFunctionoid* > returnPair( NULL,
                                                                        NULL );
    size_t commaPosition( commaSeparatedAliases.find_first_of( ',' ) );
    if( ( commaPosition > 0 )
        &&
        ( commaPosition < ( commaSeparatedAliases.size() - 2 ) ) )
    {
      returnPair.first
      = GetFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                                               commaSeparatedAliases.substr( 0,
                                                           commaPosition ) ) );
      if( returnPair.first != NULL )
      {
        returnPair.second
        = GetFunctionoid( BOL::StringParser::trimFromFrontAndBack(
                         commaSeparatedAliases.substr( commaPosition + 1 ) ) );
      }
    }
    return returnPair;
  }

  // This returns a string that should be valid Python to set the running
  // parameters in order for a given renormalization scale.
  std::string RunningParameterManager::RunningParametersAsPython() const
  {
    std::stringstream stringBuilder;
    stringBuilder << "# Running parameters, rp0, rp1, ...:\n";
    for( std::vector< ParameterFunctionoid* >::const_iterator
         runningParameter( parameterFunctionoidPointers.begin() );
         runningParameter < parameterFunctionoidPointers.end();
         ++runningParameter )
    {
      stringBuilder
      << (*runningParameter)->PythonParameterName() << " = 0.0\n";
    }
    stringBuilder << "\n"
    "# This function updates the running parameters for a new scale, through\n"
    "# the natural logarithm of the scale, lnQ:\n"
    "def UpdateRunningParameters( lnQ ):\n";
    for( std::vector< ParameterFunctionoid* >::const_iterator
         runningParameter( parameterFunctionoidPointers.begin() );
         runningParameter < parameterFunctionoidPointers.end();
         ++runningParameter )
    {
      stringBuilder
      << "    global " << (*runningParameter)->PythonParameterName() << "\n"
      << "    " << (*runningParameter)->PythonParameterEvaluation() << "\n";
    }
    stringBuilder << "\n"
    "# The running parameters are initially set for the lowest scale found\n"
    "# for the blocks in the SLHA file.\n"
    "UpdateRunningParameters( math.log( " << lowestBlockScale << " ) )\n"
    "invQSq = ( 1.0 / ( " << lowestBlockScale << " )**2 )\n"
    "# End of running parameters.\n";

    return stringBuilder.str();
  }

} /* namespace VevaciousPlusPlus */
