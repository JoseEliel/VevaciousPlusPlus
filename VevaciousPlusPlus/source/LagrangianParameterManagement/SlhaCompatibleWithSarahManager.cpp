/*
 * SlhaCompatibleWithSarahManager.cpp
 *
 *  Created on: Nov 3, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LagrangianParameterManagement/SlhaCompatibleWithSarahManager.hpp"

namespace VevaciousPlusPlus
{

  SlhaCompatibleWithSarahManager::SlhaCompatibleWithSarahManager(
                                          std::string const& validBlocksString,
                                           std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                             std::string const& fixedScaleType,
                                         std::string const& fixedScaleArgument,
                                           std::string const& maximumScaleType,
                                    std::string const& maximumScaleArgument ) :
    SlhaBlocksWithSpecialCasesManager( validBlocksString,
                                       minimumScaleType,
                                       minimumScaleArgument,
                                       fixedScaleType,
                                       fixedScaleArgument,
                                       maximumScaleType,
                                       maximumScaleArgument )
  {
	RegisterDerivedParameters(derivedparameters); //This is not optimal since the derived parameters shouldn't overlap with the SarahAliases
    InitializeSarahAliases();
  }

  SlhaCompatibleWithSarahManager::SlhaCompatibleWithSarahManager(
                                 std::set< std::string > const& validBlocksSet,
                                           std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                             std::string const& fixedScaleType,
                                         std::string const& fixedScaleArgument,
                                           std::string const& maximumScaleType,
                                    std::string const& maximumScaleArgument ) :
    SlhaBlocksWithSpecialCasesManager( validBlocksSet,
                                       minimumScaleType,
                                       minimumScaleArgument,
                                       fixedScaleType,
                                       fixedScaleArgument,
                                       maximumScaleType,
                                       maximumScaleArgument )
  {
	RegisterDerivedParameters(derivedparameters); //This is not optimal since the derived parameters shouldn't overlap with the SarahAliases
    InitializeSarahAliases();
  }

  SlhaCompatibleWithSarahManager::SlhaCompatibleWithSarahManager(
                                             std::string const& xmlFileName ) :
    SlhaBlocksWithSpecialCasesManager( xmlFileName )
  {
    RegisterDerivedParameters(derivedparameters); //This is not optimal since the derived parameters shouldn't overlap with the SarahAliases
    InitializeSarahAliases();
  }

  SlhaCompatibleWithSarahManager::~SlhaCompatibleWithSarahManager()
  {
    // This does nothing beyond what the base SlhaBlocksWithSpecialCasesManager
    // destructor does.
  }

   void SlhaCompatibleWithSarahManager::InitializeSarahAliases()
  {
    // Putting the simple block names as special cases for SARAH means that
    // a SARAH-generated model file which assumes the extra SARAH blocks will
    // still work with non-SARAH SLHA files, to the extent that the DRbar
    // values will be used as both the SARAH tree and loop values.
    MapCaseStringAndSlhaBlockToCaseString( "DsbVd",
                                           "HMIX[102]" );
    MapCaseStringAndSlhaBlockToCaseString( "DsbVu",
                                           "HMIX[103]" );
    // The constructor for the base SlhaBlocksWithSpecialCasesManager covers
    // adding the basic Bmu to the alias mapping, even though this derived
    // class over-writes what ends up as its functionoid.
    MapCaseStringAndSlhaBlockToCaseString( "Bmu",
                                           "HMIX[101]" );
    MapCaseStringAndSlhaBlockToCaseString( "muTree",
                                           "TREEHMIX[1]" );
    MapCaseStringAndSlhaBlockToCaseString( "muLoop",
                                           "LOOPHMIX[1]" );
    MapCaseStringAndSlhaBlockToCaseString( "muDelta",
                                           "DELTAHMIX[1]" );
    MapCaseStringAndSlhaBlockToCaseString( "BmuTree",
                                           "TREEHMIX[101]" );
    MapCaseStringAndSlhaBlockToCaseString( "BmuLoop",
                                           "LOOPHMIX[101]" );
    MapCaseStringAndSlhaBlockToCaseString( "BmuDelta",
                                           "DELTAHMIX[101]" );
    MapCaseStringAndSlhaBlockToCaseString( "mHdSqTree",
                                           "TREEMSOFT[ 21 ]" );
    MapCaseStringAndSlhaBlockToCaseString( "mHdSqLoop",
                                           "LOOPMSOFT[ 21 ]" );
    MapCaseStringAndSlhaBlockToCaseString( "mHdSqDelta",
                                           "DELTAMSOFT[ 21 ]" );
    MapCaseStringAndSlhaBlockToCaseString( "mHuSqTree",
                                           "TREEMSOFT[ 22 ]" );
    MapCaseStringAndSlhaBlockToCaseString( "mHuSqLoop",
                                           "LOOPMSOFT[ 22 ]" );
    MapCaseStringAndSlhaBlockToCaseString( "mHuSqDelta",
                                           "DELTAMSOFT[ 22 ]" );
  }

  // This duplicates a lot of code from RegisterUnregisteredSpecialCase, but
  // there doesn't seem to be an elegant way of using the common code as
  // there is too much entanglement with registering new parameters or not.
  double SlhaCompatibleWithSarahManager::OnceOffSpecialCase(
                                                 std::string const& caseString,
                                          double const logarithmOfScale ) const
  {
    if( ( caseString == "Bmu" )
        ||
        ( caseString == "DsbVd" )
        ||
        ( caseString == "DsbVu" ) )
    {
      // Assume that the special case is "Bmu" and then change if it is "DsbVd"
      // or "DsbVu".
      std::string sarahBlock( "HMIX[ 101 ]" );
      if( caseString == "DsbVd" )
      {
        sarahBlock = "HMIX[ 102 ]";
      }
      else if( caseString == "DsbVu" )
      {
        sarahBlock = "HMIX[ 103 ]";
      }
      LhaTwoSourceFunctionoid temporaryParameter( 0,
                                                  0,
                                                  0 );
      return temporaryParameter( OnceOffParameter( sarahBlock,
                                                   logarithmOfScale ),
             SlhaBlocksWithSpecialCasesManager::OnceOffSpecialCase( caseString,
                                                          logarithmOfScale ) );
    }
    else if( ( caseString == "muTree" )
             ||
             ( caseString == "muLoop" )
             ||
             ( caseString == "BmuTree" )
             ||
             ( caseString == "BmuLoop" )
             ||
             ( caseString == "mHdSqTree" )
             ||
             ( caseString == "mHdSqLoop" )
             ||
             ( caseString == "mHuSqTree" )
             ||
             ( caseString == "mHuSqLoop" ) )
    {
      std::string baseCase( caseString.substr( 0,
                                               caseString.size() - 4 ) );
      std::string const
      treeOrLoop( ( caseString[ caseString.size() - 1 ] == 'e' ) ?
                  "TREE" :
                  "LOOP" );
      std::string sarahBlockEntry( "error" );
      if( baseCase == "Bmu" )
      {
        sarahBlockEntry = FormatVariable( treeOrLoop + "HMIX[ 101 ]" );
      }
      else
      {
        if( baseCase == "mu" )
        {
          baseCase = "HMIX[ 1 ]";
        }
        else if( baseCase == "mHdSq" )
        {
          baseCase = "MSOFT[ 21 ]";
        }
        else if( baseCase == "mHuSq" )
        {
          baseCase = "MSOFT[ 22 ]";
        }
        sarahBlockEntry = FormatVariable( treeOrLoop + baseCase );
      }

      LhaTwoSourceFunctionoid temporaryParameter( 0,
                                                  0,
                                                  0 );
      return temporaryParameter( OnceOffBlockEntry( sarahBlockEntry,
                                                    logarithmOfScale ),
                                 OnceOffParameter( baseCase,
                                                   logarithmOfScale ) );
    }
    else if( ( caseString == "muDelta" )
             ||
             ( caseString == "BmuDelta" )
             ||
             ( caseString == "mHdSqDelta" )
             ||
             ( caseString == "mHuSqDelta" ) )
    {
      std::string baseCase( caseString.substr( 0,
                                               caseString.size() - 5 ) );
      std::string sarahBlockEntry( "error" );
      if( baseCase == "Bmu" )
      {
        sarahBlockEntry = FormatVariable( "DELTAHMIX[ 101 ]" );
      }
      else if( baseCase == "mu" )
      {
        sarahBlockEntry = FormatVariable( "DELTAHMIX[ 1 ]" );
      }
      else if( baseCase == "mHdSq" )
      {
        sarahBlockEntry = FormatVariable( "DELTAMSOFT[ 21 ]" );
      }
      else if( baseCase == "mHuSq" )
      {
        sarahBlockEntry = FormatVariable( "DELTAMSOFT[ 22 ]" );
      }

      double const deltaBlockEntry( OnceOffBlockEntry( sarahBlockEntry,
                                                       logarithmOfScale ) );
      if( deltaBlockEntry != 0.0 )
      {
        return deltaBlockEntry;
      }
      LhaDifferenceFunctionoid temporaryParameter( 0,
                                                   0,
                                                   0 );
      return temporaryParameter( OnceOffParameter( ( baseCase + "Tree" ),
                                                   logarithmOfScale ),
                                 OnceOffParameter( ( baseCase + "Loop" ),
                                                   logarithmOfScale ) );
    }
    else
    {
      return SlhaBlocksWithSpecialCasesManager::OnceOffSpecialCase( caseString,
                                                            logarithmOfScale );
    }
  }

  // This adds the parameter based on the alias given by switchString for the
  // parameter.
  std::pair< bool, size_t >
  SlhaCompatibleWithSarahManager::RegisterUnregisteredSpecialCase(
                                                std::string const& caseString )
  {
    if( ( caseString == "Bmu" )
        ||
        ( caseString == "DsbVd" )
        ||
        ( caseString == "DsbVu" ) )
    {
      // Assume that the special case is "Bmu" and then change if it is "DsbVd"
      // or "DsbVu".
      std::string sarahBlockEntry( "HMIX[ 101 ]" );
      if( caseString == "DsbVd" )
      {
        sarahBlockEntry = "HMIX[ 102 ]";
      }
      else if( caseString == "DsbVu" )
      {
        sarahBlockEntry = "HMIX[ 103 ]";
      }

      // Since what sarahBlockEntry contains is already an alias for caseString
      // (so as to have potential function files designed by SARAH work with
      // pure SLHA files, which may or may not be a good idea), we have to look
      // for the block entry functionoid directly, and create it if it doesn't
      // exist, rather than recursing through
      // LesHouchesAccordBlockEntryManager::RegisterParameter(sarahBlockEntry),
      // which will just come back here through
      // SlhaBlocksWithSpecialCasesManager::RegisterUnregisteredParameter
      // finding sarahBlockEntry as an alias for caseString, then calling the
      // derived override of RegisterUnregisteredSpecialCase, bringing us back
      // here.
      size_t const sarahIndex( RegisterBlockEntry( sarahBlockEntry ) );

      // Since the control flow got to here, the parameter to which caseString
      // refers is definitely not already registered, so it's not a problem to
      // directly register it here, yielding the pure SLHA special case.
      size_t const pureSlhaIndex(
            SlhaBlocksWithSpecialCasesManager::RegisterUnregisteredSpecialCase(
                                                         caseString ).second );

      // This will overwrite activeParametersToIndices[ caseString ] to map
      // caseString to numberOfDistinctActiveParameters as its index in the
      // values vector (corresponding now to the new SlhaTwoSourceFunctionoid),
      // rather than to the pure SLHA functionoid as was set by
      // RegisterParameter( caseString ) (or previous to that).
      return PairAddNewDerivedParameter( caseString,
                                         new LhaTwoSourceFunctionoid(
                                              numberOfDistinctActiveParameters,
                                                                    sarahIndex,
                                                             pureSlhaIndex ) );
    }
    else if( ( caseString == "muTree" )
             ||
             ( caseString == "muLoop" )
             ||
             ( caseString == "BmuTree" )
             ||
             ( caseString == "BmuLoop" )
             ||
             ( caseString == "mHdSqTree" )
             ||
             ( caseString == "mHdSqLoop" )
             ||
             ( caseString == "mHuSqTree" )
             ||
             ( caseString == "mHuSqLoop" ) )
    {
      std::string baseCase( caseString.substr( 0,
                                               caseString.size() - 4 ) );
      std::string const
      treeOrLoop( ( caseString[ caseString.size() - 1 ] == 'e' ) ?
                  "TREE" :
                  "LOOP" );
      std::string sarahBlockEntry( "error" );
      if( baseCase == "Bmu" )
      {
        sarahBlockEntry = FormatVariable( treeOrLoop + "HMIX[ 101 ]" );
      }
      else
      {
        if( baseCase == "mu" )
        {
          baseCase = "HMIX[ 1 ]";
        }
        else if( baseCase == "mHdSq" )
        {
          baseCase = "MSOFT[ 21 ]";
        }
        else if( baseCase == "mHuSq" )
        {
          baseCase = "MSOFT[ 22 ]";
        }
        sarahBlockEntry = FormatVariable( treeOrLoop + baseCase );
      }

      // Since what sarahBlockEntry contains might already be an alias for
      // caseString (so as to have potential function files designed by SARAH
      // work with pure SLHA files, which may or may not be a good idea), we
      // have to look for the block entry functionoid directly, and create it
      // if it doesn't exist, rather than recursing through
      // LesHouchesAccordBlockEntryManager::RegisterParameter(sarahBlockEntry),
      // which will just come back here through
      // SlhaBlocksWithSpecialCasesManager::RegisterUnregisteredParameter
      // finding sarahBlockEntry as an alias for caseString, then calling the
      // derived override of RegisterUnregisteredSpecialCase, bringing us back
      // here.
      size_t const sarahIndex( RegisterBlockEntry( sarahBlockEntry ) );
      size_t const pureSlhaIndex( RegisterParameter( baseCase ).second );

      // This will overwrite activeParametersToIndices[ caseString ] to map
      // caseString to numberOfDistinctActiveParameters as its index in the
      // values vector (corresponding now to the new SlhaTwoSourceFunctionoid),
      // rather than to the pure SLHA functionoid as was set by
      // RegisterParameter( caseString ) (or previous to that).
      return PairAddNewDerivedParameter( caseString,
                                         new LhaTwoSourceFunctionoid(
                                              numberOfDistinctActiveParameters,
                                                                    sarahIndex,
                                                             pureSlhaIndex ) );
    }
    else if( ( caseString == "muDelta" )
             ||
             ( caseString == "BmuDelta" )
             ||
             ( caseString == "mHdSqDelta" )
             ||
             ( caseString == "mHuSqDelta" ) )
    {
      std::string baseCase( caseString.substr( 0,
                                               caseString.size() - 5 ) );
      std::string deltaBlockEntry( "error" );
      if( baseCase == "Bmu" )
      {
        deltaBlockEntry = FormatVariable( "DELTAHMIX[ 101 ]" );
      }
      else if( baseCase == "mu" )
      {
        deltaBlockEntry = FormatVariable( "DELTAHMIX[ 1 ]" );
      }
      else if( baseCase == "mHdSq" )
      {
        deltaBlockEntry = FormatVariable( "DELTAMSOFT[ 21 ]" );
      }
      else if( baseCase == "mHuSq" )
      {
        deltaBlockEntry = FormatVariable( "DELTAMSOFT[ 22 ]" );
      }

      // Since what sarahBlockEntry contains might already be an alias for
      // caseString (so as to have potential function files designed by SARAH
      // work with pure SLHA files, which may or may not be a good idea), we
      // have to look for the block entry functionoid directly, and create it
      // if it doesn't exist, rather than recursing through
      // LesHouchesAccordBlockEntryManager::RegisterParameter(sarahBlockEntry),
      // which will just come back here through
      // SlhaBlocksWithSpecialCasesManager::RegisterUnregisteredParameter
      // finding sarahBlockEntry as an alias for caseString, then calling the
      // derived override of RegisterUnregisteredSpecialCase, bringing us back
      // here.
      size_t const deltaIndex( RegisterBlockEntry( deltaBlockEntry ) );
      size_t const
      treeIndex( RegisterParameter( ( baseCase + "Tree" ) ).second );
      size_t const
      loopIndex( RegisterParameter( ( baseCase + "Loop" ) ).second );
      size_t const differenceIndex( AddNewDerivedParameter( caseString,
                                                  new LhaDifferenceFunctionoid(
                                              numberOfDistinctActiveParameters,
                                                                     loopIndex,
                                                               treeIndex ) ) );

      // This will overwrite activeParametersToIndices[ caseString ] to map
      // caseString to numberOfDistinctActiveParameters as its index in the
      // values vector (corresponding now to the new SlhaTwoSourceFunctionoid),
      // rather than to the pure SLHA functionoid as was set by
      // RegisterParameter( caseString ) (or previous to that).
      return PairAddNewDerivedParameter( caseString,
                                         new LhaTwoSourceFunctionoid(
                                              numberOfDistinctActiveParameters,
                                                                    deltaIndex,
                                                           differenceIndex ) );
    }
    else
    {
      return
      SlhaBlocksWithSpecialCasesManager::RegisterUnregisteredSpecialCase(
                                                                  caseString );
    }
  }

  //Registers a new derived Parameter given by the arguments in the Modelfile.
  //Currently only supports IFNONZERO[ ... ] and SLHAblock parameters
  void SlhaCompatibleWithSarahManager::RegisterDerivedParameters(std::vector<std::pair<std::string,std::string>> derivedparameters)
    {
    for( auto it = derivedparameters.begin(); it != derivedparameters.end(); it++) //first get all the parameter definitions
	{
        if(((*it).second).find("IFNONZERO")==std::string::npos)
            {
           LesHouchesAccordBlockEntryManager::RegisterNewParameter(LesHouchesAccordBlockEntryManager::CreateNewBlockEntry((*it).second),(*it).first);
            }
	}
    for( auto it = derivedparameters.begin(); it != derivedparameters.end(); it++) //TwoSource
      {
        if(((*it).second).find("IFNONZERO")!=std::string::npos)
            {
			std::string parnames(it->second);
			parnames = parnames.substr(parnames.find('[')+1, (parnames.find(']') - parnames.find('[') - 1));
			std::pair< bool, size_t > ifnonzero1 = LesHouchesAccordBlockEntryManager::RegisterParameter(parnames.substr(0,parnames.find(',')));
            std::pair< bool, size_t > ifnonzero2 = LesHouchesAccordBlockEntryManager::RegisterParameter(parnames.substr(parnames.find(',') + 1));
            if(ifnonzero1.first && ifnonzero2.first)
            {
                AddNewDerivedParameter((*it).first, new LhaTwoSourceFunctionoid(numberOfDistinctActiveParameters,ifnonzero1.second,ifnonzero2.second));

            }
            else
            {
                std::stringstream errorBuilder;
                errorBuilder
                << "RegisterDerivedParameter does not recognize the given IFNONZERO parameters, please check if the given parameters are correctly defined." << std::endl;
                throw std::runtime_error( errorBuilder.str() );
            }
            }
      }
	}

} /* namespace VevaciousPlusPlus */
