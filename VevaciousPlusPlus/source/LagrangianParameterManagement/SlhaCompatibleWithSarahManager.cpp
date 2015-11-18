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
                                      std::string const& fixedScaleArgument ) :
    SlhaBlocksWithSpecialCasesManager( validBlocksString,
                                       minimumScaleType,
                                       minimumScaleArgument,
                                       fixedScaleType,
                                       fixedScaleArgument )
  {
    InitializeSarahAliases();
  }

  SlhaCompatibleWithSarahManager::SlhaCompatibleWithSarahManager(
                                 std::set< std::string > const& validBlocksSet,
                                           std::string const& minimumScaleType,
                                       std::string const& minimumScaleArgument,
                                             std::string const& fixedScaleType,
                                      std::string const& fixedScaleArgument ) :
    SlhaBlocksWithSpecialCasesManager( validBlocks,
                                       minimumScaleType,
                                       minimumScaleArgument,
                                       fixedScaleType,
                                       fixedScaleArgument )
  {
    InitializeSarahAliases();
  }

  SlhaCompatibleWithSarahManager::~SlhaCompatibleWithSarahManager()
  {
    // This does nothing beyond what the base SlhaBlocksWithSpecialCasesManager
    // destructor does.
  }


  // This adds all the valid aliases to aliasesToSwitchStrings.
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
    MapCaseStringAndSlhaBlockToCaseString( "BmuTree",
                                           "TREEHMIX[101]" );
    MapCaseStringAndSlhaBlockToCaseString( "BmuLoop",
                                           "LOOPHMIX[101]" );
    MapCaseStringAndSlhaBlockToCaseString( "mHdSqTree",
                                           "TREEMSOFT[ 21 ]" );
    MapCaseStringAndSlhaBlockToCaseString( "mHdSqLoop",
                                           "LOOPMSOFT[ 21 ]" );
    MapCaseStringAndSlhaBlockToCaseString( "mHuSqTree",
                                           "TREEMSOFT[ 22 ]" );
    MapCaseStringAndSlhaBlockToCaseString( "mHuSqLoop",
                                           "LOOPMSOFT[ 22 ]" );
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
      SlhaTwoSourceFunctionoid temporaryParameter( 0,
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

      SlhaTwoSourceFunctionoid temporaryParameter( 0,
                                                   0,
                                                   0 );
      return
      temporaryParameter( OnceOffBlockEntry( ( treeOrLoop + baseCase ),
                                             logarithmOfScale ),
                          OnceOffParameter( baseCase,
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
    // debugging:
    /**/std::cout << std::endl << "debugging:"
    << std::endl
    << "SlhaCompatibleWithSarahManager::RegisterUnregisteredSpecialCase( \""
    << caseString << "\" ) called.";
    std::cout << std::endl;/**/

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
      return AddNewDerivedParameter( caseString,
                                     new SlhaTwoSourceFunctionoid(
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

      // debugging:
      /**/std::cout << std::endl << "debugging:"
      << std::endl
      << "caseString = \"" << caseString << "\", baseCase = \"" << baseCase
      << "\", treeOrLoop + baseCase = \"" << ( treeOrLoop + baseCase )
      << "\", sarahBlockEntry = \"" << sarahBlockEntry << "\"";
      std::cout << std::endl;/**/

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
      return AddNewDerivedParameter( caseString,
                                     new SlhaTwoSourceFunctionoid(
                                              numberOfDistinctActiveParameters,
                                                                   sarahIndex,
                                                             pureSlhaIndex ) );
    }
    else
    {
      return
      SlhaBlocksWithSpecialCasesManager::RegisterUnregisteredSpecialCase(
                                                                  caseString );
    }
  }

} /* namespace VevaciousPlusPlus */
