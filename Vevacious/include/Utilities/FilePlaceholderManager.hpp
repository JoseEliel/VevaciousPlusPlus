/*
 * FilePlaceholderManager.hpp
 *
 *  Created on: Dec 15, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef FILEPLACEHOLDERMANAGER_HPP_
#define FILEPLACEHOLDERMANAGER_HPP_

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <dirent.h>

namespace VevaciousPlusPlus
{
  // This struct is just a convenient grouping of data for the following class.
  struct FilenameTriple
  {
    std::string inputFile;
    std::string placeholderFile;
    std::string outputFile;
  };

  // This class manages triples of filenames with the intention that each
  //triple consists of a name for an input file, a name for a placeholder
  // file to indicate that the output file is currently being worked on, and a
  // name for an output file.
  class FilePlaceholderManager
  {
  public:
    FilePlaceholderManager( std::string const inputSuffix = "",
                            std::string const placeholderSuffix = "",
                            std::string const outputSuffix = "" ) :
                            inputSuffix( inputSuffix ),
                            placeholderSuffix( placeholderSuffix ),
                            outputSuffix( outputSuffix ),
                            filenameTriples(),
                            whichTriple(),
                            currentTriple(),
                            lastPlaceholder( "" ),
                            directoryPointer( NULL ),
                            structPointer( NULL ),
                            currentFilename( "" ) {}

    ~FilePlaceholderManager() { CloseDirectory(); }


    // This takes all the names of all the files in the directory given by
    // inputDirectory which end in inputSuffix, and places each in a
    // FilenameTriple with the placeholder filename given by replacing
    // inputSuffix with placeholderSuffix, replacing the directory path with
    // placeholderDirectory, and the same for the output file using
    // outputSuffix and outputDirectory.
    void PrepareFilenames( std::string const& inputDirectory,
                           std::string const& placeholderDirectory,
                           std::string const& outputDirectory );

    // This takes all the names of all the files in baseFilenames and places
    // each in a FilenameTriple with inputSuffix appended as the input
    // filename, with placeholderSuffix appended as the placeholder filename,
    // replacing the directory path with placeholderDirectory, and with
    // outputSuffix appended as the output file, replacing the directory path
    // with outputDirectory.
    void PrepareFilenames( std::vector< std::string > const& baseFilenames,
                           std::string const inputDirectory,
                           std::string const& placeholderDirectory,
                           std::string const& outputDirectory );

    // This looks to find the first FilenameTriple in filenameTriples which
    // has both a placeholder filename and an output filename which do not yet
    // exist in the file system, and returns true if there was such a triple.
    // If deleteLastPlaceholder is true, the previous placeholder is also
    // deleted.
    bool HoldNextPlace( bool const deleteLastPlaceholder = true );

    std::string const& CurrentInput() const { return whichTriple->inputFile; }

    std::string const& CurrentPlaceholder() const
    { return whichTriple->placeholderFile; }

    std::string const& CurrentOutput() const
    { return whichTriple->outputFile; }


  protected:
    std::string const inputSuffix;
    std::string const placeholderSuffix;
    std::string const outputSuffix;
    std::vector< FilenameTriple > filenameTriples;
    std::vector< FilenameTriple >::const_iterator whichTriple;
    FilenameTriple currentTriple;
    std::string lastPlaceholder;
    DIR* directoryPointer;
    // struct dirent* structPointer;
    dirent* structPointer;
    std::string currentFilename;


    // This tries to close the directory pointed to by directoryPointer if it
    // is not NULL, throwing an exception if this was not possible.
    inline void CloseDirectory();

    // This opens the given directory putting the result in directoryPointer,
    // throwing an exception if this was not possible.
    void OpenDirectory( std::string const& inputDirectory );

    // This reads the "filenames" in directoryPointer until it finds a name
    // that is not ".", "..", or too short to contain inputSuffix, or until
    // there is nothing left in directoryPointer. It returns true if it finds a
    // file, false if not.
    bool ReadNextFilenameInDirectory();

    // This tries to run mkdir with the -p argument, throwing an exception if
    // it could not.
    void EnsureDirectoryExists( std::string const& directoryName );

    // This deletes the file with the given name, throwing an exception if it
    // cannot.
    void DeleteFile( std::string const& fileName );

    // This tries to run systemCommand, throwing an exception if it cannot.
    void RunSystemCommand( std::string const& systemCommand );

    // This tries to find the next input file without its corresponding output
    // file or placeholder, and returns true if it found such a file.
    bool FindNextPlace();

    // This returns true if a file with the name fileName exists, determined by
    // trying to open it.
    bool FileExists( std::string const& fileName );
  };





  // This takes all the names of all the files in the directory given by
  // inputDirectory which end in inputSuffix, and places each in a
  // FilenameTriple with the placeholder filename given by replacing
  // inputSuffix with placeholderSuffix, replacing the directory path with
  // placeholderDirectory, and the same for the output file using
  // outputSuffix and outputDirectory.
  inline void
  FilePlaceholderManager::PrepareFilenames( std::string const& inputDirectory,
                                       std::string const& placeholderDirectory,
                                           std::string const& outputDirectory )
  {
    std::vector< std::string > validBaseFilenames;
    OpenDirectory( inputDirectory );
    while( ReadNextFilenameInDirectory() )
    {
      if( currentFilename.compare( ( currentFilename.size()
                                     - inputSuffix.size() ),
                                   inputSuffix.size(),
                                   inputSuffix ) == 0 )
      {
        validBaseFilenames.push_back( currentFilename.substr( 0,
                           ( currentFilename.size() - inputSuffix.size() ) ) );
      }
    }
    PrepareFilenames( validBaseFilenames,
                      inputDirectory,
                      placeholderDirectory,
                      outputDirectory );
    CloseDirectory();
  }

  // This takes all the names of all the files in baseFilenames and places each
  // in a FilenameTriple with inputSuffix appended as the input filename, with
  // placeholderSuffix appended as the placeholder filename, replacing the
  // directory path with placeholderDirectory, and with outputSuffix appended
  // as the output file, replacing the directory path with outputDirectory.
  inline void FilePlaceholderManager::PrepareFilenames(
                               std::vector< std::string > const& baseFilenames,
                                              std::string const inputDirectory,
                                       std::string const& placeholderDirectory,
                                           std::string const& outputDirectory )
  {
    filenameTriples.clear();
    EnsureDirectoryExists( placeholderDirectory );
    EnsureDirectoryExists( outputDirectory );
    for( std::vector< std::string >::const_iterator
         baseFilename( baseFilenames.begin() );
         baseFilenames.end() > baseFilename;
         ++baseFilename )
    {
      currentTriple.inputFile.assign( inputDirectory + "/"
                                      + (*baseFilename) + inputSuffix );
      currentTriple.placeholderFile.assign( placeholderDirectory + "/"
                                  + (*baseFilename) + placeholderSuffix );
      currentTriple.outputFile.assign( outputDirectory + "/"
                                       + (*baseFilename) + outputSuffix );
      filenameTriples.push_back( currentTriple );
    }
    whichTriple = filenameTriples.begin();
  }

  // This looks to find the first FilenameTriple in filenameTriples which
  // has both a placeholder filename and an output filename which do not yet
  // exist in the file system, and returns true if there was such a triple.
  // If deleteLastPlaceholder is true, the previous placeholder is also
  // deleted.
  inline bool
  FilePlaceholderManager::HoldNextPlace( bool const deleteLastPlaceholder )
  {
    bool const foundNextPlace( FindNextPlace() );
    if( deleteLastPlaceholder )
    {
      DeleteFile( lastPlaceholder );
    }
    if( foundNextPlace )
    {
      lastPlaceholder.assign( whichTriple->placeholderFile );
      std::ofstream placeholderStream( lastPlaceholder.c_str() );
      placeholderStream << "This is a placeholder!" << std::endl;
      placeholderStream.close();
      return true;
    }
    else
    {
      return false;
    }
  }

  // This tries to close the directory pointed to by directoryPointer if it is
  // not NULL, throwing an exception if this was not possible.
  inline void FilePlaceholderManager::CloseDirectory()
  {
    if( ( directoryPointer != NULL )
        &&
        ( closedir( directoryPointer ) == -1 ) )
    {
      throw std::runtime_error( "Could not close directory!" );
    }
    directoryPointer = NULL;
  }

  // This opens the given directory putting the result in directoryPointer,
  // throwing an exception if this was not possible.
  inline void
  FilePlaceholderManager::OpenDirectory( std::string const& inputDirectory )
  {
    directoryPointer = opendir( inputDirectory.c_str() );
    if( directoryPointer == NULL )
    {
      throw std::runtime_error( "Could not read directory!" );
    }
  }

  // This reads the "filenames" in directoryPointer until it finds a name that
  // is not ".", "..", or too short to contain inputSuffix, or until there is
  // nothing left in directoryPointer. It returns true if it finds a file,
  // false if not.
  inline bool FilePlaceholderManager::ReadNextFilenameInDirectory()
  {
    structPointer = readdir( directoryPointer );
    if( structPointer == NULL )
    {
      return false;
    }
    currentFilename.assign( structPointer->d_name );
    if( ( currentFilename == "." )
        ||
        ( currentFilename == ".." )
        ||
        !( inputSuffix.size() < currentFilename.size() ) )
    {
      return ReadNextFilenameInDirectory();
    }
    return true;
  }

  // This tries to run mkdir with the -p argument, throwing an exception if it
  // could not.
  inline void FilePlaceholderManager::EnsureDirectoryExists(
                                             std::string const& directoryName )
  {
    std::stringstream commandBuilder;
    commandBuilder << "mkdir -p " << directoryName;
    RunSystemCommand( commandBuilder.str() );
  }

  // This tries to run systemCommand, throwing an exception if it cannot.
  inline void
  FilePlaceholderManager::RunSystemCommand( std::string const& systemCommand )
  {
    if( system( systemCommand.c_str() ) == -1 )
    {
      std::stringstream errorBuilder;
      errorBuilder
      << "Could not execute \"" << systemCommand << "\"!";
      throw std::runtime_error( errorBuilder.str() );
    }
  }

  // This deletes the file with the given name, throwing an exception if it
  // cannot.
  inline void FilePlaceholderManager::DeleteFile( std::string const& fileName )
  {
    std::stringstream commandBuilder;
    commandBuilder << "rm " << fileName;
    RunSystemCommand( commandBuilder.str() );
  }

  // This tries to find the next input file without its corresponding output
  // file or placeholder, and returns true if it found such a file.
  inline bool FilePlaceholderManager::FindNextPlace()
  {
    while( ( whichTriple != filenameTriples.end() )
           &&
           ( FileExists( whichTriple->outputFile )
             ||
             FileExists( whichTriple->placeholderFile )
             ||
             !(FileExists( whichTriple->inputFile )) ) )
    {
      ++whichTriple;
    }
    return ( whichTriple != filenameTriples.end() );
  }

  // This returns true if a file with the name fileName exists, determined by
  // trying to open it.
  inline bool
  FilePlaceholderManager::FileExists( std::string const& fileName )
  {
    std::ifstream fileStream( fileName.c_str() );
    if( fileStream.good() )
    {
      fileStream.close();
      return true;
    }
    else
    {
      return false;
    }
  }

} /* namespace BOL */
#endif /* FILEPLACEHOLDERMANAGER_HPP_ */
