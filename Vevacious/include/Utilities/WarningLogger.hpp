/*
 * WarningLogger.hpp
 *
 *  Created on: Dec 17, 2015
 *      Author: bol
 */

#ifndef WARNINGLOGGER_HPP_
#define WARNINGLOGGER_HPP_

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <cstddef>
#include <iomanip>

namespace VevaciousPlusPlus
{

  class WarningLogger
  {
  public:
    static void
    SetWarningRecord( std::vector< std::string >* const warningDestination );

    // This prints the warning to std::cout and also stores it for later
    // recall.
    static void LogWarning( std::string const& warningMessage );


  private:
    static std::vector< std::string >* warningMessages;
  };




  inline void WarningLogger::SetWarningRecord(
                         std::vector< std::string >* const warningDestination )
  {
    warningMessages = warningDestination;
  }

  // This prints the warning to std::cout and also stores it for later
  // recall.
  inline void WarningLogger::LogWarning( std::string const& warningMessage )
  {
    if( warningMessages != NULL )
    {
      warningMessages->push_back( warningMessage );
    }
    std::cout << "Warning: " << warningMessage << std::endl;
  }

}

#endif /* WARNINGLOGGER_HPP_ */
