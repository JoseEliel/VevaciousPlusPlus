/*
 * VersionInformation.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef VERSIONINFORMATION_HPP_
#define VERSIONINFORMATION_HPP_

#include <string>

namespace VevaciousPlusPlus
{

  class VersionInformation
  {
  public:
    static std::string CurrentVersion() { return "1.0.00.alpha001"; }
    static std::string CurrentCitation() { return "[none as yet]"; }
  };

} /* namespace VevaciousPlusPlus */
#endif /* VERSIONINFORMATION_HPP_ */
