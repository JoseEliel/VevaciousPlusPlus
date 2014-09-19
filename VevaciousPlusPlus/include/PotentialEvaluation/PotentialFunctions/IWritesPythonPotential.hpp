/*
 * IWritesPythonPotential.hpp
 *
 *  Created on: Apr 16, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef PYTHONCONVERTIBLEPOTENTIAL_HPP_
#define PYTHONCONVERTIBLEPOTENTIAL_HPP_

#include "CommonIncludes.hpp"

namespace VevaciousPlusPlus
{

  // This is just an abstract base class that is performing the role of an
  // interface.
  class IWritesPythonPotential
  {
  public:
    IWritesPythonPotential(){}
    virtual ~IWritesPythonPotential(){}


    // This should write the potential as
    // def PotentialFunction( fv ): return ...
    // in pythonFilename for fv being an array of floating-point numbers in the
    // same order as they are for the field configurations as internal to this
    // C++ code.
    virtual void WriteAsPython( std::string const pythonFilename ) const = 0;
  };

} /* namespace VevaciousPlusPlus */
#endif /* PYTHONCONVERTIBLEPOTENTIAL_HPP_ */
