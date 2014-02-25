/*
 * VevaciousPlusPlusMain.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */


#include "../include/StandardIncludes.hpp"
#include "../include/VevaciousPlusPlus.hpp"

int main( int argumentCount,
          char** argumentCharArrays )
{
  // This takes care of the command-line arguments.
  BOL::ArgumentParser argumentParser( argumentCount,
                                      argumentCharArrays,
                                      "init",
                                      "VevaciousPlusPlusInitialization.xml" );

  // To create the VevaciousPlusPlus object, one must know what kind of
  // strategies one wants to use to minimize the potential and to calculate
  // tunneling between vacua. These strategies are chosen by passing instances
  // of classes derived from abstract base classes to the VevaciousPlusPlus
  // constructor. Since constructing such instances requires one to provide the
  // code enabling them to perform their tasks, one has to start with an
  // instance of a class derived from the PotentialFunction abstract base class
  // and then pass that to the classes derived from the PotentialMinimizer and
  // TunnelingCalculator abstract base classes, unless one wants to out in some
  // hacking effort.

  // The default strategies are:
  // - minimizing the potential by a combination of homotopy continuation to
  //   find starting points with gradient minimization, as implemented in the
  //   HomotopyContinuationAndGradient class, which requires an instance of a
  //   class derived from the HomotopyContinuationReadyPotential abstract base
  //   class;
  // - calculating tunneling by evaluating the bounce action along paths
  //   between the vacua parameterized by splines, as implemented in the
  //   BounceWithSplines class, which requires an instance of a class derived
  //   from the PotentialFunction abstract base clases.
  // Hence creating an instance of the MassCorrectedPotential class, which
  // inherits from both HomotopyContinuationReadyPotential and
  // PotentialFunction (through HomotopyContinuationReadyPotential itself
  // already being derived from PotentialFunction), forms the basis of the
  // components of the VevaciousPlusPlus object.

  // The MassCorrectedPotential constructor takes a string with the name of the
  // model file (including the path):
  VevaciousPlusPlus::MassCorrectedPotential
  potentialFunction( argumentParser.fromTag( "model",
                                             "./ModelFiles/SM.vin" ) );

  // Now the HomotopyContinuationAndGradient object and the BounceWithSplines
  // object can be constructed:
  VevaciousPlusPlus::HomotopyContinuationAndGradient
  potentialMinimizer( potentialFunction );
  VevaciousPlusPlus::BounceWithSplines
  tunnelingCalculator( potentialFunction );

  // Create the VevaciousTwo object, telling it where to find its settings,
  // such as the model to use and the way to calculate the tunneling time.
  VevaciousPlusPlus::VevaciousPlusPlus vevaciousPlusPlus( argumentParser,
                                                          potentialMinimizer,
                                                         tunnelingCalculator );

  // Solve a parameter point:
  std::string slhaFile( argumentParser.fromTag( "slha",
                                                "slha.out" ) );
  vevaciousPlusPlus.RunPoint( slhaFile );

  // Write the results:
  vevaciousPlusPlus.WriteXmlResults( argumentParser.fromTag( "output",
                                                    ( slhaFile + ".vout" ) ) );
  vevaciousPlusPlus.WriteSlhaResults( slhaFile );

  // this was a triumph! I'm making a note here:
  return EXIT_SUCCESS;
}






