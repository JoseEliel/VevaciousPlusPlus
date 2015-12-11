/*
 * ThermalFunctions.hpp
 *
 *  Created on: Mar 11, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef THERMALFUNCTIONS_HPP_
#define THERMALFUNCTIONS_HPP_

#include <string>
#include <sstream>

namespace VevaciousPlusPlus
{

  class ThermalFunctions
  {
  public:
    static double BosonicJ( double const squareRatio );
    static double FermionicJ( double const squareRatio );

    static std::string JFunctionsAsPython();

  private:
    // LOTS OF DANGER HERE! Unfortunately we cannot rely on the user having a
    // C++11-compliant compiler, so we are forced to use const arrays of
    // doubles, and have to know how big they are and what they represent.
    static double const bosonMinusOneToMinusTwelve[ 111 ];
    static double const bosonZeroToMinusOne[ 101 ];
    static double const bosonZeroToPlusOne[ 101 ];
    static double const bosonPlusOneToPlusOneHundred[ 100 ];

    // Since the fermion masses-squared are the eigenvalues of mass matrices
    // that could be diagonalized by unitary transformations on the left and
    // right, they must be positive semi-definite. However, numerical issues
    // could easily lead to small but negative eigenvalues of the mass-squared
    // matrix, so it's just easier to leave the interpolation in.
    static double const fermionMinusOneToMinusTwelve[ 111 ];
    static double const fermionZeroToMinusOne[ 101 ];
    static double const fermionZeroToPlusOne[ 101 ];
    static double const fermionPlusOneToPlusOneHundred[ 100 ];

    // -1 to -12 (element [0] is -1, [111] is -12), in steps of 0.1.
    static double BosonMinusOneToMinusTwelve( double const squareRatio );

    // 0 to -1 (element [0] is 0, [101] is -1), in steps of 0.01.
    static double BosonZeroToMinusOne( double const squareRatio );

    // 0 to 1 (element [0] is 0, [101] is 1), in steps of 0.01.
    static double BosonZeroToPlusOne( double const squareRatio );

    // 1 to 100 (element [0] is 1, [100] is 100), in steps of 1.0.
    static double BosonPlusOneToPlusOneHundred( double const squareRatio );

    // -1 to -12 (element [0] is -1, [111] is -12), in steps of 0.1.
    static double FermionMinusOneToMinusTwelve( double const squareRatio );

    // 0 to -1 (element [0] is 0, [101] is -1), in steps of 0.01.
    static double FermionZeroToMinusOne( double const squareRatio );

    // 0 to 1 (element [0] is 0, [101] is 1), in steps of 0.01.
    static double FermionZeroToPlusOne( double const squareRatio );

    // 1 to 100 (element [0] is 1, [100] is 100), in steps of 1.0.
    static double FermionPlusOneToPlusOneHundred( double const squareRatio );
  };




  inline double ThermalFunctions::BosonicJ( double const squareRatio )
  {
    if( squareRatio <= -12.0 )
    {
      return 0.0;
    }
    else if( squareRatio <= -1.0 )
    {
      return BosonMinusOneToMinusTwelve( squareRatio );
    }
    else if( squareRatio < 0.0 )
    {
      return BosonZeroToMinusOne( squareRatio );
    }
    else if( squareRatio < 1.0 )
    {
      return BosonZeroToPlusOne( squareRatio );
    }
    else if( squareRatio < 100.0 )
    {
      return BosonPlusOneToPlusOneHundred( squareRatio );
    }
    else
    {
      return 0.0;
    }
  }

  // -1 to -12 (element [0] is -1, [111] is -12), in steps of 0.1.
  inline double ThermalFunctions::FermionicJ( double const squareRatio )
  {
    if( squareRatio <= -12.0 )
    {
      return 0.0;
    }
    else if( squareRatio <= -1.0 )
    {
      return FermionMinusOneToMinusTwelve( squareRatio );
    }
    else if( squareRatio < 0.0 )
    {
      return FermionZeroToMinusOne( squareRatio );
    }
    else if( squareRatio < 1.0 )
    {
      return FermionZeroToPlusOne( squareRatio );
    }
    else if( squareRatio < 100.0 )
    {
      return FermionPlusOneToPlusOneHundred( squareRatio );
    }
    else
    {
      return 0.0;
    }
  }

  // -1 to -12 (element [0] is -1, [111] is -12), in steps of 0.1.
  inline double
  ThermalFunctions::BosonMinusOneToMinusTwelve( double const squareRatio )
  {
    double scaledRatio( -10.0 * ( squareRatio + 1.0 ) );
    size_t floorIndex( static_cast< size_t >( scaledRatio ) );
    return ( bosonMinusOneToMinusTwelve[ floorIndex ]
             + ( 0.1 * ( scaledRatio - floorIndex )
                     * ( bosonMinusOneToMinusTwelve[ floorIndex + 1 ]
                         - bosonMinusOneToMinusTwelve[ floorIndex ] ) ) );
  }

  // 0 to -1 (element [0] is 0, [101] is -1), in steps of 0.01.
  inline double
  ThermalFunctions::BosonZeroToMinusOne( double const squareRatio )
  {
    double scaledRatio( -100.0 * squareRatio );
    size_t floorIndex( static_cast< size_t >( scaledRatio ) );
    return ( bosonZeroToMinusOne[ floorIndex ]
             + ( 0.01 * ( scaledRatio - floorIndex )
                      * ( bosonZeroToMinusOne[ floorIndex + 1 ]
                          - bosonZeroToMinusOne[ floorIndex ] ) ) );
  }

  // 0 to 1 (element [0] is 0, [101] is 1), in steps of 0.01.
  inline double
  ThermalFunctions::BosonZeroToPlusOne( double const squareRatio )
  {
    double scaledRatio( 100.0 * squareRatio );
    size_t floorIndex( static_cast< size_t >( scaledRatio ) );
    return ( bosonZeroToPlusOne[ floorIndex ]
             + ( 0.01 * ( scaledRatio - floorIndex )
                      * ( bosonZeroToPlusOne[ floorIndex + 1 ]
                          - bosonZeroToPlusOne[ floorIndex ] ) ) );
  }

  // 1 to 100 (element [0] is 1, [100] is 100), in steps of 1.0.
  inline double
  ThermalFunctions::BosonPlusOneToPlusOneHundred( double const squareRatio )
  {
    double scaledRatio( squareRatio - 1.0 );
    size_t floorIndex( static_cast< size_t >( scaledRatio ) );
    return ( bosonPlusOneToPlusOneHundred[ floorIndex ]
             + ( ( scaledRatio - floorIndex )
                 * ( bosonPlusOneToPlusOneHundred[ floorIndex + 1 ]
                     - bosonPlusOneToPlusOneHundred[ floorIndex ] ) ) );
  }

  // -1 to -12 (element [0] is -1, [111] is -12), in steps of 0.1.
  inline double
  ThermalFunctions::FermionMinusOneToMinusTwelve( double const squareRatio )
  {
    double scaledRatio( -10.0 * ( squareRatio + 1.0 ) );
    size_t floorIndex( static_cast< size_t >( scaledRatio ) );
    return ( fermionMinusOneToMinusTwelve[ floorIndex ]
             + ( 0.1 * ( scaledRatio - floorIndex )
                     * ( fermionMinusOneToMinusTwelve[ floorIndex + 1 ]
                         - fermionMinusOneToMinusTwelve[ floorIndex ] ) ) );
  }

  // 0 to -1 (element [0] is 0, [101] is -1), in steps of 0.01.
  inline double
  ThermalFunctions::FermionZeroToMinusOne( double const squareRatio )
  {
    double scaledRatio( -100.0 * squareRatio );
    size_t floorIndex( static_cast< size_t >( scaledRatio ) );
    return ( fermionZeroToMinusOne[ floorIndex ]
             + ( 0.01 * ( scaledRatio - floorIndex )
                      * ( fermionZeroToMinusOne[ floorIndex + 1 ]
                          - fermionZeroToMinusOne[ floorIndex ] ) ) );
  }

  // 0 to 1 (element [0] is 0, [101] is 1), in steps of 0.01.
  inline double
  ThermalFunctions::FermionZeroToPlusOne( double const squareRatio )
  {
    double scaledRatio( 100.0 * squareRatio );
    size_t floorIndex( static_cast< size_t >( scaledRatio ) );
    return ( fermionZeroToPlusOne[ floorIndex ]
             + ( 0.01 * ( scaledRatio - floorIndex )
                      * ( fermionZeroToPlusOne[ floorIndex + 1 ]
                          - fermionZeroToPlusOne[ floorIndex ] ) ) );
  }

  // 1 to 100 (element [0] is 1, [100] is 100), in steps of 1.0.
  inline double
  ThermalFunctions::FermionPlusOneToPlusOneHundred( double const squareRatio )
  {
    double scaledRatio( squareRatio - 1.0 );
    size_t floorIndex( static_cast< size_t >( scaledRatio ) );
    return ( fermionPlusOneToPlusOneHundred[ floorIndex ]
             + ( ( scaledRatio - floorIndex )
                 * ( fermionPlusOneToPlusOneHundred[ floorIndex + 1 ]
                     - fermionPlusOneToPlusOneHundred[ floorIndex ] ) ) );
  }

} /* namespace VevaciousPlusPlus */
#endif /* THERMALFUNCTIONS_HPP_ */
