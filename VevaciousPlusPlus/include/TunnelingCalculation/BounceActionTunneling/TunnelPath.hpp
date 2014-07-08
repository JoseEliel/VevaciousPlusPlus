/*
 * TunnelPath.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef TUNNELPATH_HPP_
#define TUNNELPATH_HPP_

namespace VevaciousPlusPlus
{

  class TunnelPath
  {
  public:
    TunnelPath( size_t const numberOfFields,
                double const temperatureValue = 0.0 );
    virtual ~TunnelPath();

    size_t NumberOfFields() const{ return numberOfFields; }

    // This should fill fieldConfiguration with the values that the fields
    // should have when the path auxiliary is given by auxiliaryValue.
    void PutOnPathAt( std::vector< double >& fieldConfiguration,
                      double const auxiliaryValue ) const = 0;

    // This should return the dot product with itself of the derivative of the
    // field vector with respect to the path auxiliary evaluated at
    // auxiliaryValue.
    double SlopeSquared( double const auxiliaryValue ) const = 0;

    // This should return the dot product of the first derivative of the field
    // vector with the second derivative, both with respect to the path
    // auxiliary, evaluated at auxiliaryValue.
    double SlopeDotAcceleration( double const auxiliaryValue ) const = 0;

    double TemperatureValue() const{ return temperatureValue; }

    bool NonZeroTemperature() const{ return nonZeroTemperature; }

    // This is for debugging.
    std::string AsDebuggingString() const = 0;

    // This is for debugging.
    std::string FieldsString( double const auxiliaryValue ) const;


  protected:
    size_t const numberOfFields;

    void SetTemperature( double const temperatureValue )
    { this->temperatureValue = temperatureValue;
      nonZeroTemperature = ( temperatureValue > 0.0 ); }


  private:
    double temperatureValue;
    bool nonZeroTemperature;
  };




  // This is for debugging.
  inline std::string
  TunnelPath::FieldsString( double const auxiliaryValue ) const
  {
    std::vector< double > fieldConfiguration( numberOfFields );
    PutOnPathAt( fieldConfiguration,
                 auxiliaryValue );
    std::stringstream returnStream;
    for( size_t fieldIndex( 0 );
         fieldIndex < numberOfFields;
         ++fieldIndex )
    {
      if( fieldIndex > 0 )
      {
        returnStream << ", ";
      }
      returnStream
      << "f[" << fieldIndex << "] = " << fieldConfiguration[ fieldIndex ];
    }
    return returnStream.str();
  }

} /* namespace VevaciousPlusPlus */
#endif /* TUNNELPATH_HPP_ */
