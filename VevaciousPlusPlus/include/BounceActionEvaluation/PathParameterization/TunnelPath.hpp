/*
 * TunnelPath.hpp
 *
 *  Created on: Jul 1, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef TUNNELPATH_HPP_
#define TUNNELPATH_HPP_

#include "CommonIncludes.hpp"
#include "Eigen/Dense"

namespace VevaciousPlusPlus
{

  class TunnelPath
  {
  public:
    TunnelPath( size_t const numberOfFields,
                std::vector< double > const& pathParameterization,
                double const pathTemperature = 0.0 );
    virtual ~TunnelPath();


    std::vector< double >&
    PathParameterization(){ return pathParameterization; }

    std::vector< double > const&
    PathParameterization() const{ return pathParameterization; }

    size_t NumberOfFields() const{ return numberOfFields; }

    // This should fill fieldConfiguration with the values that the fields
    // should have when the path auxiliary is given by auxiliaryValue.
    virtual void PutOnPathAt( std::vector< double >& fieldConfiguration,
                              double const auxiliaryValue ) const = 0;

    // This should fill fieldConfiguration with the values that the fields
    // should have when the path auxiliary is given by auxiliaryValue.
    virtual void PutOnPathAt( Eigen::VectorXd& fieldConfiguration,
                              double const auxiliaryValue ) const = 0;

    // This should return the dot product with itself of the derivative of the
    // field vector with respect to the path auxiliary evaluated at
    // auxiliaryValue.
    virtual double SlopeSquared( double const auxiliaryValue ) const = 0;

    // This should return the dot product of the first derivative of the field
    // vector with the second derivative, both with respect to the path
    // auxiliary, evaluated at auxiliaryValue.
    virtual double
    SlopeDotAcceleration( double const auxiliaryValue ) const = 0;

    double TemperatureValue() const{ return pathTemperature; }

    bool NonZeroTemperature() const{ return nonZeroTemperature; }

    // This is for debugging.
    virtual std::string AsDebuggingString() const = 0;

    // This is for debugging.
    std::string FieldsString( double const auxiliaryValue ) const;


  protected:
    size_t const numberOfFields;
    std::vector< double > pathParameterization;

    void SetTemperature( double const temperatureValue )
    { this->pathTemperature = temperatureValue;
      nonZeroTemperature = ( temperatureValue > 0.0 ); }


  private:
    double pathTemperature;
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
