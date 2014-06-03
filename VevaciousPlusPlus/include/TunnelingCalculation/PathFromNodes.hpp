/*
 * PathFromNodes.hpp
 *
 *  Created on: Jun 2, 2014
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef PATHFROMNODES_HPP_
#define PATHFROMNODES_HPP_

namespace VevaciousPlusPlus
{

  class PathFromNodes
  {
  public:
    PathFromNodes( size_t const numberOfFields,
                   size_t const referenceFieldIndex,
                   size_t const numberOfVaryingPathNodes );
    virtual
    ~PathFromNodes();


    // This turns a flattened matrix of numbers parameterizing the path from
    // the false vacuum to the true vacuum through field space. It assumes that
    // numberOfVaryingPathNodes nodes of numberOfParameterizationFields field
    // values (in the plane where the reference field is 0) are given by
    // pathParameterization, and sets fieldsAsPolynomials and
    // fieldDerivativesAsPolynomials appropriately.
    void operator()( std::vector< double > const& pathParameterization,
                     std::vector< double > const& straightPath,
                     double const straightPathInverseLengthSquared,
                     std::vector< double > const& falseVacuumConfiguration,
                     std::vector< SimplePolynomial >& fieldsAsPolynomials,
        std::vector< SimplePolynomial >& fieldDerivativesAsPolynomials ) const;

    // This sets pathParameterization to be repeated nodes of stepSizeFraction
    // times straightPath, less the reference field. There is no return value
    // optimization because we cannot be sure that poor physicist users will
    // have access to a C++11-compliant compiler.
    void
    InitialStepsForMinuit( std::vector< double >& pathParameterization,
                           std::vector< double > const& straightPath,
                           double const stepSizeFraction ) const;

    void SetReferenceField( size_t const referenceFieldIndex )
    { this->referenceFieldIndex = referenceFieldIndex; }

    size_t NumberOfVaryingPathNodes() const{ return numberOfVaryingPathNodes; }

    size_t NumberOfParameterizationFields() const
    { return numberOfParameterizationFields; }

    size_t ParameterizationSize() const
    { return ( numberOfParameterizationFields * numberOfVaryingPathNodes ); }



  protected:
    size_t const numberOfFields;
    size_t referenceFieldIndex;
    size_t const numberOfParameterizationFields;
    size_t const numberOfVaryingPathNodes;
    double const pathStepSize;
    Eigen::MatrixXd const pathStepInversion;

    // This sets up a ( numberOfVaryingPathNodes + 1 )-square matrix M where
    // M_ij = (i * pathStepSize)^j, then returns its inverse.
    Eigen::MatrixXd CreatePathStepPowersInverse() const;
  };




  // This sets pathParameterization to be repeated nodes of stepSizeFraction
  // times straightPath, less the reference field. There is no return value
  // optimization because we cannot be sure that poor physicist users will
  // have access to a C++11-compliant compiler.
  inline void PathFromNodes::InitialStepsForMinuit(
                                   std::vector< double >& pathParameterization,
                                     std::vector< double > const& straightPath,
                                          double const stepSizeFraction ) const
  {
    size_t actualFieldIndex( 0 );
    for( size_t nodeIndex( 0 );
         nodeIndex < numberOfVaryingPathNodes;
         nodeIndex += numberOfVaryingPathNodes )
    {
      for( size_t parameterizationFieldIndex( 0 );
           parameterizationFieldIndex < numberOfParameterizationFields;
           ++parameterizationFieldIndex )
      {
        actualFieldIndex = parameterizationFieldIndex;
        if( parameterizationFieldIndex >= referenceFieldIndex )
        {
          ++actualFieldIndex;
        }
        pathParameterization[ nodeIndex + parameterizationFieldIndex ]
        = ( stepSizeFraction * straightPath[ actualFieldIndex ] );
      }
    }
  }

  // This sets up a ( numberOfVaryingPathNodes + 1 )-square matrix M where
  // M_ij = (i * pathStepSize)^j, then returns its inverse.
  inline Eigen::MatrixXd
  PathFromNodes::CreatePathStepPowersInverse() const
  {
    Eigen::MatrixXd pathStepPowers( ( numberOfVaryingPathNodes + 1 ),
                                    ( numberOfVaryingPathNodes + 1 ) );
    pathStepPowers( numberOfVaryingPathNodes,
                    numberOfVaryingPathNodes ) = 1.0;
    for( unsigned int rowIndex( 0 );
         rowIndex < numberOfVaryingPathNodes;
         ++rowIndex )
    {
      pathStepPowers( numberOfVaryingPathNodes,
                      rowIndex ) = 1.0;
      for( unsigned int columnIndex( 0 );
           columnIndex <= numberOfVaryingPathNodes;
           ++columnIndex )
      {
        pathStepPowers( rowIndex,
                        columnIndex )
        = pow( ( (double)( rowIndex + 1 ) * pathStepSize ),
               ( columnIndex + 1 ) );
      }
    }
    return pathStepPowers.inverse();
  }

} /* namespace VevaciousPlusPlus */
#endif /* PATHFROMNODES_HPP_ */
