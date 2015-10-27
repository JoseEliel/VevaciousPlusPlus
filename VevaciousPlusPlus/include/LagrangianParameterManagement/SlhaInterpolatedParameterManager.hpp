/*
 * SlhaInterpolatedParameterManager.hpp
 *
 *  Created on: Oct 23, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef SLHAINTERPOLATEDPARAMETERMANAGER_HPP_
#define SLHAINTERPOLATEDPARAMETERMANAGER_HPP_

#include "CommonIncludes.hpp"
#include "LagrangianParameterManager.hpp"
#include "SlhaInterpolatedParameterFunctionoid.hpp"

namespace VevaciousPlusPlus
{

  class SlhaInterpolatedParameterManager : public LagrangianParameterManager
  {
  public:
    SlhaInterpolatedParameterManager();
    virtual ~SlhaInterpolatedParameterManager();


    // This should return the value of the requested parameter at the requested
    // scale (exp( logarithmOfScale )) as an alternative to adding it to the
    // list of parameters evaluated by ParameterValues. It is almost certain to
    // be much slower if used to obtain parameters repeatedly for different
    // scales at the same parameter point, but is more efficient if there are
    // some parameters which do not need to be evaluated for the potential but
    // still depend on Lagrangian parameters, for example when evaluating the
    // VEVs of the DSB vacuum for the parameter point.
    virtual double OnceOffParameter( std::string const& parameterName,
                                     double const logarithmOfScale ) = 0;

    // This checks to see if the parameter name has already been registered,
    // and if so, returns true paired with the index to its functionoid. If
    // not, it checks to see if it corresponds to a valid special case, and if
    // so, it adds the special case functionoid to the set of active
    // functionoids and returns true paired with the index of the newly active
    // functionoid. If not, then it checks to see if the parameter name
    // corresponds to an entry in a valid SLHA block, and if so, adds the
    // normal SLHA functionoid to the set of active functionoids and returns
    // true paired with the index of the newly active functionoid. If the
    // parameter name does not fall into any of the above cases, then false is
    // returned, paired with -1, rolling over to the maximum value of size_t as
    // it is unsigned.
    virtual std::pair< bool, size_t >
    RegisterParameter( std::string const& parameterName );

    // This returns a vector of the values of the Lagrangian parameters
    // evaluated at the given scale, ordered so that the indices given out by
    // RegisterParameter correctly match the parameter with its element in the
    // returned vector.
    virtual std::vector< double >
    ParameterValues( double const logarithmOfScale ) const
    { The order in which the parameters are filled is important as e.g.
    Bmu is mAsq/sinbeta cosbeta which will depend on tanbeta which is a derived
    functionoid so will have to have been registered earlier than Bmu and be
    evaluated earlier. }

    // This should write a function in the form
    // def LagrangianParameters( Q ): return ...
    // to return an array of the values of the Lagrangian parameters evaluated
    // at the scale Q, in the order in which a call to ParametersAtScale would
    // return them internal to this C++ code. By default it produces Python
    // code which throws an exception, as a derived class may not need to write
    // Python as the CosmoTransitions Python code may not be needed.
    virtual std::string ParametersAsPython() const
    { return std::string( "def LagrangianParameters( Q ): raise"
                          " NotImplementedError( \"A C++ class derived from"
                          " Vevacious::LagrangainParameterManager did not"
                          " override the ParametersAsPython() method to"
                          " provide valid Python code to return the values of"
                          " the Lagrangian parameters at the scale Q as an"
                          " array.\")"); }

    //
    SlhaInterpolatedParameterFunctionoid const& FunctionoidAtIndex


  protected:
    size_t numberOfDistinctActiveParameters;
    std::map< std::string, size_t > activeParametersToIndices;
    std::vector< SlhaInterpolatedParameterFunctionoid >
    activeNormalSlhaParameters;
    std::vector< SlhaDerivedFunctionoid* > activeSpecialCases;
    std::map< std::string, SlhaDerivedFunctionoid* > validSpecialCases;

    // This should prepare the LagrangianParameterManager for a new parameter
    // point. It is expected  that newInput is the name of a file containing
    // input required to evaluate the Lagrangian parameters (even the values
    // of the parameters themselves, such as given by an SLHA file), but
    // could itself be an encoding of the new input directly if a derived class
    // needs it. It could even just be ignored, if the derived class takes on
    // the role of deciding parameter input itself. Regardless, there should be
    // some code to prepare the parameters so that code objects which depend on
    // them can be informed that there is a new parameter point.
    virtual void PrepareNewParameterPoint( std::string const& newInput ) = 0;

    // This finds the block name part of parameterName, and either returns a
    // pointer to that block if valid, or a NULL pointer if not. If the block
    // is valid, indexString is set to be the index string from parameterName.
    LHPC::SLHA::SparseManyIndexedBlock< double > const*
    BlockFromParameterName( std::string const& parameterName,
                            std::string& indexString ) const;
  };

} /* namespace VevaciousPlusPlus */

#endif /* SLHAINTERPOLATEDPARAMETERMANAGER_HPP_ */
