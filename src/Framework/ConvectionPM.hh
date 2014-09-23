// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ConvectionPM_hh
#define COOLFluiD_Framework_ConvectionPM_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PhysicalModelImpl.hh"
#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a Convective physical model.
/// @author Andrea Lani
template <class CONVTERM>
class ConvectionPM : public PhysicalModelImpl {

public:

  /// Constructor without arguments
  ConvectionPM(const std::string& name) :
    PhysicalModelImpl(name),
    _convectiveTerm(new CONVTERM("ConvTerm"))
  {
  }

  /// Default destructor
  virtual ~ConvectionPM()
  {
  }

  /// Get the convective name
  virtual std::string getConvectiveName() const = 0;

  /// Get the diffusive name
  virtual std::string getDiffusiveName() const = 0;
  
  /// @return the space dimension of the SubSystem
  virtual CFuint getDimension() const = 0;

  /// @return the number of equations of the SubSystem
  virtual CFuint getNbEquations() const = 0;

  /// Set the reference values
  virtual void setReferenceValues() = 0;

  /// Set the reference value for time
  virtual void setReferenceTime() = 0;

  /// Set the physical data
  virtual void computePhysicalData()
  {
    _convectiveTerm->setupPhysicalData();
  }

  /// Get the convective term
  Common::SafePtr<BaseTerm> getConvectiveTerm() const
  {
    return _convectiveTerm.get();
  }

  /// Get the diffusive term
  Common::SafePtr<BaseTerm> getDiffusiveTerm() const
  {
    return CFNULL;
  }

  /// Get the source term
  Common::SafePtr<BaseTerm> getSourceTerm() const
  {
    return CFNULL;
  }

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    PhysicalModelImpl::configure(args);
    configureNested ( _convectiveTerm.get(), args );
  }

  /// Get the concrete convective term
  Common::SafePtr<CONVTERM> getConvTerm() const
  {
    return _convectiveTerm.get();
  }

private:

  /// Private Copy Constructor
  ConvectionPM(const ConvectionPM& v);

  /// Private Assignment operator
  const ConvectionPM& operator=(const ConvectionPM& v);

  /// convective term
  std::auto_ptr<CONVTERM> _convectiveTerm;

}; // end of class ConvectionPM

//////////////////////////////////////////////////////////////////////////////

  } // namespace ConvectionPM

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ConvectionPM_hh
