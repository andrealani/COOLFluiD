// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ConvectionReactionPM_hh
#define COOLFluiD_Framework_ConvectionReactionPM_hh

//////////////////////////////////////////////////////////////////////////////



#include "PhysicalModelImpl.hh"
#include "BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a ConvectionReactionPM.
/// @author Andrea Lani
template <class CONVTERM, class REACTERM>
class ConvectionReactionPM : public PhysicalModelImpl {

public:

  /// Constructor without arguments
  ConvectionReactionPM(const std::string& name) :
    PhysicalModelImpl(name),
    _convectiveTerm(new CONVTERM("ConvTerm")),
    _sourceTerm(new REACTERM("SourceTerm"))
  {
  }

  /// Default destructor
  virtual ~ConvectionReactionPM()
  {
  }

  /// Get the convective name
  virtual std::string getConvectiveName() const = 0;

  /// Get the diffusive name
  virtual std::string getDiffusiveName() const = 0;

  /// Get the reaction name
  virtual std::string getSourceName() const = 0;

  /// @return the space dimension of the SubSystem
  virtual CFuint getDimension() const = 0;

  /// @return the number of equations of the SubSystem
  virtual CFuint getNbEquations() const = 0;

  /// Set the reference values
  virtual void setReferenceValues() = 0;

  /// Set the physical data
  virtual void computePhysicalData()
  {
    _convectiveTerm->setupPhysicalData();
    _sourceTerm->setupPhysicalData();
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
    return  _sourceTerm.get();
  }

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args )
  {
    PhysicalModelImpl::configure(args);
    configureNested ( _convectiveTerm.get(), args );
    configureNested ( _sourceTerm.get(), args );
  }

  /// Get the concrete convective term
  Common::SafePtr<CONVTERM> getConvTerm() const
  {
    return _convectiveTerm.get();
  }

  /// Get the concrete source (reaction) term
  Common::SafePtr<REACTERM> getSrcTerm() const
  {
    return _sourceTerm.get();
  }

private:

  /// Private Copy Constructor
  ConvectionReactionPM(const ConvectionReactionPM& v);

  /// Private Assignment operator
  const ConvectionReactionPM& operator=(const ConvectionReactionPM& v);

  /// convective term
  std::auto_ptr<CONVTERM> _convectiveTerm;

  /// reaction term
  std::auto_ptr<REACTERM> _sourceTerm;

}; // end of class ConvectionReactionPM

//////////////////////////////////////////////////////////////////////////////

  } // namespace ConvectionReactionPM

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ConvectionReactionPM_hh
