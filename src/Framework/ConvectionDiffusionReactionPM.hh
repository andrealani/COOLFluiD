// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_ConvectionDiffusionReactionPM_hh
#define COOLFluiD_Framework_ConvectionDiffusionReactionPM_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PhysicalModelImpl.hh"
#include "Framework/BaseTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a ConvectionDiffusionReactionPM.
/// @author Andrea Lani
template <class CONVTERM, class DIFFTERM, class REACTERM>
class ConvectionDiffusionReactionPM : public PhysicalModelImpl {

public:

  /// Constructor without arguments
  ConvectionDiffusionReactionPM(const std::string& name) :
    PhysicalModelImpl(name),
    _convectiveTerm(new CONVTERM("ConvTerm")),
    _diffusiveTerm(new DIFFTERM("DiffTerm")),
    _sourceTerm(new REACTERM("SourceTerm"))
  {
  }

  /// Default destructor
  virtual ~ConvectionDiffusionReactionPM()
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
    _diffusiveTerm->setupPhysicalData();
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
    return _diffusiveTerm.get();
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
    configureNested ( _diffusiveTerm.get(), args );
    configureNested ( _sourceTerm.get(), args );
  }

  /// Get the concrete convective term
  Common::SafePtr<CONVTERM> getConvTerm() const
  {
    return _convectiveTerm.get();
  }

  /// Get the concrete diffusive term
  Common::SafePtr<DIFFTERM> getDiffTerm() const
  {
    return _diffusiveTerm.get();
  }

  /// Get the concrete source (reaction) term
  Common::SafePtr<REACTERM> getSrcTerm() const
  {
    return _sourceTerm.get();
  }

private:

  /// Private Copy Constructor
  ConvectionDiffusionReactionPM(const ConvectionDiffusionReactionPM& v);

  /// Private Assignment operator
  const ConvectionDiffusionReactionPM& operator=(const ConvectionDiffusionReactionPM& v);

  /// convective term
  std::auto_ptr<CONVTERM> _convectiveTerm;

  /// diffusive term
  std::auto_ptr<DIFFTERM> _diffusiveTerm;

  /// reaction term
  std::auto_ptr<REACTERM> _sourceTerm;

}; // end of class ConvectionDiffusionReactionPM

//////////////////////////////////////////////////////////////////////////////

  } // namespace ConvectionDiffusionReactionPM

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_ConvectionDiffusionReactionPM_hh
