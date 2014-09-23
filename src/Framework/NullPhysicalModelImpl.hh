// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullPhysicalModelImpl_hh
#define COOLFluiD_Framework_NullPhysicalModelImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "PhysicalModelImpl.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class State;

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a NullPhysicalModelImpl.
/// @author Andrea Lani
/// @author Tiago Quintino

class Framework_API NullPhysicalModelImpl : public PhysicalModelImpl {
public:

  /// Constructor without arguments
  NullPhysicalModelImpl(const std::string& name);

  /// Default destructor
  ~NullPhysicalModelImpl();

  /// Set the max number of states data
  void setMaxNbStatesData(const CFuint maxNbStatesData)
  {
  }

  /// Get the convective name
  std::string getConvectiveName() const;

  /// Get the diffusive name
  std::string getDiffusiveName() const;

  /// Get the source name
  std::string getSourceName() const;

  /// @return the space dimension of the SubSystem
  CFuint getDimension() const;

  /// @return the number of equations of the SubSystem
  CFuint getNbEquations() const;

  /// Check if this state is in a valid state
  bool validate(const State& state) const;

  /// Get the convective term
  Common::SafePtr<BaseTerm> getConvectiveTerm() const
  {
    return CFNULL;
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

  /// Configure
  virtual void configure ( Config::ConfigArgs& args );

private:

  /// Set the reference values
  void setReferenceValues();

  /// Set the physical data
  void computePhysicalData();

}; // end of class NullPhysicalModelImpl

//////////////////////////////////////////////////////////////////////////////

  } // namespace NullPhysicalModelImpl

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullPhysicalModelImpl_hh
