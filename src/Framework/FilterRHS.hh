// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_FilterRHS_hh
#define COOLFluiD_Framework_FilterRHS_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
#include "Framework/NumericalStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a configurable filter for residual values
/// @author Andrea Lani
class Framework_API FilterRHS : public Framework::NumericalStrategy {
public: // typedefs
  
  typedef Environment::ConcreteProvider<FilterRHS,1> PROVIDER;
  typedef const std::string& ARG1;

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  FilterRHS(const std::string& name);

  /// Destructor.
  virtual ~FilterRHS();

  /// Configure the data from the supplied arguments.
  /// @param args missing documentation
  virtual void configure ( Config::ConfigArgs& args );

  /// Filter the given value of the residual
  /// @param iEq    equation ID
  /// @param value  value to filter
  virtual void filter (CFuint iVar, CFreal& dU) = 0;
  
  /// Gets the Class name
  static std::string getClassName() { return "FilterRHS"; }

  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}

protected: // data

  /// array of flags telling if the variables must be filtered
  std::vector<bool> m_maskIDs;

}; // end of class FilterRHS

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

Framework_Factory(FilterRHS) // define the factoty instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_FilterRHS_hh
