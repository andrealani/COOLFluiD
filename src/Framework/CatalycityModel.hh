// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_CatalycityModel_hh
#define COOLFluiD_Framework_CatalycityModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PhysicalPropertyLibrary.hh"
#include "MathTools/RealVector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class PhysicalChemicalLibrary;
    
//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a catalycity library
/// @author Andrea Lani
class Framework_API CatalycityModel : public Framework::PhysicalPropertyLibrary {
public:
  
  typedef Environment::ConcreteProvider<CatalycityModel,1> PROVIDER;
  typedef const std::string& ARG1;
  
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);
  
  /// Constructor without arguments
  CatalycityModel(const std::string& name);
  
  /// Default destructor
  virtual ~CatalycityModel();
  
  /// set up private data
  virtual void setup();
  
  /// Configures this configurable object.
  virtual void configure ( Config::ConfigArgs& args );
  
  /// Compute the mass production
  /// @param temp  roto-translational temperature
  /// @param rho   density [kg/m^3]
  /// @param ys    array of species mass fractions
  /// @param mp    array of species mass production due to catalycity [kg m^-3 s^-1]
  virtual void computeMassProduction(const CFreal temp, 
				     const CFreal rho, 
				     const RealVector& ys, 
				     RealVector& mp) = 0;
  
  /// Gets the Class name
  static std::string getClassName() { return "CatalycityModel"; }
  
protected:
  
  /// pointer to the physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
  
}; // end of class CatalycityModel
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_CatalycityModel_hh
