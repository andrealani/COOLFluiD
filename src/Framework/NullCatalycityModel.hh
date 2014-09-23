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

#ifndef COOLFluiD_Framework_NullCatalycityModel_hh
#define COOLFluiD_Framework_NullCatalycityModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/CatalycityModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a catalycity library
/// @author Andrea Lani
class Framework_API NullCatalycityModel : public Framework::CatalycityModel {
public:
  
  /// Constructor without arguments
  NullCatalycityModel(const std::string& name);
  
  /// Default destructor
  ~NullCatalycityModel();
  
  /// Compute the mass production
  /// @param temp  roto-translational temperature
  /// @param rho   density [kg/m^3]
  /// @param ys    array of species mass fractions
  /// @param mp    array of species mass production due to catalycity [kg m^-3 s^-1]
  void computeMassProduction(const CFreal temp, 
			     const CFreal rho, 
			     const RealVector& ys, 
			     RealVector& mp);
  
}; // end of class NullCatalycityModel
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullCatalycityModel_hh
