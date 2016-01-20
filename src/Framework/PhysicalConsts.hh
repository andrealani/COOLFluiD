// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_PhysicalConsts_hh
#define COOLFluiD_Framework_PhysicalConsts_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Common/NonInstantiable.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// Provides an a set of static functions for physical constants
/// @author Andrea Lani
class PhysicalConsts : public Common::NonInstantiable<PhysicalConsts> {
public:
  
  /// Universal constant of gas [J/mol/K]
  static CFreal UnivRgas() {return 8.31447215;} 
  
  /// Boltzmann constant [J/K]
  static CFreal Boltzmann() {return 1.38065042e-23;} 
  
  /// Planck constant [J s]
  static CFreal Planck() {return 6.62606896e-34;} 
  
  /// Avogadro constant [1/mol]
  static CFreal Avogadro() {return 6.02214179e23;} 
  
  /// Speed of light [m/s]
  static CFreal LightSpeed() {return  299792458.;} 
  
  /// Stephan-Boltzmann constant [W m−2 K−4]
  static CFreal StephanBolzmann() {return  5.670373e-8;} 
  
  /// Magnetic permeability in vacuum [Volt·s/(Ampere·m)]
  static CFreal VacuumPermeability() {return 4.*MathTools::MathConsts::CFrealPi()*1e-7;} 
  
  /// Magnetic permittivity in vacuum [F/m]
  static CFreal VacuumPermittivity() {return 1/(LightSpeed()*LightSpeed()*VacuumPermeability());}
  
  /// Electron's mass [Kg]
  static CFreal ElectronMass() {return 9.1093821545e-31;} 

  /// Proton's mass [Kg]
  static CFreal ProtonMass() {return 1.67262177774e-27;}
  
  /// Hydrogen's mass [Kg]
  static CFreal HydrogenMass() {return 1.67262177774e-27+9.1093821545e-31;}

  /// Hydrogen's mass [C]
  static CFreal ElectronCharge() {return 1.602176620898e-19;}
};
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_PhysicalConsts_hh
