// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Numerics_NewtonMethod_CombustionTurbUpdateSol_hh
#define COOLFluiD_Numerics_NewtonMethod_CombustionTurbUpdateSol_hh

//////////////////////////////////////////////////////////////////////////////

#include "NewtonIteratorData.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
  /// This class represents a solution updater based on TurbUpdateSol, but
  /// valid for rho_i, u, v, T, k, omega variables instead of P, u, v, T, k, omega
  /// Author: Alessandro Mazzetti
	
class CombustionTurbUpdateSol : public NewtonIteratorCom {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit CombustionTurbUpdateSol(const std::string& name);

  /// Destructor.
  virtual ~CombustionTurbUpdateSol()
  {
  }

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Execute Processing actions
  virtual void execute();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected:

  // handle to states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  // handle to rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  // handle to update coefficient
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// the relaxation parameter vector
  std::vector<CFreal> _alpha;
  
  //free stream values of the turbulent variables
  CFreal _KInlet;
  CFreal _OmegaInlet;

  //clip value for temperature
  CFreal _TClip;
  
  /// pointer to the physical chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;

  ///// mass species fractions in the internal cell
  //RealVector _ysIn;
  
  /// Array containing an index for negative rho_i
  RealVector negativeRho_i;

}; // class CombustionTurbUpdateSol

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_NewtonMethod_CombustionTurbUpdateSol_hh
