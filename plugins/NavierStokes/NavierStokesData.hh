#ifndef COOLFluiD_Physics_NavierStokes_NavierStokesData_hh
#define COOLFluiD_Physics_NavierStokes_NavierStokesData_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealVector.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Physics { 

    namespace NavierStokes { 
    
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the data for the NavierStokes physical model
 * 
 * @author Andrea Lani
 */
class NavierStokesData {

public:
  
  /**
   * Overloading of the assignment operator
   */
  const NavierStokesData& operator=(const NavierStokesData& data) 
  {
    dynViscosity = data.dynViscosity; // dynamic viscosity
    thermConductivity = data.thermConductivity; // thermal conductivity
    reynolds = data.reynolds; // Reynolds number
    
    return *this;
  }
  
  /**
   * Set up private data
   */
  void setup()
  {
  }
  
  CFreal dynViscosity; // dynamic viscosity
  
  CFreal thermConductivity; // thermal conductivity
  
  CFreal reynolds; // Reynolds number
  
};

//////////////////////////////////////////////////////////////////////////////
 
    } // namespace NavierStokes
 
  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NavierStokesData_hh
