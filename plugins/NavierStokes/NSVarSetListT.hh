// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Physics_NavierStokes_VarSetListT_hh
#define COOLFluiD_Physics_NavierStokes_VarSetListT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"
#include "Framework/VarSetTransformerT.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////
    
/// This class binds together a list of diffusive variable sets (VS)
/// @author Ray Vandenhoeck
template <typename V1, typename V2>
class NSVarSetListT {
public:
  /// enumerators for dimension, equations number, physical data array size 
  enum {DIM=V1::DIM, NBEQS=V1::NBEQS, DATASIZE=V1::DATASIZE};
  
  typedef typename V1::PTERM PTERM;  /// typedefs for shared convective data  
  typedef typename V1::DTERM DTERM;  /// typedefs for shared diffusive data  
  typedef V1 SOLUTION_VS;            /// typedef for solution variable sets 
  typedef V2 UPDATE_VS;              /// typedef for update variable sets 
  
  // Constructor
  HOST_DEVICE NSVarSetListT(typename DTERM::template DeviceConfigOptions<NOTYPE>* dcoNS,typename PTERM::template DeviceConfigOptions<NOTYPE>* dco) 
  {
    m_solutionVS.setModelData(dcoNS, dco);
    m_updateVS.setModelData(dcoNS, dco);
    m_up2Sol.setModelData(dco);
  }
  
  // Destructor
  HOST_DEVICE ~NSVarSetListT() {}
  
  // Get the solution variable set
  HOST_DEVICE SOLUTION_VS* getSolutionVS() {return &m_solutionVS;} 
  
  // Get the update variable set
  HOST_DEVICE UPDATE_VS* getUpdateVS() {return &m_updateVS;} 
  
  // Get the update to solution variable transformer
  HOST_DEVICE Framework::VarSetTransformerT<UPDATE_VS, SOLUTION_VS, NOTYPE>* 
  getUpdateToSolution() {return &m_up2Sol;} 
  
private:
  
  /// solution variable set
  SOLUTION_VS m_solutionVS;
  
  /// update variable set
  UPDATE_VS   m_updateVS;
  
  /// variable set transformer from update to solution
  Framework::VarSetTransformerT<UPDATE_VS, SOLUTION_VS, NOTYPE> m_up2Sol;
};
    
//////////////////////////////////////////////////////////////////////////////

     } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_NavierStokes_NSVarSetListT_hh
