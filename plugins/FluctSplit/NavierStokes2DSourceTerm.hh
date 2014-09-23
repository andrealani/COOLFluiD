#ifndef COOLFluiD_Numerics_FluctSplit_NavierStokes2DSourceTerm_hh
#define COOLFluiD_Numerics_FluctSplit_NavierStokes2DSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeSourceTermFSM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD { 
  
  namespace Framework {
    class State; 
  }
  
  namespace Physics {
    namespace NavierStokes {
      class NavierStokesVarSet;
    }
  }
  
 
    
    namespace FluctSplit {
      
      class InwardNormalsData;
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents the axisymmetric source term for the Navier-Stokes
 * equations
 * 
 * @author Andrea Lani
 */

template <class UPDATEVAR>
class NavierStokes2DSourceTerm : public ComputeSourceTermFSM {
public:
  
  /**
   * Constructor
   * @see NavierStokes2D
   */
  NavierStokes2DSourceTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~NavierStokes2DSourceTerm();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Compute the source term
   */ 
  virtual void computeSourceFSM(Framework::GeometricEntity *const cell,
				RealVector& source,
				const FluctSplit::InwardNormalsData& normalsData);  
  
protected: // data
  
  // acquaintance of the diffusive variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> m_diffVar;
  
  // acquaintance of the update convective var set
  Common::SafePtr<UPDATEVAR> m_updateVar;
  
  // average radius in the cell
  CFreal m_avRadius;
  
  // cell volume
  CFreal m_cellVolume;
  
  // array of cell states
  std::vector<RealVector*> m_states;
  
  // 2D array (matrix) of values (rho, u, v, T)
  RealMatrix m_values;
  
}; // end of class NavierStokes2DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit
  


} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes2DSourceTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_NavierStokes2DSourceTerm_hh
