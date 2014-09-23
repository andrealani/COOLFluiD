#ifndef COOLFluiD_Numerics_FluctSplit_MHD2DPrimSourceTerm_hh
#define COOLFluiD_Numerics_FluctSplit_MHD2DPrimSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/ComputeSourceTermFSM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD { 
    
  namespace Physics {
    namespace MHD {
      class MHD2DVarSet;
    }
  }
  
 
    
    namespace FluctSplit {
      
      class InwardNormalsData;
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a MHD physical model 2D for conservative 
 * variables
 * 
 * @author Andrea Lani
 *
 *
 *
 */
class MHD2DSourceTerm : public ComputeSourceTermFSM {
  
public:
    
  /**
   * Constructor
   * @see MHD2D
   */
  MHD2DSourceTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  ~MHD2DSourceTerm();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  void setup();
  
  /**
   * Compute the source term
   */ 
  virtual void computeSourceFSM(Framework::GeometricEntity *const cell,
			                          RealVector& source,
			                          const FluctSplit::InwardNormalsData& normalsData);  
  
 private: // data
  
  /// Corresponding variable set
  Common::SafePtr<Physics::MHD::MHD2DVarSet> _varSet;
  
  /// physical data 
  RealVector m_physicalData;
  
}; // end of class MHD2DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit
  


} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_MHD2DSourceTerm_hh
