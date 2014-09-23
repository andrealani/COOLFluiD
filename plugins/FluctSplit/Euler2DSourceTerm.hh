#ifndef COOLFluiD_Numerics_FluctSplit_Euler2DSourceTerm_hh
#define COOLFluiD_Numerics_FluctSplit_Euler2DSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeSourceTermFSM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD { 
  
  namespace Framework {
    class State; 
  }
  
  namespace Physics {
    namespace NavierStokes {
      class Euler2DVarSet;
    }
  }
  
 
    
    namespace FluctSplit {
      
      class InwardNormalsData;
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a Euler physical model 2D for conservative 
 * variables
 * 
 * @author Andrea Lani
 */
class Euler2DSourceTerm : public ComputeSourceTermFSM {
public:
  
  /**
   * Constructor
   * @see Euler2D
   */
  Euler2DSourceTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  ~Euler2DSourceTerm();
  
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
  
  /// corresponding variable set
  Common::SafePtr<Physics::NavierStokes::Euler2DVarSet> _varSet;
  
  /// vector to store temporary result
  RealVector _temp;
  
 }; // end of class Euler2DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit
  


} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_Euler2DSourceTerm_hh
