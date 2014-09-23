#ifndef COOLFluiD_Numerics_FluctSplit_ScalarLinearadv2DSourceTerm_hh
#define COOLFluiD_Numerics_FluctSplit_ScalarLinearadv2DSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "ComputeSourceTermFSM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD { 
  
  namespace Framework {
    class State; 
  }
    
 
    
    namespace FluctSplit {
      
      class InwardNormalsData;
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a source term for LinearizedEuler2D MONOPOLE
 * 
 * @author Nadege Villledieu
 */
class ScalarLinearadv2DSourceTerm : public ComputeSourceTermFSM {
public:
  
  /**
   * Constructor
   * @see LinEuler2D
   */
  ScalarLinearadv2DSourceTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  ~ScalarLinearadv2DSourceTerm();
  
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

 static void defineConfigOptions(Config::OptionList& options);

private: // data
  /// Coefficient of the source term, if it is positive the source term is a productive source
  /// otherwise it is a dissipative one
  CFreal m_alpha;

  
 }; // end of class ScalarLinearadv2DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit
  


} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ScalarLinearadv2DSourceTermMonopole_hh
