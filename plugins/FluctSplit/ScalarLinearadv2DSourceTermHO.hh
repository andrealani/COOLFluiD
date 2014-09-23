#ifndef COOLFluiD_Numerics_FluctSplit_ScalarLinearadv2DSourceTermHO_hh
#define COOLFluiD_Numerics_FluctSplit_ScalarLinearadv2DSourceTermHO_hh

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
class ScalarLinearadv2DSourceTermHO : public ComputeSourceTermFSM {
public:
  
  /**
   * Constructor
   * @see LinEuler2D
   */
  ScalarLinearadv2DSourceTermHO(const std::string& name);
  
  /**
   * Default destructor
   */
  ~ScalarLinearadv2DSourceTermHO();
  
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
  RealVector source_value;
  
   /// quadrature points per sub-element
  RealVector qd0;
  RealVector qd1;
  RealVector qd2;
  RealVector wqd;


  /// states at quadrature points
  std::vector<Framework::State*> qdstates;

  CFreal m_alpha;

 }; // end of class ScalarLinearadv2DSourceTermHO

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit
  


} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_ScalarRotationadv2DSourceTermMonopole_hh
