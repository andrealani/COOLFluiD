#ifndef COOLFluiD_Numerics_FluctSplit_ScalarUnsteady2DSourceTermHO_hh
#define COOLFluiD_Numerics_FluctSplit_ScalarUnsteady2DSourceTermHO_hh

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
class ScalarUnsteady2DSourceTermHO : public ComputeSourceTermFSM {
public:
  
  /**
   * Constructor
   * @see LinEuler2D
   */
  ScalarUnsteady2DSourceTermHO(const std::string& name);
  
  /**
   * Default destructor
   */
  ~ScalarUnsteady2DSourceTermHO();
  
  /**
   * Set up private data and data of the aggregated classes 
   * in this command before processing phase
   */
  void setup();
  
  /**
   * Compute the source term
   */ 
  void computeSourceFSM(Framework::GeometricEntity *const cell,
				RealVector& source,
				const FluctSplit::InwardNormalsData& normalsData);  

  static void defineConfigOptions(Config::OptionList& options);
private: // data
   /// Coefficient of the source term, if it is positive the source term is a productive source
  /// otherwise it is a dissipative one
  CFreal m_alpha;
  RealVector source_value;
  
   /// quadrature points per sub-element
  RealVector qd0;
  RealVector qd1;
  RealVector qd2;
  RealVector wqd;


  /// states at quadrature points
  std::vector<Framework::State*> qdstates;

 }; // end of class LinEuler2DSourceTermSTM

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit
  


} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LinearizedEuler2DSourceTermMonopole_hh
