#ifndef COOLFluiD_Numerics_FluctSplitNEQ_TCNEQAxiSourceTerm_hh
#define COOLFluiD_Numerics_FluctSplitNEQ_TCNEQAxiSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/NavierStokes2DSourceTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD { 
    
  namespace Framework {
    class PhysicalChemicalLibrary;
  }
    
 
        
    namespace FluctSplitNEQ {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a source term for Thermo Chemical NEQ physical model
 * 
 * @author Andrea Lani
 *
 */
template <class UPDATEVAR>
class TCNEQAxiSourceTerm : public FluctSplit::NavierStokes2DSourceTerm<UPDATEVAR> {
  
public:
    
  /**
   * Constructor
   */
  TCNEQAxiSourceTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  ~TCNEQAxiSourceTerm();
  
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
  
private: // data
    
  /// physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> m_library;
    
  /// array to store the production/destruction term
  RealVector m_omega;
  
  /// array to store the mass fractions
  RealVector m_ys;
  
  /// vibrational temperatures
  RealVector m_tvDim;
  
  /// array to store the production/destruction term
  RealVector m_omegaTv;
      
}; // end of class TCNEQAxiSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplitNEQ
  


} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "TCNEQAxiSourceTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplitNEQ_TCNEQAxiSourceTerm_hh
