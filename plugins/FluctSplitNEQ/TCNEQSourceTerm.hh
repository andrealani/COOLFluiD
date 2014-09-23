#ifndef COOLFluiD_Numerics_FluctSplitNEQ_TCNEQSourceTerm_hh
#define COOLFluiD_Numerics_FluctSplitNEQ_TCNEQSourceTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/ComputeSourceTermFSM.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD { 
    
  namespace Framework {
    class State; 
    class GeometricEntity;
    class PhysicalChemicalLibrary;
  }
    
 

    namespace FluctSplit {
      class InwardNormalsData;
    }
    
    namespace FluctSplitNEQ {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a source term for Thermo Chemical NEQ physical model
 * 
 * @author Andrea Lani
 *
 */
template <class UPDATEVAR >
class TCNEQSourceTerm : public FluctSplit::ComputeSourceTermFSM {
  
public:
    
  /**
   * Constructor
   */
  TCNEQSourceTerm(const std::string& name);
  
  /**
   * Default destructor
   */
  ~TCNEQSourceTerm();
  
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
  
  /**
   * Compute the source term
   */ 
  void computeSourceFSM(Framework::GeometricEntity *const cell,
			std::vector<RealVector>& source,
			const FluctSplit::InwardNormalsData& normalsData);  
private:
  
  /**
   * Compute the TCNEQ source term starting from physical data
   */ 
  void computeTCNEQSource(CFuint iState,
			  const RealVector& pData,
			  RealVector& sourceBkp);
  
 private: // data
  
  /// corresponding variable set
  Common::SafePtr<UPDATEVAR> _varSet;
  
  /// physical-chemical library
  Common::SafePtr<Framework::PhysicalChemicalLibrary> _library;
  
  /// array storing all the physical data
  RealVector _pData;
  
  /// array to store the production/destruction term
  RealVector _omega;
  
  /// array to store the mass fractions
  RealVector _ys;
  
  /// vibrational temperatures
  RealVector _tvDim;
  
  /// array to store the production/destruction term
  RealVector _omegaTv;
  
  /// backup of the source per state
  std::vector<RealVector> _sourceVecBkp;
    
}; // end of class TCNEQSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplitNEQ
  


} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "TCNEQSourceTerm.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplitNEQ_TCNEQSourceTerm_hh
