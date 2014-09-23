#ifndef COOLFluiD_Numerics_FluctSplit_LinearizedEuler3DSourceTermValiant_hh
#define COOLFluiD_Numerics_FluctSplit_LinearizedEuler3DSourceTermValiant_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/ComputeSourceTermFSM.hh"
#include "LinEuler/LinEulerTerm.hh"
#include "LinEuler/LinEuler3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class State; 
  }
  namespace FluctSplit {

    class InwardNormalsData;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a source term for LinearizedEuler3D MONOPOLE for the VALIANT project
 * 
 * @author  Lilla Edit Koloszar
 */
class LinearizedEuler3DSourceTermValiant : public ComputeSourceTermFSM {
public:
  
  /**
   * Constructor
   * @see LinEuler3D
   */
  LinearizedEuler3DSourceTermValiant(const std::string& name);
  
  /**
   * Default destructor
   */
  ~LinearizedEuler3DSourceTermValiant();
  
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
  Common::SafePtr<Physics::LinearizedEuler::LinEuler3DVarSet> _varSet;

  /// acquaintance of the PhysicalModel
  Common::SafePtr<Physics::LinearizedEuler::LinEulerTerm> _model;
  
  /// vector to store temporary result
  RealVector _temp;

  /// width of the source
  CFreal m_alpha;

  /// amplitude of the source
  CFreal m_eps;

  /// frequency of the source
  CFreal m_freq;

  /// location vector of the source
  std::vector<CFreal> _sourceloc;

 }; // end of class LinEuler3DSourceTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_LinearizedEuler3DSourceTermValiant_hh
