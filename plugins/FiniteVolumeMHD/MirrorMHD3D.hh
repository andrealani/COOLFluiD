#ifndef COOLFluiD_Numerics_FiniteVolume_MirrorMHD3D_hh
#define COOLFluiD_Numerics_FiniteVolume_MirrorMHD3D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD3DVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the mirror bc
   *
   * @author Mehmet Sarp Yalim
   * @author Andrea Lani
   *
   */
class MirrorMHD3D : public FVMCC_BC {

public:
  
  /**
   * Constructor
   */
  MirrorMHD3D(const std::string& name);

  /**
   * Default destructor
   */
  ~MirrorMHD3D();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);
  
 private:

  /// physical model var set
  Common::SafePtr<Physics::MHD::MHD3DVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;
  
  /// physical model data
  RealVector _dataGhostState;
  
}; // end of class MirrorMHD3D

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MirrorMHD3D_hh
