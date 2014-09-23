#ifndef COOLFluiD_Numerics_FiniteVolume_MirrorMHD3DProjectionTanakaLinearDensityExtrapolation_hh
#define COOLFluiD_Numerics_FiniteVolume_MirrorMHD3DProjectionTanakaLinearDensityExtrapolation_hh

//////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_BC.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that applies the mirror bc in 3D
   * for projection scheme according to Tanaka's approach with linear
   * density extrapolation to mimic the exponential density distribution when
   * taking the gravity of the planet into account
   *
   * @author Mehmet Sarp Yalim
   * @author Andrea Lani
   *
   */
class MirrorMHD3DProjectionTanakaLinearDensityExtrapolation : public FVMCC_BC {

public:

  /**
   * Constructor
   */
  MirrorMHD3DProjectionTanakaLinearDensityExtrapolation(const std::string& name);

  /**
   * Default destructor
   */
  ~MirrorMHD3DProjectionTanakaLinearDensityExtrapolation();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
 private:
  
  /// socket for stencil
  Framework::DataSocketSink<std::vector<Framework::State*> > socket_stencil;
  
  /// vector of IDs of the neighbour cells aligned with the boundary cells 
  //  in the normal direction to the boundary faces (i.e. planet surface)
  std::vector<CFuint> _alignedNgbourStateID;

  /// vector of distances of aligned neighbour cells from the origin
  std::vector<CFreal> _rNgbourState;

  /// vector of places of aligned neighbour cells in the stencil of boundary
  //  cells
  std::vector<CFuint> _placeInStencil;
  
  /// physical model var set
  Common::SafePtr<Physics::MHD::MHD3DProjectionVarSet> _varSet;

  /// physical model data
  RealVector _dataInnerState;

  /// physical model data
  RealVector _dataGhostState;

}; // end of class MirrorMHD3DProjectionTanakaLinearDensityExtrapolation

//////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_MirrorMHD3DProjectionTanakaLinearDensityExtrapolation_hh
