#ifndef COOLFluiD_FluxReconstructionEntropyMHD_BCSuperInletEntropyMHD_hh
#define COOLFluiD_FluxReconstructionEntropyMHD_BCSuperInletEntropyMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include <vector>

#include "Common/CFMap.hh"
#include "Common/ConnectivityTable.hh"
#include "Common/SafePtr.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/RealVector.hh"
#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/FaceToCellGEBuilder.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Inlet BC for 3D entropy-variable MHD (EntropyMHD).
 *
 * Enforces (in primitive variables, mirror form ghost = 2 face - int):
 *   rho      = rho_fixed (Dirichlet)
 *   sigma    = sigmaFixed = pFixed*rhoFixed^(1-gamma) (Dirichlet)
 *   Br       = Br_dipole  (Dirichlet; prescribed background potential field)
 *   B_tan    = B_tan_int  (zero-gradient, COCONUT-style)
 *   V_tan    = 0          (no tangential flow at the base)
 *   V_r      = alpha * Vr_int  (alpha in [0,1] is an open/closed field mask)
 *   phi      : zero-gradient absorbing (ghost = +int) per Derigs 2018;
 *              pair with EntropyMHDSourceTerm.DednerDamping = true
 *
 * Magnetogram support:
 *   When UseMagnetogram=true, Br is read from a .dat file (COCONUT format:
 *   Cartesian point cloud on the unit sphere with non-dim Br) instead of
 *   the analytic dipole formula. Interpolation via inverse-distance
 *   weighting from NbClosestPoints nearest neighbors.
 *
 * Open/closed field mask (alpha):
 *   Uses runtime radiality = |Br_face| / |B_face| from the BC-consistent
 *   face B (Br from dipole Dirichlet, B_tan from interior). Adapts to
 *   any background field shape. Smoothstep between rClosed and rOpen.
 *   Setting UseOpenClosedMask=false reverts to alpha=1 everywhere.
 *
 * @author Rayan Dhib
 */
class BCSuperInletEntropyMHD : public BCStateComputer {
public:

  static void defineConfigOptions(Config::OptionList& options);

  BCSuperInletEntropyMHD(const std::string& name);
  virtual ~BCSuperInletEntropyMHD();

  virtual void setup();
  virtual void unsetup();

  void computeGhostStates(const std::vector<Framework::State*>& intStates,
                           std::vector<Framework::State*>& ghostStates,
                           const std::vector<RealVector>& normals,
                           const std::vector<RealVector>& coords);

  void computeGhostGradients(const std::vector<std::vector<RealVector*> >& intGrads,
                              std::vector<std::vector<RealVector*> >& ghostGrads,
                              const std::vector<RealVector>& normals,
                              const std::vector<RealVector>& coords);

private:

  /// fixed density (non-dimensional)
  CFreal m_rhoFixed;

  /// fixed pressure (non-dimensional)
  CFreal m_pFixed;

  /// magnetic dipole moment components (non-dimensional)
  CFreal m_mX;
  CFreal m_mY;
  CFreal m_mZ;

  /// precomputed sigmaFixed = pFixed * rhoFixed^(1-gamma)
  CFreal m_sigmaFixed;

  /// Open/closed field mask toggle (true = mask Vr, false = extrapolate Vr everywhere)
  bool m_useOpenClosedMask;

  /// Radiality threshold below which a point is fully in the closed-field region
  CFreal m_rClosed;

  /// Radiality threshold above which a point is fully in the open-field region
  CFreal m_rOpen;

  /// Fixed radial velocity at inlet (non-dimensional). When > 0, V_face uses
  /// this value instead of extrapolating from interior. Removes BC feedback.
  CFreal m_VrFixed;

  /// Magnetogram toggle and config
  bool m_useMagnetogram;
  std::string m_magnetogramFile;
  CFuint m_nbClosestPts;

  /// Magnetogram point cloud (Cartesian unit-sphere coords + non-dim Br)
  std::vector<CFreal> m_magX;
  std::vector<CFreal> m_magY;
  std::vector<CFreal> m_magZ;
  std::vector<CFreal> m_magBr;

  /// Precomputed Br at every inlet flux point, flat-indexed by
  /// localFaceID * m_nbrFaceFlxPntsMax + localFluxID.
  /// Filled once in setup() and read O(1) in computeGhostStates().
  std::vector<CFreal> m_flxPntBrs;

  /// Map from face global ID to local face ID (0..nInletFaces-1)
  Common::CFMap<CFuint, CFuint> m_globalToLocalTRSFaceID;

  /// Max number of flux points on any inlet face (for flat indexing)
  CFuint m_nbrFaceFlxPntsMax;

  /// Read magnetogram .dat file (COCONUT format)
  void readMagnetogram();

  /// Interpolate Br at (x,y,z) from the magnetogram point cloud via IDW
  CFreal interpolateBrMagnetogram(const CFreal x, const CFreal y, const CFreal z) const;

  /// Precompute m_flxPntBrs by walking the Inlet TRS once. Called from setup()
  /// when UseMagnetogram=true. Mirrors BCInletHyperPoisson's setup pattern -
  /// the magnetogram is static, so interpolation must not repeat every iter.
  void precomputeInletBr();

}; // end of class BCSuperInletEntropyMHD

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionEntropyMHD_BCSuperInletEntropyMHD_hh
