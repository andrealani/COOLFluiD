#include <fstream>
#include <algorithm>

#include "FluxReconstructionEntropyMHD/FluxReconstructionEntropyMHD.hh"
#include "FluxReconstructionEntropyMHD/BCSuperInletEntropyMHD.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/TopologicalRegionSet.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
#include "Common/SafePtr.hh"
#include "EntropyMHD/EntropyMHDTerm.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCSuperInletEntropyMHD, FluxReconstructionSolverData, BCStateComputer, FluxReconstructionEntropyMHDModule>
  BCSuperInletEntropyMHDProvider("SuperInletEntropyMHD");

//////////////////////////////////////////////////////////////////////////////

void BCSuperInletEntropyMHD::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("rhoFixed", "Fixed density at inlet (non-dimensional).");
  options.addConfigOption< CFreal >("pFixed", "Fixed pressure at inlet (non-dimensional).");
  options.addConfigOption< CFreal >("mX", "X-component of magnetic dipole moment (non-dimensional).");
  options.addConfigOption< CFreal >("mY", "Y-component of magnetic dipole moment (non-dimensional).");
  options.addConfigOption< CFreal >("mZ", "Z-component of magnetic dipole moment (non-dimensional).");
  options.addConfigOption< bool >("UseOpenClosedMask",
    "Enable open/closed field mask on Vr (default true). "
    "When false, Vr is extrapolated everywhere with V_tan=0.");
  options.addConfigOption< CFreal >("ROpenClosedClosed",
    "Radiality (|Br_face|/|B_face|) below which a point is considered fully "
    "in the closed-field region (V_face=0). Default 0.67 ~ 24 deg latitude for a dipole.");
  options.addConfigOption< CFreal >("ROpenClosedOpen",
    "Radiality above which a point is considered fully in the open-field region "
    "(V_face = Vr_int * r_hat). Default 0.83 ~ 36 deg latitude for a dipole.");
  options.addConfigOption< CFreal >("VrFixed",
    "Fixed radial velocity at inlet (non-dimensional). When > 0, uses this "
    "value instead of extrapolating Vr from interior. Removes BC feedback. Default 0 (= extrapolate).");
  options.addConfigOption< bool >("UseMagnetogram",
    "Read Br from a magnetogram .dat file instead of the analytic dipole formula. Default false.");
  options.addConfigOption< std::string >("MagnetogramFile",
    "Path to the magnetogram .dat file (COCONUT format: Cartesian point cloud on unit sphere).");
  options.addConfigOption< CFuint >("NbClosestPoints",
    "Number of closest points for inverse-distance-weighted Br interpolation. Default 8.");
}

//////////////////////////////////////////////////////////////////////////////

BCSuperInletEntropyMHD::BCSuperInletEntropyMHD(const std::string& name)
  : BCStateComputer(name)
{
  addConfigOptionsTo(this);

  m_rhoFixed = 1.0;
  setParameter("rhoFixed", &m_rhoFixed);

  m_pFixed = 0.03851;
  setParameter("pFixed", &m_pFixed);

  m_mX = 0.0;
  setParameter("mX", &m_mX);
  m_mY = 0.0;
  setParameter("mY", &m_mY);
  m_mZ = 0.0;
  setParameter("mZ", &m_mZ);

  m_useOpenClosedMask = true;
  setParameter("UseOpenClosedMask", &m_useOpenClosedMask);

  m_rClosed = 0.67;
  setParameter("ROpenClosedClosed", &m_rClosed);

  m_rOpen = 0.83;
  setParameter("ROpenClosedOpen", &m_rOpen);

  m_VrFixed = 0.0;
  setParameter("VrFixed", &m_VrFixed);

  m_useMagnetogram = false;
  setParameter("UseMagnetogram", &m_useMagnetogram);

  m_magnetogramFile = "";
  setParameter("MagnetogramFile", &m_magnetogramFile);

  m_nbClosestPts = 8;
  setParameter("NbClosestPoints", &m_nbClosestPts);

  m_nbrFaceFlxPntsMax = 0;
}

//////////////////////////////////////////////////////////////////////////////

BCSuperInletEntropyMHD::~BCSuperInletEntropyMHD()
{
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperInletEntropyMHD::setup()
{
  BCStateComputer::setup();

  m_needsSpatCoord = true;

  // Precompute sigmaFixed = pFixed * rhoFixed^(1-gamma)
  const CFreal gamma = Framework::PhysicalModelStack::getActive()->getImplementor()
    ->getConvectiveTerm().d_castTo<Physics::EntropyMHD::EntropyMHDTerm>()->getGamma();
  m_sigmaFixed = m_pFixed * std::pow(m_rhoFixed, 1.0 - gamma);

  if (m_useMagnetogram)
  {
    readMagnetogram();
    precomputeInletBr();
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperInletEntropyMHD::unsetup()
{
  BCStateComputer::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperInletEntropyMHD::computeGhostStates(
    const vector<State*>& intStates,
    vector<State*>& ghostStates,
    const vector<RealVector>& normals,
    const vector<RealVector>& coords)
{
  // BC structure:
  //  - phi: zero-gradient absorbing (ghost = +int) per Derigs 2018.
  //    Pair with the Dedner parabolic damping in EntropyMHDSourceTerm
  //    (DednerDamping = true) and a refSpeed >= max fast magnetosonic
  //    speed, otherwise the cleaning channel is bandwidth-limited.
  //  - rho, V, B, sigma: mirror-form Dirichlet (ghost = 2 face - int) so
  //    the LF Riemann central flux sees (int + ghost)/2 = face exactly.
  //  - alpha (open/closed mask): per-flux-point smoothstep on the local
  //    radiality |Br|/|B|. Each flux point evaluates its own mask value;
  //    no face-averaging. A face-averaged alpha misbehaves on faces that
  //    straddle the helmet-streamer transition: the single face value
  //    blends open and closed flux points together, feeding closed-side
  //    flux points a fractional-open velocity ghost, which drives spurious
  //    recirculation at the inlet.

  const CFuint nbrFlxPnts = intStates.size();
  const CFreal BmagEps = 1.0e-12;

  // Resolve the base offset for the precomputed Br cache once per face.
  // m_face is set by BCStateComputer::setFace() just before the call.
  CFuint faceLocalID = 0;
  if (m_useMagnetogram)
  {
    faceLocalID = m_globalToLocalTRSFaceID.find(m_face->getID());
  }

  for (CFuint iState = 0; iState < nbrFlxPnts; ++iState)
  {
    State& intState   = *intStates[iState];
    State& ghostState = *ghostStates[iState];

    const CFreal x = coords[iState][XX];
    const CFreal y = coords[iState][YY];
    const CFreal z = coords[iState][ZZ];
    const CFreal r = std::sqrt(x*x + y*y + z*z);
    const CFreal rInv = 1.0 / r;
    const CFreal rx = x * rInv;
    const CFreal ry = y * rInv;
    const CFreal rz = z * rInv;

    // Ghost state layout: (rho, u, v, w, Bx, By, Bz, sigma, phi)

    // ---- B: Dirichlet on Br (prescribed), zero-gradient on B_tan ----
    CFreal BrPrescribed;
    if (m_useMagnetogram)
    {
      BrPrescribed = m_flxPntBrs[faceLocalID * m_nbrFaceFlxPntsMax + iState];
    }
    else
    {
      const CFreal mDotr = m_mX*x + m_mY*y + m_mZ*z;
      const CFreal r4Inv = 1.0 / (r * r * r * r);
      BrPrescribed = 2.0 * mDotr * r4Inv;
    }

    const CFreal BrInt = intState[4]*rx + intState[5]*ry + intState[6]*rz;
    const CFreal deltaBr = BrPrescribed - BrInt;
    const CFreal Bx_face = intState[4] + deltaBr * rx;
    const CFreal By_face = intState[5] + deltaBr * ry;
    const CFreal Bz_face = intState[6] + deltaBr * rz;
    ghostState[4] = 2.0 * Bx_face - intState[4];
    ghostState[5] = 2.0 * By_face - intState[5];
    ghostState[6] = 2.0 * Bz_face - intState[6];

    const CFreal Bmag_face =
        std::sqrt(Bx_face*Bx_face + By_face*By_face + Bz_face*Bz_face);

    // ---- Per-flux-point alpha from the local radiality |Br|/|B| ----
    CFreal alpha = 1.0;
    if (m_useOpenClosedMask)
    {
      const CFreal radiality = (Bmag_face > BmagEps)
          ? std::abs(BrPrescribed) / Bmag_face
          : 0.0;

      if (m_rOpen > m_rClosed)
      {
        const CFreal tRaw = (radiality - m_rClosed) / (m_rOpen - m_rClosed);
        const CFreal t = std::max(0.0, std::min(1.0, tRaw));
        alpha = t * t * (3.0 - 2.0 * t);
      }
      else
      {
        alpha = (radiality >= m_rClosed) ? 1.0 : 0.0;
      }
    }

    // ---- rho, sigma: mirror-form Dirichlet ----
    // ghost = 2 * face - int makes (int + ghost)/2 = face exactly so the LF
    // central flux sees the prescribed Dirichlet state at this flux point.
    ghostState[0] = 2.0 * m_rhoFixed - intState[0];
    ghostState[7] = 2.0 * m_sigmaFixed - intState[7];

    // ---- Velocity: radial with V_tan = 0, Vr masked by alpha ----
    const CFreal Vr_raw = (m_VrFixed > 0.0) ? m_VrFixed
        : (intState[1]*rx + intState[2]*ry + intState[3]*rz);
    const CFreal Vr_target = alpha * Vr_raw;
    const CFreal Vx_face = Vr_target * rx;
    const CFreal Vy_face = Vr_target * ry;
    const CFreal Vz_face = Vr_target * rz;

    ghostState[1] = 2.0 * Vx_face - intState[1];
    ghostState[2] = 2.0 * Vy_face - intState[2];
    ghostState[3] = 2.0 * Vz_face - intState[3];

    // phi: absorbing (zero-gradient, ghost = +int) per Derigs 2018.
    ghostState[8] = intState[8];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperInletEntropyMHD::computeGhostGradients(
    const vector<vector<RealVector*> >& intGrads,
    vector<vector<RealVector*> >& ghostGrads,
    const vector<RealVector>& normals,
    const vector<RealVector>& coords)
{
  // B rows (4-6): Dirichlet-consistent mirror-normal (Br is pinned, so
  // the face-averaged gradient cancels its normal component). All other
  // rows: zero-gradient. At eta=0 the diffusive flux is identically zero,
  // so this rule is inert for the non-resistive bit-identity test.
  const CFuint nbrFlxPnts = intGrads.size();

  for (CFuint iState = 0; iState < nbrFlxPnts; ++iState)
  {
    const RealVector& n = normals[iState];
    const CFuint nbrVars = intGrads[iState].size();

    for (CFuint iVar = 0; iVar < nbrVars; ++iVar)
    {
      const RealVector& gInt = *intGrads[iState][iVar];
      RealVector& gGhost     = *ghostGrads[iState][iVar];

      if (iVar >= 4 && iVar <= 6)
      {
        const CFreal gDotN = gInt[XX]*n[XX] + gInt[YY]*n[YY] + gInt[ZZ]*n[ZZ];
        gGhost[XX] = gInt[XX] - 2.0 * gDotN * n[XX];
        gGhost[YY] = gInt[YY] - 2.0 * gDotN * n[YY];
        gGhost[ZZ] = gInt[ZZ] - 2.0 * gDotN * n[ZZ];
      }
      else
      {
        gGhost = gInt;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperInletEntropyMHD::readMagnetogram()
{
  // COCONUT magnetogram format:
  //   Line 1: nSurfaces (always 1)
  //   Line 2: !SURFACENAME  nPoints
  //   Lines 3+: x y z Br   (Cartesian on unit sphere, Br non-dimensionalized)
  cf_assert(!m_magnetogramFile.empty());

  std::ifstream file(m_magnetogramFile.c_str());
  if (!file.is_open())
  {
    CFLog(ERROR, "BCSuperInletEntropyMHD: cannot open magnetogram file: " << m_magnetogramFile << "\n");
    cf_assert(false);
  }

  CFuint nSurfaces = 0;
  file >> nSurfaces;

  std::string surfaceName;
  CFuint nPoints = 0;
  file >> surfaceName >> nPoints;

  cf_assert(nPoints > 0);

  m_magX.resize(nPoints);
  m_magY.resize(nPoints);
  m_magZ.resize(nPoints);
  m_magBr.resize(nPoints);

  for (CFuint i = 0; i < nPoints; ++i)
  {
    file >> m_magX[i] >> m_magY[i] >> m_magZ[i] >> m_magBr[i];
  }

  file.close();

  CFLog(INFO, "BCSuperInletEntropyMHD: loaded magnetogram with "
        << nPoints << " points from " << m_magnetogramFile << "\n");
}

//////////////////////////////////////////////////////////////////////////////

CFreal BCSuperInletEntropyMHD::interpolateBrMagnetogram(
    const CFreal x, const CFreal y, const CFreal z) const
{
  // Project to unit sphere (magnetogram cloud lives at r=1)
  const CFreal r = std::sqrt(x*x + y*y + z*z);
  const CFreal rInv = 1.0 / r;
  const CFreal xs = x * rInv;
  const CFreal ys = y * rInv;
  const CFreal zs = z * rInv;

  const CFuint nPts = m_magX.size();
  const CFuint nClosest = std::min(m_nbClosestPts, (CFuint)nPts);

  // Compute squared distances to all cloud points
  // Then partial-sort to find the nClosest nearest neighbors
  std::vector<std::pair<CFreal, CFuint> > dist2idx(nPts);
  for (CFuint i = 0; i < nPts; ++i)
  {
    const CFreal dx = xs - m_magX[i];
    const CFreal dy = ys - m_magY[i];
    const CFreal dz = zs - m_magZ[i];
    dist2idx[i] = std::make_pair(dx*dx + dy*dy + dz*dz, i);
  }

  std::partial_sort(dist2idx.begin(), dist2idx.begin() + nClosest, dist2idx.end());

  // Inverse-distance weighting
  CFreal wSum  = 0.0;
  CFreal BrSum = 0.0;
  for (CFuint k = 0; k < nClosest; ++k)
  {
    const CFreal dist = std::sqrt(dist2idx[k].first);
    const CFreal w = 1.0 / std::max(dist, (CFreal)1.0e-15);
    wSum  += w;
    BrSum += w * m_magBr[dist2idx[k].second];
  }

  return BrSum / wSum;
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperInletEntropyMHD::precomputeInletBr()
{
  using namespace Framework;

  // Locate the inlet TRS by name (m_trsNames[0] is set by BCStateComputer from
  // Simulator.SubSystem.Flow.Inlet.applyTRS in the CFcase).
  cf_assert(!m_trsNames.empty());
  const std::string inletName = m_trsNames[0];

  std::vector<Common::SafePtr<TopologicalRegionSet> > trsList =
      MeshDataStack::getActive()->getTrsList();
  Common::SafePtr<TopologicalRegionSet> thisTRS = CFNULL;
  for (CFuint iTRS = 0; iTRS < trsList.size(); ++iTRS)
  {
    if (trsList[iTRS]->getName() == inletName)
    {
      thisTRS = trsList[iTRS];
      break;
    }
  }
  cf_assert(thisTRS.isNotNull());

  // Face builder (face-to-cell geometric entity pool). The "second" face
  // builder is reserved for BC-side use; the primary one is owned by the
  // convective boundary command during residual assembly.
  Common::SafePtr<GeometricEntityPool<FaceToCellGEBuilder> > faceBuilder =
      getMethodData().getSecondFaceBuilder();

  Common::SafePtr<TopologicalRegionSet> cellTRS =
      MeshDataStack::getActive()->getTrs("InnerCells");

  FaceToCellGEBuilder::GeoData& faceData = faceBuilder->getDataGE();
  faceData.cellsTRS   = cellTRS;
  faceData.facesTRS   = thisTRS;
  faceData.isBoundary = true;

  // FR element data: we need face flx point local coords. For prisms both
  // triag (top/bottom) and quad (sides) face types exist with different
  // numbers of flx points, so we select per face via the "per type" table.
  std::vector<FluxReconstructionElementData*>& frLocalData =
      getMethodData().getFRLocalData();
  cf_assert(!frLocalData.empty());

  const CFGeoShape::Type elemShape = frLocalData[0]->getShape();
  const CFPolyOrder::Type order    = frLocalData[0]->getPolyOrder();

  Common::SafePtr<std::vector<RealVector> > defaultFlxLocalCoords =
      frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  Common::SafePtr<std::vector<std::vector<RealVector> > > flxCoordsPerType =
      frLocalData[0]->getFaceFlxPntsLocalCoordsPerType();

  // Max flx points across face types (quad side for prisms, single type otherwise).
  if (elemShape == CFGeoShape::PRISM)
  {
    m_nbrFaceFlxPntsMax = (order + 1) * (order + 1);
  }
  else
  {
    m_nbrFaceFlxPntsMax = defaultFlxLocalCoords->size();
  }

  const CFuint nbInletFaces = thisTRS->getLocalNbGeoEnts();
  m_flxPntBrs.assign(nbInletFaces * m_nbrFaceFlxPntsMax, 0.0);

  CFLog(INFO, "BCSuperInletEntropyMHD: precomputing Br at " << nbInletFaces
        << " inlet faces x up to " << m_nbrFaceFlxPntsMax
        << " flx pnts/face (magnetogram has " << m_magX.size() << " points)\n");

  for (CFuint iFace = 0; iFace < nbInletFaces; ++iFace)
  {
    faceData.idx = iFace;
    GeometricEntity* face = faceBuilder->buildGE();
    const CFuint faceGlobalID = face->getID();
    m_globalToLocalTRSFaceID.insert(faceGlobalID, iFace);

    // Pick correct flx pnt local coords depending on the face shape.
    const std::vector<RealVector>* flxLocalCoords = CFNULL;
    if (elemShape == CFGeoShape::PRISM)
    {
      flxLocalCoords = (face->getShape() == CFGeoShape::TRIAG)
          ? &((*flxCoordsPerType)[0])
          : &((*flxCoordsPerType)[1]);
    }
    else
    {
      flxLocalCoords = &(*defaultFlxLocalCoords);
    }
    const CFuint nbrFaceFlxPnts = flxLocalCoords->size();
    cf_assert(nbrFaceFlxPnts <= m_nbrFaceFlxPntsMax);

    for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
    {
      const RealVector worldCoord =
          face->computeCoordFromMappedCoord((*flxLocalCoords)[iFlx]);
      const CFreal Br = interpolateBrMagnetogram(
          worldCoord[XX], worldCoord[YY], worldCoord[ZZ]);
      m_flxPntBrs[iFace * m_nbrFaceFlxPntsMax + iFlx] = Br;
    }

    faceBuilder->releaseGE();
  }

  CFLog(INFO, "BCSuperInletEntropyMHD: Br precompute done ("
        << nbInletFaces * m_nbrFaceFlxPntsMax << " cached values)\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
