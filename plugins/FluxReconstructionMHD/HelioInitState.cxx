#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "FluxReconstructionMHD/HelioInitState.hh"

#include "MHD/MHD3DProjectionVarSet.hh"
#include "MHD/MHDTerm.hh"

#include "Environment/DirPaths.hh"

#include <fstream>
#include <cmath>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<HelioInitState, FluxReconstructionSolverData, FluxReconstructionMHDModule>
HelioInitStateProvider("HelioInitState");

//////////////////////////////////////////////////////////////////////////////

void HelioInitState::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("FileNameTw",
    "Path to boundary surface data file (.dat).");
  options.addConfigOption< CFuint >("NbClosestPoints",
    "Number of closest points for angular interpolation (default: 6).");
  options.addConfigOption< CFreal >("InnerRadius",
    "Inner boundary radius in mesh units (0 = auto-detect from file, default: 0).");
  options.addConfigOption< CFreal >("VelocityExponent",
    "Radial velocity scaling exponent: Vr ~ r^alpha (default: 0 = constant Vr).");
  options.addConfigOption< bool >("IncludeParkerSpiral",
    "Wind Bphi via Parker spiral formula (default: true).");
  options.addConfigOption< CFreal >("SunRotationRate",
    "Non-dimensional solar rotation rate Omega*l0/V0 (default: 0.0409 for helio normalization).");
  options.addConfigOption< std::string >("InputVar",
    "Input variable set name (default: empty = same as update var).");
}

//////////////////////////////////////////////////////////////////////////////

HelioInitState::HelioInitState(const std::string& name) :
  FluxReconstructionSolverCom(name),
  m_fileNameTw(""),
  m_nbClosestPoints(6),
  m_innerRadius(0.0),
  m_velExponent(0.0),
  m_parkerSpiral(true),
  m_omegaSun(0.0409),
  m_inputVarStr(""),
  m_surfaces(),
  m_rBoundary(0.0),
  m_nbrEqs(0),
  m_dim(0),
  m_varSet(CFNULL),
  m_inputToUpdateVar(),
  m_initPntsStates(),
  m_inputState(CFNULL),
  m_initPntCoords(),
  m_tmpNode()
{
  addConfigOptionsTo(this);

  setParameter("FileNameTw",          &m_fileNameTw);
  setParameter("NbClosestPoints",     &m_nbClosestPoints);
  setParameter("InnerRadius",         &m_innerRadius);
  setParameter("VelocityExponent",    &m_velExponent);
  setParameter("IncludeParkerSpiral", &m_parkerSpiral);
  setParameter("SunRotationRate",     &m_omegaSun);
  setParameter("InputVar",            &m_inputVarStr);
}

//////////////////////////////////////////////////////////////////////////////

HelioInitState::~HelioInitState()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    HelioInitState::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void HelioInitState::configure(Config::ConfigArgs& args)
{
  CFAUTOTRACE;

  FluxReconstructionSolverCom::configure(args);

  // get the physical model
  std::string namespc = getMethodData().getNamespace();
  SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(namespc);
  SafePtr<PhysicalModel> physModel =
    PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  // get the name of the update variable set
  std::string updateVarStr = getMethodData().getUpdateVarStr();

  // default input vars to update vars if not specified
  if (m_inputVarStr.empty()) {
    m_inputVarStr = updateVarStr;
  }

  // create the transformer from input to update variables
  std::string provider =
    VarSetTransformer::getProviderName(
      physModel->getNameImplementor(), m_inputVarStr, updateVarStr);

  m_inputToUpdateVar =
    Environment::Factory<VarSetTransformer>::getInstance()
    .getProvider(provider)->create(physModel->getImplementor());
  cf_assert(m_inputToUpdateVar.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void HelioInitState::setup()
{
  CFAUTOTRACE;
  FluxReconstructionSolverCom::setup();

  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim    = PhysicalModelStack::getActive()->getDim();

  m_initPntCoords.resize(m_dim);
  m_tmpNode.resize(m_dim);

  // set up variable transformer
  const CFuint maxNbStatesInCell =
    MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  m_inputToUpdateVar->setup(maxNbStatesInCell);

  m_initPntsStates.resize(maxNbStatesInCell);
  for (CFuint iState = 0; iState < maxNbStatesInCell; ++iState)
  {
    m_initPntsStates[iState] = new State();
  }
  m_inputState = new State();

  m_varSet = getMethodData().getUpdateVar();

  // validate file name
  if (m_fileNameTw.empty()) {
    throw Common::BadValueException(FromHere(),
      "HelioInitState: FileNameTw not set. Specify the boundary surface data file.");
  }

  // read surface data
  CFLog(INFO, "HelioInitState: Reading surface data from " << m_fileNameTw << "\n");
  readSurfaceData(m_surfaces, m_fileNameTw);

  // auto-detect inner boundary radius if not specified
  if (m_innerRadius <= 0.0) {
    m_rBoundary = MathTools::MathConsts::CFrealMax();
    for (CFuint is = 0; is < m_surfaces.size(); ++is) {
      const SurfaceData& sf = *m_surfaces[is];
      const CFuint nbPoints = sf.rho.size();
      for (CFuint ip = 0; ip < nbPoints; ++ip) {
        CFreal r2 = 0.0;
        for (CFuint d = 0; d < m_dim; ++d) {
          r2 += sf.xyz(ip, d) * sf.xyz(ip, d);
        }
        m_rBoundary = std::min(m_rBoundary, std::sqrt(r2));
      }
    }
    CFLog(INFO, "HelioInitState: Auto-detected inner radius = "
          << m_rBoundary << "\n");
  } else {
    m_rBoundary = m_innerRadius;
    CFLog(INFO, "HelioInitState: Using configured inner radius = "
          << m_rBoundary << "\n");
  }

  if (m_rBoundary < 1.0e-14) {
    throw Common::BadValueException(FromHere(),
      "HelioInitState: Inner boundary radius is zero or negative.");
  }

  // log configuration
  CFLog(INFO, "HelioInitState: NbClosestPoints = " << m_nbClosestPoints << "\n");
  CFLog(INFO, "HelioInitState: VelocityExponent = " << m_velExponent << "\n");
  CFLog(INFO, "HelioInitState: IncludeParkerSpiral = " << m_parkerSpiral << "\n");
  CFLog(INFO, "HelioInitState: SunRotationRate = " << m_omegaSun << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void HelioInitState::unsetup()
{
  CFAUTOTRACE;

  for (CFuint iState = 0; iState < m_initPntsStates.size(); ++iState)
  {
    deletePtr(m_initPntsStates[iState]);
  }
  m_initPntsStates.resize(0);

  deletePtr(m_inputState);

  for (CFuint is = 0; is < m_surfaces.size(); ++is)
  {
    deletePtr(m_surfaces[is]);
  }
  m_surfaces.resize(0);

  FluxReconstructionSolverCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void HelioInitState::readSurfaceData(std::vector<SurfaceData*>& surfaces,
                                     const std::string& fileName)
{
  boost::filesystem::path fname;
  if (Common::StringOps::startsWith(fileName, "."))
  {
    fname = boost::filesystem::path(fileName);
  }
  else if (boost::filesystem::path(fileName).is_absolute())
  {
    fname = boost::filesystem::path(fileName);
  }
  else
  {
    fname = Environment::DirPaths::getInstance().getBaseDir() /
            boost::filesystem::path(fileName);
  }

  std::ifstream fin(fname.string().c_str());
  if (!fin) {
    throw Common::FilesystemException(FromHere(),
      "HelioInitState: Could not open file: " + fname.string());
  }

  // Format:
  // NumberOfSurfaces
  // SURFACE_NAME1 NumberOfPoints
  // x y z Rho u v w Bx By Bz P
  // ...

  CFuint nbSurf = 0;
  fin >> nbSurf;
  surfaces.resize(nbSurf);

  for (CFuint is = 0; is < nbSurf; ++is) {
    std::string nameSurf = "";
    fin >> nameSurf;
    CFuint nbPoints = 0;
    fin >> nbPoints;

    CFLog(INFO, "HelioInitState: Surface " << is << "/" << nbSurf
          << ", name = " << nameSurf << ", nbPoints = " << nbPoints << "\n");

    SurfaceData* sf = new SurfaceData();
    sf->xyz.resize(nbPoints, m_dim);
    sf->rho.resize(nbPoints);
    sf->u.resize(nbPoints);
    sf->v.resize(nbPoints);
    sf->w.resize(nbPoints);
    sf->Bx.resize(nbPoints);
    sf->By.resize(nbPoints);
    sf->Bz.resize(nbPoints);
    sf->p.resize(nbPoints);

    for (CFuint ip = 0; ip < nbPoints; ++ip) {
      for (CFuint iDim = 0; iDim < m_dim; ++iDim) {
        fin >> sf->xyz(ip, iDim);
      }
      fin >> sf->rho[ip];
      fin >> sf->u[ip];
      fin >> sf->v[ip];
      fin >> sf->w[ip];
      fin >> sf->Bx[ip];
      fin >> sf->By[ip];
      fin >> sf->Bz[ip];
      fin >> sf->p[ip];
    }

    surfaces[is] = sf;
  }

  CFLog(INFO, "HelioInitState: Finished reading " << nbSurf << " surface(s)\n");
}

//////////////////////////////////////////////////////////////////////////////

void HelioInitState::interpolateSurface(const RealVector& queryCoords,
                                        CFreal& rho, CFreal& u, CFreal& v,
                                        CFreal& w, CFreal& Bx, CFreal& By,
                                        CFreal& Bz, CFreal& p)
{
  // compute unit direction vector of query point
  const CFreal qx = queryCoords[0];
  const CFreal qy = queryCoords[1];
  const CFreal qz = (m_dim == DIM_3D) ? queryCoords[2] : 0.0;
  const CFreal qr = std::sqrt(qx*qx + qy*qy + qz*qz);

  CFreal qhat_x = 0.0, qhat_y = 0.0, qhat_z = 0.0;
  if (qr > 1.0e-14) {
    qhat_x = qx / qr;
    qhat_y = qy / qr;
    qhat_z = qz / qr;
  }

  ClosestPointData closestPoint;
  closestPoint.surfaceIDs.resize(m_nbClosestPoints);
  closestPoint.pointsIDs.resize(m_nbClosestPoints);
  closestPoint.r.resize(m_nbClosestPoints);
  closestPoint.reset();

  const CFuint nbSurf = m_surfaces.size();
  bool exactMatch = false;

  for (CFuint is = 0; is < nbSurf && (!exactMatch); ++is)
  {
    const SurfaceData& sf = *m_surfaces[is];
    const CFuint nbPoints = sf.rho.size();

    for (CFuint ip = 0; ip < nbPoints; ++ip)
    {
      // compute unit direction of boundary point
      CFreal px = sf.xyz(ip, 0);
      CFreal py = sf.xyz(ip, 1);
      CFreal pz = (m_dim == DIM_3D) ? sf.xyz(ip, 2) : 0.0;
      CFreal pr = std::sqrt(px*px + py*py + pz*pz);

      if (pr < 1.0e-14) continue;

      CFreal phat_x = px / pr;
      CFreal phat_y = py / pr;
      CFreal phat_z = pz / pr;

      // angular distance on unit sphere
      CFreal cosAngle = qhat_x*phat_x + qhat_y*phat_y + qhat_z*phat_z;
      cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
      CFreal angDist = std::acos(cosAngle);

      if (angDist < 1.0e-10)
      {
        // exact angular match
        rho = sf.rho[ip];
        u   = sf.u[ip];
        v   = sf.v[ip];
        w   = sf.w[ip];
        Bx  = sf.Bx[ip];
        By  = sf.By[ip];
        Bz  = sf.Bz[ip];
        p   = sf.p[ip];
        exactMatch = true;
        break;
      }

      // insert into k-nearest sorted list
      CFint counter = -1;
      for (CFuint n = 0; n < m_nbClosestPoints; ++n)
      {
        if (angDist < closestPoint.r[n]) {
          counter++;
        }
      }

      if (counter >= 0)
      {
        for (CFuint i = 0; i < static_cast<CFuint>(counter); ++i) {
          closestPoint.regressionFromTo(i + 1, i);
        }
        closestPoint.surfaceIDs[counter] = is;
        closestPoint.pointsIDs[counter]  = ip;
        closestPoint.r[counter]          = angDist;
      }
    }
  }

  if (!exactMatch)
  {
    // inverse-distance weighted interpolation
    CFreal wRho = 0.0, wU = 0.0, wV = 0.0, wW = 0.0;
    CFreal wBx = 0.0, wBy = 0.0, wBz = 0.0, wP = 0.0;
    CFreal sumWeights = 0.0;

    for (CFuint n = 0; n < m_nbClosestPoints; ++n)
    {
      const CFint idxs_s = closestPoint.surfaceIDs[n];
      const CFint idxp_s = closestPoint.pointsIDs[n];
      if (idxs_s < 0 || idxp_s < 0) continue; // unused slot

      const CFuint idxs = static_cast<CFuint>(idxs_s);
      const SurfaceData& sf = *m_surfaces[idxs];
      cf_assert(closestPoint.r[n] > 0.0);
      const CFreal weight = 1.0 / closestPoint.r[n];
      sumWeights += weight;

      const CFuint idxp = static_cast<CFuint>(idxp_s);
      wRho += weight * sf.rho[idxp];
      wU   += weight * sf.u[idxp];
      wV   += weight * sf.v[idxp];
      wW   += weight * sf.w[idxp];
      wBx  += weight * sf.Bx[idxp];
      wBy  += weight * sf.By[idxp];
      wBz  += weight * sf.Bz[idxp];
      wP   += weight * sf.p[idxp];
    }

    if (sumWeights > 0.0) {
      rho = wRho / sumWeights;
      u   = wU   / sumWeights;
      v   = wV   / sumWeights;
      w   = wW   / sumWeights;
      Bx  = wBx  / sumWeights;
      By  = wBy  / sumWeights;
      Bz  = wBz  / sumWeights;
      p   = wP   / sumWeights;
    } else {
      // fallback: should not happen
      rho = 1.0; u = 0.0; v = 0.0; w = 0.0;
      Bx = 0.0; By = 0.0; Bz = 0.0; p = 1.0;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void HelioInitState::applyParkerDecay(const RealVector& coords,
                                      CFreal r, CFreal theta,
                                      CFreal& rho, CFreal& u, CFreal& v,
                                      CFreal& w, CFreal& Bx, CFreal& By,
                                      CFreal& Bz, CFreal& p)
{
  const CFreal r_bc = m_rBoundary;

  // avoid division by zero at origin; no decay inside boundary
  if (r < 1.0e-14 || r_bc < 1.0e-14) return;
  if (r <= r_bc) return;

  const CFreal r_ratio = r / r_bc;       // > 1 for points outside boundary
  const CFreal inv_r_ratio = r_bc / r;   // < 1

  // --- Spherical unit vectors at this point ---
  const CFreal x = coords[0];
  const CFreal y = coords[1];
  const CFreal z = (m_dim == DIM_3D) ? coords[2] : 0.0;

  const CFreal sinTheta = std::sin(theta);
  const CFreal cosTheta = std::cos(theta);
  const CFreal rxy = std::sqrt(x*x + y*y);

  CFreal sinPhi = 0.0, cosPhi = 0.0;
  if (rxy > 1.0e-14) {
    cosPhi = x / rxy;
    sinPhi = y / rxy;
  } else {
    cosPhi = 1.0;
    sinPhi = 0.0;
  }

  // e_r     = (sinTheta*cosPhi, sinTheta*sinPhi, cosTheta)
  // e_theta = (cosTheta*cosPhi, cosTheta*sinPhi, -sinTheta)
  // e_phi   = (-sinPhi, cosPhi, 0)

  // --- Decompose velocity into spherical components ---
  const CFreal Vr     =  u*sinTheta*cosPhi + v*sinTheta*sinPhi + w*cosTheta;
  const CFreal Vtheta =  u*cosTheta*cosPhi + v*cosTheta*sinPhi - w*sinTheta;
  const CFreal Vphi   = -u*sinPhi          + v*cosPhi;

  // --- Decompose B into spherical components ---
  const CFreal Br     =  Bx*sinTheta*cosPhi + By*sinTheta*sinPhi + Bz*cosTheta;
  const CFreal Btheta =  Bx*cosTheta*cosPhi + By*cosTheta*sinPhi - Bz*sinTheta;
  const CFreal Bphi_bc = -Bx*sinPhi         + By*cosPhi;

  // --- Apply decay ---

  // Density: mass conservation rho * Vr * r^2 = const
  rho *= (inv_r_ratio * inv_r_ratio);

  // Velocity: Vr ~ r^alpha (default: constant), angular components decay as r_bc/r
  const CFreal Vr_new     = Vr * std::pow(r_ratio, m_velExponent);
  const CFreal Vtheta_new = Vtheta * inv_r_ratio;
  const CFreal Vphi_new   = Vphi * inv_r_ratio;

  // Magnetic field: Br ~ r^-2 (div B = 0), Btheta ~ r^-2
  const CFreal Br_new     = Br * (inv_r_ratio * inv_r_ratio);
  const CFreal Btheta_new = Btheta * (inv_r_ratio * inv_r_ratio);

  // Bphi: Parker spiral or simple decay
  CFreal Bphi_new;
  if (m_parkerSpiral && std::abs(sinTheta) > 1.0e-10)
  {
    // Parker spiral: Bphi = Br(r0) * (r0/r)^2 * (-Omega * r * sinTheta / Vr)
    const CFreal Vr_safe = (std::abs(Vr_new) > 1.0e-10) ? Vr_new : 1.0e-10;
    Bphi_new = Br * (inv_r_ratio * inv_r_ratio) *
               (-m_omegaSun * r * sinTheta / Vr_safe);
  }
  else
  {
    // simple decay
    Bphi_new = Bphi_bc * (inv_r_ratio * inv_r_ratio);
  }

  // Pressure: polytropic p ~ rho^gamma => p ~ r^(-2*gamma)
  SafePtr<MHDTerm> model = PhysicalModelStack::getActive()->getImplementor()
    ->getConvectiveTerm().d_castTo<MHDTerm>();
  const CFreal gamma = model->getGamma();
  p *= std::pow(inv_r_ratio, 2.0 * gamma);

  // --- Recompose to Cartesian ---
  u = Vr_new*sinTheta*cosPhi + Vtheta_new*cosTheta*cosPhi - Vphi_new*sinPhi;
  v = Vr_new*sinTheta*sinPhi + Vtheta_new*cosTheta*sinPhi + Vphi_new*cosPhi;
  w = Vr_new*cosTheta        - Vtheta_new*sinTheta;

  Bx = Br_new*sinTheta*cosPhi + Btheta_new*cosTheta*cosPhi - Bphi_new*sinPhi;
  By = Br_new*sinTheta*sinPhi + Btheta_new*cosTheta*sinPhi + Bphi_new*cosPhi;
  Bz = Br_new*cosTheta        - Btheta_new*sinTheta;
}

//////////////////////////////////////////////////////////////////////////////

void HelioInitState::executeOnTrs()
{
  CFAUTOTRACE;

  CFLog(INFO, "HelioInitState: Initializing heliospheric field...\n");

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

  // get the ElementTypeData
  SafePtr< vector<ElementTypeData> > elemType =
    MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrElemTypes = elemType->size();

  // get inner cells TRS
  SafePtr<TopologicalRegionSet> trs =
    MeshDataStack::getActive()->getTrs("InnerCells");

  // prepares to loop over cells
  SafePtr< GeometricEntityPool<StdTrsGeoBuilder> > geoBuilder =
    getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;

  CFuint stateCount = 0;

  // loop over element types
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbrElems = (*elemType)[iElemType].getNbElems();
    CFuint cellIdx = (*elemType)[iElemType].getStartIdx();

    // local coordinates of solution points
    SafePtr< vector< RealVector > > locInitPntCoords =
      frLocalData[iElemType]->getSolPntsLocalCoords();

    const CFuint nbrStates = frLocalData[iElemType]->getNbrOfSolPnts();

    // loop over elements
    for (CFuint iElem = 0; iElem < nbrElems; ++iElem, ++cellIdx)
    {
      // build the GeometricEntity
      geoData.idx = cellIdx;
      GeometricEntity *const cell = geoBuilder->buildGE();

      // get the states
      vector<State*>* solPntStates = cell->getStates();
      cf_assert(solPntStates->size() == nbrStates);

      // loop over solution points
      for (CFuint iPnt = 0; iPnt < nbrStates; ++iPnt)
      {
        // compute global coordinates
        m_initPntCoords = cell->computeCoordFromMappedCoord((*locInitPntCoords)[iPnt]);

        const CFreal x = m_initPntCoords[0];
        const CFreal y = m_initPntCoords[1];
        const CFreal z = (m_dim == DIM_3D) ? m_initPntCoords[2] : 0.0;
        const CFreal r = std::sqrt(x*x + y*y + z*z);
        const CFreal theta = std::atan2(std::sqrt(x*x + y*y), z);

        // interpolate boundary values at this angular position
        CFreal rho, u_val, v_val, w_val, Bx_val, By_val, Bz_val, p_val;
        interpolateSurface(m_initPntCoords, rho, u_val, v_val, w_val,
                           Bx_val, By_val, Bz_val, p_val);

        // apply radial decay from boundary to current radius
        applyParkerDecay(m_initPntCoords, r, theta,
                         rho, u_val, v_val, w_val,
                         Bx_val, By_val, Bz_val, p_val);

        // pack into input state (primitive MHD3DProjection)
        (*m_inputState)[0] = rho;
        (*m_inputState)[1] = u_val;
        (*m_inputState)[2] = v_val;
        (*m_inputState)[3] = w_val;
        (*m_inputState)[4] = Bx_val;
        (*m_inputState)[5] = By_val;
        (*m_inputState)[6] = Bz_val;
        (*m_inputState)[7] = p_val;
        (*m_inputState)[8] = 0.0;  // phi (GLM) = 0

        // transform to update variables and adimensionalize
        *m_initPntsStates[iPnt] = *m_inputToUpdateVar->transform(m_inputState);
        m_varSet->setAdimensionalValues(*m_initPntsStates[iPnt],
                                        *(*solPntStates)[iPnt]);

        stateCount++;
      }

      // release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }

  CFLog(INFO, "HelioInitState: Initialized " << stateCount
        << " solution points with radial extrapolation\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
