#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/State.hh"

#include "Common/BadValueException.hh"
#include "MathTools/MathConsts.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/BasePositivityIDP.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

void BasePositivityIDP::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("RelDensityFloor",
    "Density floor as a fraction of the cell mean density (default 1e-2). Dimensionless, so it needs no retuning between non-dimensionalizations.");
  options.addConfigOption< CFreal >("AbsDensityFloor",
    "Absolute density floor in code units (default 0, i.e. disabled).");
  options.addConfigOption< CFreal >("RelPressureFloor",
    "Pressure floor as a fraction of the cell mean pressure (default 1e-2).");
  options.addConfigOption< CFreal >("AbsPressureFloor",
    "Absolute pressure floor in code units (default 0, i.e. disabled).");
  options.addConfigOption< CFreal >("Delta",
    "Safety margin keeping the floor strictly below the cell mean, which is what makes theta lie in [0,1] by construction (default 1e-6).");
  options.addConfigOption< CFuint >("MaxVerifyIters",
    "Maximum bisection steps in the verify corrector (default 5).");
  options.addConfigOption< CFuint >("ShowRate",
    "Print frequency for limiter statistics (default 1).");
  options.addConfigOption< bool >("EnableDensityConstraint",
    "Enforce the density floor (default true).");
  options.addConfigOption< bool >("EnablePressureConstraint",
    "Enforce the pressure floor (default true). Set false for the density-only validation stage.");
  options.addConfigOption< std::string >("ScaleMode",
    "Which conservative components are scaled: Hydro (freeze B and phi, div B untouched), Full (scale everything), or Auto (try Hydro, fall back to Full). Default Auto.");
}

//////////////////////////////////////////////////////////////////////////////

BasePositivityIDP::BasePositivityIDP(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_posPrev("posPrev"),
  socket_volumes("volumes"),
  socket_outputPP("outputPP"),
  m_cellBuilder(CFNULL),
  m_cell(CFNULL),
  m_cellStates(CFNULL),
  m_updateToSolutionVecTrans(CFNULL),
  m_tmpState(CFNULL),
  m_tmpCoord(CFNULL),
  m_nbrEqs(0),
  m_dim(0),
  m_order(0),
  m_nbrSolPnts(0),
  m_nbrEvalPnts(0),
  m_nbrZSPnts(0),
  m_nbrLobatto(0),
  m_lobattoW1(0.),
  m_evalPntPolyVals(),
  m_cellAvgSolCoefs(CFNULL),
  m_consSolOrig(),
  m_consSol(),
  m_consEval(),
  m_consMean(),
  m_consBase(),
  m_tmpUpdate(),
  m_scaleMode(AUTO),
  m_nbLimited(0),
  m_nbMeanFailures(0),
  m_nbSequential(0),
  m_nbHydroFallbacks(0),
  m_nbVerifyBisections(0),
  m_nbVerifyToMean(0),
  m_minTheta(1.),
  m_minRho(0.),
  m_minP(0.),
  m_minPOut(0.)
{
  addConfigOptionsTo(this);

  m_relDensityFloor = 1.0e-2;
  setParameter("RelDensityFloor", &m_relDensityFloor);

  m_absDensityFloor = 0.0;
  setParameter("AbsDensityFloor", &m_absDensityFloor);

  m_relPressureFloor = 1.0e-2;
  setParameter("RelPressureFloor", &m_relPressureFloor);

  m_absPressureFloor = 0.0;
  setParameter("AbsPressureFloor", &m_absPressureFloor);

  m_delta = 1.0e-6;
  setParameter("Delta", &m_delta);

  m_maxVerifyIters = 5;
  setParameter("MaxVerifyIters", &m_maxVerifyIters);

  m_showrate = 1;
  setParameter("ShowRate", &m_showrate);

  m_enableDensity = true;
  setParameter("EnableDensityConstraint", &m_enableDensity);

  m_enablePressure = true;
  setParameter("EnablePressureConstraint", &m_enablePressure);

  m_scaleModeStr = "Auto";
  setParameter("ScaleMode", &m_scaleModeStr);
}

//////////////////////////////////////////////////////////////////////////////

BasePositivityIDP::~BasePositivityIDP()
{
}

//////////////////////////////////////////////////////////////////////////////

void BasePositivityIDP::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
BasePositivityIDP::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_posPrev);
  result.push_back(&socket_volumes);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
BasePositivityIDP::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_outputPP);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void BasePositivityIDP::setup()
{
  CFAUTOTRACE;

  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim    = PhysicalModelStack::getActive()->getDim();

  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);

  // One evaluation point set is built, for one element type.
  // So throw rather than silently check the wrong points.
  if (frLocalData.size() != 1)
  {
    throw BadValueException (FromHere(),
      "BasePositivityIDP: mixed element grids are not supported (one evaluation point set is built).");
  }

  m_order      = static_cast<CFuint>(frLocalData[0]->getPolyOrder());
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();

  m_cellAvgSolCoefs = frLocalData[0]->getCellAvgSolCoefs();

  m_updateToSolutionVecTrans = getMethodData().getUpdateToSolutionVecTrans();

  // The metric weighted cell mean needs |J| per solution point, which StdSetup
  // stores in the volumes socket, but only when this flag is on. Without it the
  // socket is empty and the mean would read past the end.
  if (!getMethodData().createVolumesSocket())
  {
    throw BadValueException (FromHere(),
      "BasePositivityIDP needs the volumes socket for the metric weighted cell mean. Set ComputeVolumeForEachState = true.");
  }

  if      (m_scaleModeStr == "Hydro") { m_scaleMode = HYDRO; }
  else if (m_scaleModeStr == "Full")  { m_scaleMode = FULL;  }
  else if (m_scaleModeStr == "Auto")  { m_scaleMode = AUTO;  }
  else
  {
    throw BadValueException (FromHere(),
      "BasePositivityIDP: ScaleMode must be Hydro, Full or Auto, got " + m_scaleModeStr);
  }

  buildEvalPointSet();

  m_consSolOrig.resize(m_nbrSolPnts);
  m_consSol.resize(m_nbrSolPnts);
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_consSolOrig[iSol].resize(m_nbrEqs);
    m_consSol[iSol].resize(m_nbrEqs);
  }
  m_consEval.resize(m_nbrEvalPnts);
  for (CFuint iPnt = 0; iPnt < m_nbrEvalPnts; ++iPnt)
  {
    m_consEval[iPnt].resize(m_nbrEqs);
  }
  m_consMean.resize(m_nbrEqs);
  m_consBase.resize(m_nbrEqs);
  m_tmpUpdate.resize(m_nbrEqs);

  // scratch State for the update to conservative transform. Coordinates are
  // set to a dummy node because some variable sets query them.
  m_tmpState = new State();
  m_tmpState->resize(m_nbrEqs);
  m_tmpCoord = new Node(false);
  m_tmpCoord->resize(m_dim);
  *m_tmpCoord = 0.0;
  m_tmpState->setSpaceCoordinates(m_tmpCoord);

  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  CFuint nbrCells = 0;
  for (CFuint iType = 0; iType < elemType->size(); ++iType)
  {
    nbrCells += (*elemType)[iType].getNbElems();
  }
  DataHandle< CFreal > output = socket_outputPP.getDataHandle();
  output.resize(nbrCells*m_nbrSolPnts);
  output = 0.0;

  setupPhysics();

  CFLog(NOTICE, "PositivityIDP: p = " << m_order
        << ", sol pnts = " << m_nbrSolPnts
        << ", eval pnts = " << m_nbrEvalPnts
        << " (" << m_nbrZSPnts << " Zhang-Shu interior)"
        << ", scale mode = " << m_scaleModeStr << "\n");

  if (m_nbrZSPnts > 0)
  {
    CFLog(NOTICE, "PositivityIDP: Zhang-Shu rule N = " << m_nbrLobatto
          << ", w1 = " << m_lobattoW1
          << " -> an explicit SSP run needs sum_d (dt*lambda_d/h_d) <= " << m_lobattoW1 << "\n");
  }

  if (!m_enablePressure)
  {
    CFLog(NOTICE, "PositivityIDP: WARNING pressure constraint DISABLED (density-only validation stage)\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void BasePositivityIDP::buildEvalPointSet()
{
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

  SafePtr< vector< RealVector > > solPntCoords = frLocalData[0]->getSolPntsLocalCoords();
  SafePtr< vector< RealVector > > flxPntCoords = frLocalData[0]->getFlxPntsLocalCoords();

  const CFuint nbrFlxPnts = flxPntCoords->size();

  // Evaluation point list:
  //   solution points   what the scheme samples for the discontinuous flux
  //   flux points       what the Riemann solver samples
  //   Zhang-Shu nodes   what the cell mean update argument needs
  // The first two must be admissible or the wave speeds produce NaN, which is
  // why they are enforced even though the Zhang-Shu argument does not ask for
  // the solution points.
  vector< RealVector > evalCoords;

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    evalCoords.push_back((*solPntCoords)[iSol]);
  }
  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    evalCoords.push_back((*flxPntCoords)[iFlx]);
  }

  // Points of the Gauss-Lobatto rule the Zhang-Shu argument needs: exactness of
  // degree >= p requires 2N-3 >= p, i.e. N >= (p+3)/2. In integer arithmetic
  // that is N = (p+4)/2, floored at 2.
  m_nbrLobatto = (m_order + 4)/2;
  if (m_nbrLobatto < 2) m_nbrLobatto = 2;

  // first Gauss-Lobatto weight, normalized so the weights sum to one
  m_lobattoW1 = 1.0/(static_cast<CFreal>(m_nbrLobatto)*static_cast<CFreal>(m_nbrLobatto - 1));

  // Interior Gauss-Lobatto abscissae, i.e. all nodes except -1 and +1. The
  // endpoints are already covered by the flux points. Values match
  // FluxReconstructionMethod/Lobatto.cxx, which returns the (ORDER+1) point rule.
  vector< CFreal > lobattoInterior;
  switch (m_nbrLobatto)
  {
    case 2:
      break; // endpoints only
    case 3:
      lobattoInterior.push_back(0.0);
      break;
    case 4:
      lobattoInterior.push_back(-0.4472135954999579392818);
      lobattoInterior.push_back(+0.4472135954999579392818);
      break;
    case 5:
      lobattoInterior.push_back(-0.6546536707079771437983);
      lobattoInterior.push_back(0.0);
      lobattoInterior.push_back(+0.6546536707079771437983);
      break;
    case 6:
      lobattoInterior.push_back(-0.765055323929464692851);
      lobattoInterior.push_back(-0.2852315164806450963142);
      lobattoInterior.push_back(+0.2852315164806450963142);
      lobattoInterior.push_back(+0.765055323929464692851);
      break;
    default:
      CFLog(WARN, "PositivityIDP: no Gauss-Lobatto table for N = " << m_nbrLobatto
            << " (p = " << m_order << "). Zhang-Shu interior nodes omitted; "
            << "enforcement falls back to solution and flux points only.\n");
      break;
  }

  // Directions along which the element is a tensor product, and along which the
  // Zhang-Shu construction therefore applies.
  const CFGeoShape::Type shape = frLocalData[0]->getShape();
  vector< CFuint > zsDirs;
  switch (shape)
  {
    case CFGeoShape::QUAD:
      zsDirs.push_back(KSI); zsDirs.push_back(ETA);
      break;
    case CFGeoShape::HEXA:
      zsDirs.push_back(KSI); zsDirs.push_back(ETA); zsDirs.push_back(ZTA);
      break;
    case CFGeoShape::PRISM:
      // Only the extrusion is a tensor product. In plane nothing is missing from
      // the point set: the triangle solution points are the Williams-Shunn points
      // and the cell average coefficients are their quadrature weights, all
      // positive, so the mean is a convex combination of values already enforced
      // here. What is missing is a boundary inclusive rule, i.e. Gauss-Lobatto
      // lines under the collapsed transform, which is what yields the explicit
      // CFL constant on a triangle.
      zsDirs.push_back(ZTA);
      CFLog(NOTICE, "PositivityIDP: prism element. Zhang-Shu nodes in the extrusion "
            "direction; in plane the solution points are the WS quadrature points of "
            "the cell mean, but the CFL constant there is not rigorous.\n");
      break;
    default:
      CFLog(WARN, "PositivityIDP: element shape has no tensor product direction. "
            "Enforcement uses solution and flux points only and the Zhang-Shu CFL "
            "constant does not apply.\n");
      break;
  }

  // For each tensor direction the Zhang-Shu nodes are
  //   {interior Lobatto in that direction} x {solution point coords elsewhere}.
  // The base points are the solution points whose coordinate in that direction
  // is minimal, which picks exactly one representative per combination of the
  // remaining coordinates. This needs neither the 1D tables nor any assumption
  // about solution point ordering, and it works for the prism, whose triangular
  // node coordinates are private to its element data class.
  for (CFuint iDir = 0; iDir < zsDirs.size(); ++iDir)
  {
    const CFuint d = zsDirs[iDir];

    CFreal dMin = MathTools::MathConsts::CFrealMax();
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      dMin = min(dMin, (*solPntCoords)[iSol][d]);
    }

    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      // exact comparison is intended, the coordinates come from a single table
      if ((*solPntCoords)[iSol][d] != dMin) continue;

      for (CFuint iLob = 0; iLob < lobattoInterior.size(); ++iLob)
      {
        RealVector pnt = (*solPntCoords)[iSol];
        pnt[d] = lobattoInterior[iLob];
        evalCoords.push_back(pnt);
        ++m_nbrZSPnts;
      }
    }
  }

  m_nbrEvalPnts = evalCoords.size();

  // Interpolation from solution points to every evaluation point. For the
  // solution points this returns the Kronecker delta by the Lagrange property,
  // so a single uniform matrix covers the whole set.
  m_evalPntPolyVals = frLocalData[0]->getSolPolyValsAtNode(evalCoords);

  cf_assert(m_evalPntPolyVals.size() == m_nbrEvalPnts);
}

//////////////////////////////////////////////////////////////////////////////

void BasePositivityIDP::unsetup()
{
  CFAUTOTRACE;

  if (m_tmpState != CFNULL)
  {
    m_tmpState->resetSpaceCoordinates();
    delete m_tmpCoord;
    m_tmpCoord = CFNULL;
    delete m_tmpState;
    m_tmpState = CFNULL;
  }
}

//////////////////////////////////////////////////////////////////////////////

void BasePositivityIDP::computeCellConsStates()
{
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    // transform() returns a pointer to a buffer overwritten by the next call,
    // so copy straight away
    const State* const consState =
      m_updateToSolutionVecTrans->transform((*m_cellStates)[iSol]);

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_consSolOrig[iSol][iEq] = (*consState)[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BasePositivityIDP::computeMetricMean()
{
  DataHandle< CFreal > volumes = socket_volumes.getDataHandle();

  // The discretely conserved quantity is sum_j c_j |J_j| U_j, with c_j the
  // normalized reference space quadrature weights and |J_j| the metric
  // Jacobian. StdSetup stores |J_j| in the volumes socket, sign flip included.
  m_consMean = 0.0;
  CFreal wSum = 0.0;

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    const CFreal w = (*m_cellAvgSolCoefs)[iSol]
                     * volumes[(*m_cellStates)[iSol]->getLocalID()];
    wSum += w;
    m_consMean += w*m_consSolOrig[iSol];
  }

  cf_assert(wSum > 0.0);
  m_consMean /= wSum;
}

//////////////////////////////////////////////////////////////////////////////

void BasePositivityIDP::interpolateToEvalPnts(const std::vector< RealVector >& consSol,
                                              std::vector< RealVector >& consEval) const
{
  for (CFuint iPnt = 0; iPnt < m_nbrEvalPnts; ++iPnt)
  {
    consEval[iPnt] = 0.0;
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      consEval[iPnt] += m_evalPntPolyVals[iPnt][iSol]*consSol[iSol];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BasePositivityIDP::computeTraceStates(std::vector< RealVector >& consEval)
{
  // What the scheme actually samples is Cons(interp_Upd(x)), not interp_Cons(x),
  // whenever the update variables are not the conservative set. The corrector
  // must check the former, so interpolate the stored states first and transform
  // afterwards.
  for (CFuint iPnt = 0; iPnt < m_nbrEvalPnts; ++iPnt)
  {
    m_tmpUpdate = 0.0;
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_tmpUpdate += m_evalPntPolyVals[iPnt][iSol]*(*((*m_cellStates)[iSol]));
    }

    m_tmpState->copyData(m_tmpUpdate);
    const State* const consState = m_updateToSolutionVecTrans->transform(m_tmpState);

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      consEval[iPnt][iEq] = (*consState)[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BasePositivityIDP::buildHydroBaseState(const RealVector& consPnt)
{
  // Base state of the Hydro scaling at this point: the mean for the scaled
  // components, the local value for the frozen ones. This is what theta = 0
  // gives, so it is both the chord base for the pressure constraint and the
  // admissibility test that decides the Auto fallback.
  m_consBase = consPnt;

  const std::vector< CFuint >& idx = scaledIndices(HYDRO);
  for (CFuint i = 0; i < idx.size(); ++i)
  {
    m_consBase[idx[i]] = m_consMean[idx[i]];
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal BasePositivityIDP::thetaDensity(const CFreal rhoBase,
                                       const CFreal rho,
                                       const CFreal epsRho) const
{
  if (rho >= epsRho) return 1.0;

  // rho is affine along the scaling segment, so this is exact. epsRho is kept
  // strictly below the mean by construction, so the denominator is positive
  // whenever this branch is taken.
  const CFreal denom = rhoBase - rho;
  cf_assert(denom > 0.0);
  return (rhoBase - epsRho)/denom;
}

//////////////////////////////////////////////////////////////////////////////

CFreal BasePositivityIDP::thetaPressure(const CFreal pBase,
                                        const CFreal p,
                                        const CFreal epsP) const
{
  if (p >= epsP) return 1.0;

  // p is concave on {rho > 0}, so along the scaling segment it lies above its
  // chord. Requiring the chord to clear the floor is therefore sufficient.
  // This is a bound and not the exact root of p(theta) = epsP, so it can only
  // under-estimate theta, which is the safe direction.
  const CFreal denom = pBase - p;
  if (denom <= 0.0) return 1.0;
  return (pBase - epsP)/denom;
}

//////////////////////////////////////////////////////////////////////////////

void BasePositivityIDP::applyThetaFromOrig(const CFreal theta, const ScaleMode mode)
{
  const std::vector< CFuint >& idx = scaledIndices(mode);

  // Always rebuild from the original state so repeated corrector trials do not
  // compound. Scalings about the same mean compose multiplicatively, so a
  // sequentially computed theta1*theta2 applied once here is identical to
  // applying the two stages in turn.
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_consSol[iSol] = m_consSolOrig[iSol];
    for (CFuint i = 0; i < idx.size(); ++i)
    {
      const CFuint iEq = idx[i];
      m_consSol[iSol][iEq] = m_consMean[iEq]
                             + theta*(m_consSolOrig[iSol][iEq] - m_consMean[iEq]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BasePositivityIDP::scaleEvalStates(const CFreal theta, const ScaleMode mode)
{
  const std::vector< CFuint >& idx = scaledIndices(mode);

  // The interpolation reproduces constants, so interpolating the scaled
  // solution point values is the same as scaling the interpolated values about
  // the same mean. That saves a full re-interpolation in the sequential path.
  for (CFuint iPnt = 0; iPnt < m_nbrEvalPnts; ++iPnt)
  {
    for (CFuint i = 0; i < idx.size(); ++i)
    {
      const CFuint iEq = idx[i];
      m_consEval[iPnt][iEq] = m_consMean[iEq]
                              + theta*(m_consEval[iPnt][iEq] - m_consMean[iEq]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BasePositivityIDP::writeBackStates()
{
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    consToUpdate(m_consSol[iSol], m_tmpUpdate);
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      (*((*m_cellStates)[iSol]))[iEq] = m_tmpUpdate[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

bool BasePositivityIDP::allAdmissible(const std::vector< RealVector >& consEval,
                                      const CFreal epsRho,
                                      const CFreal epsP) const
{
  for (CFuint iPnt = 0; iPnt < m_nbrEvalPnts; ++iPnt)
  {
    CFreal rho, p, B2;
    constraintsAtPoint(consEval[iPnt], rho, p, B2);

    if (m_enableDensity  && rho < epsRho) return false;
    if (m_enablePressure && p   < epsP  ) return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////

void BasePositivityIDP::execute()
{
  CFTRACEBEGIN;

  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  DataHandle< CFreal > output  = socket_outputPP.getDataHandle();
  DataHandle< CFreal > posPrev = socket_posPrev.getDataHandle();

  const bool hasArtVisc = getMethodData().hasArtificialViscosity();
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  if (hasArtVisc)
  {
    posPrev = MathTools::MathConsts::CFrealMax();
  }

  m_nbLimited          = 0;
  m_nbMeanFailures     = 0;
  m_nbSequential       = 0;
  m_nbHydroFallbacks   = 0;
  m_nbVerifyBisections = 0;
  m_nbVerifyToMean     = 0;
  m_minTheta           = 1.0;
  m_minRho             = MathTools::MathConsts::CFrealMax();
  m_minP               = MathTools::MathConsts::CFrealMax();
  m_minPOut            = MathTools::MathConsts::CFrealMax();

  const CFuint nbrElemTypes = elemType->size();

  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();

    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();
      m_cellStates = m_cell->getStates();

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        output[(*m_cellStates)[iSol]->getLocalID()] = 0.0;
      }

      computeCellConsStates();
      computeMetricMean();

      CFreal rhoMean, pMean, B2Mean;
      constraintsAtPoint(m_consMean, rhoMean, pMean, B2Mean);

      // An inadmissible cell mean is a scheme failure, not a limiter case.
      // Report it loudly instead of shifting the mean, which would inject mass
      // or energy and hide the real problem.
      if ((m_enableDensity && rhoMean <= 0.0) || (m_enablePressure && pMean <= 0.0))
      {
        ++m_nbMeanFailures;

        // Detail lines are for diagnosis, not for every stage of every step:
        // with a multi stage time integrator this command runs several times
        // per iteration, so gate them on the show rate and cap them per call.
        if (iter % m_showrate == 0 && m_nbMeanFailures <= 3)
        {
          CFLog(WARN, "PositivityIDP: INADMISSIBLE CELL MEAN, cell " << m_cell->getID()
                << " at " << (*m_cellStates)[0]->getCoordinates()
                << " : rho = " << rhoMean << ", p = " << pMean << "\n");
        }
        m_cellBuilder->releaseGE();
        continue;
      }

      // Local floors. The outer min keeps the floor strictly below the mean,
      // which is what puts theta in [0,1] without any clamping.
      const CFreal epsRho = min((1.0 - m_delta)*rhoMean,
                                max(m_absDensityFloor, m_relDensityFloor*rhoMean));
      const CFreal epsP   = min((1.0 - m_delta)*pMean,
                                max(m_absPressureFloor, m_relPressureFloor*pMean));

      interpolateToEvalPnts(m_consSolOrig, m_consEval);

      // Scan the pre-limiting state: global minima for reporting, the LLAV
      // feedback, and whether any raw density is non positive.
      CFreal cellPosPrev = MathTools::MathConsts::CFrealMax();
      CFreal cellMinP    = MathTools::MathConsts::CFrealMax();
      bool needsLim  = false;
      bool anyRhoNeg = false;

      for (CFuint iPnt = 0; iPnt < m_nbrEvalPnts; ++iPnt)
      {
        CFreal rho, p, B2;
        constraintsAtPoint(m_consEval[iPnt], rho, p, B2);

        m_minRho = min(m_minRho, rho);
        if (rho > 0.0) { m_minP = min(m_minP, p); cellMinP = min(cellMinP, p); }

        if (rho <= 0.0) anyRhoNeg = true;
        if (m_enableDensity  && rho < epsRho) needsLim = true;
        if (m_enablePressure && rho > 0.0 && p < epsP) needsLim = true;

        if (epsRho > 0.0) cellPosPrev = min(cellPosPrev, rho/epsRho);
        if (epsP   > 0.0 && rho > 0.0) cellPosPrev = min(cellPosPrev, p/epsP);
      }

      // LLAV/OB adds viscosity where positivity is at risk, so it must see the
      // pre-limiting state. Otherwise limiting a cell hides the fact that the
      // cell was in trouble in the first place.
      if (hasArtVisc)
      {
        posPrev[m_cell->getID()] = cellPosPrev;
      }

      if (!needsLim)
      {
        m_minPOut = min(m_minPOut, cellMinP);
        m_cellBuilder->releaseGE();
        continue;
      }

      // Scaling mode for this cell.
      ScaleMode mode = m_scaleMode;
      if (mode == AUTO)
      {
        mode = HYDRO;
        if (m_enablePressure)
        {
          // Hydro freezes B, so its base state carries the local field with the
          // mean hydrodynamic part. That base is not guaranteed admissible: it
          // fails where the local field is much stronger than the mean field,
          // which is exactly the low beta cells of interest.
          for (CFuint iPnt = 0; iPnt < m_nbrEvalPnts && mode == HYDRO; ++iPnt)
          {
            buildHydroBaseState(m_consEval[iPnt]);
            CFreal rhoB, pB, B2B;
            constraintsAtPoint(m_consBase, rhoB, pB, B2B);
            if (pB <= epsP) mode = FULL;
          }
        }
        if (mode == FULL) ++m_nbHydroFallbacks;
      }

      // Stage 1, density. Affine in theta, so valid whatever the sign of the
      // raw densities.
      CFreal theta = 1.0;
      if (m_enableDensity)
      {
        for (CFuint iPnt = 0; iPnt < m_nbrEvalPnts; ++iPnt)
        {
          CFreal rho, p, B2;
          constraintsAtPoint(m_consEval[iPnt], rho, p, B2);
          theta = min(theta, thetaDensity(rhoMean, rho, epsRho));
        }
      }

      if (m_enablePressure)
      {
        if (anyRhoNeg)
        {
          // With a non positive raw density the algebraic pressure is
          // meaningless (the -|m|^2/2rho term flips sign, so p can come out
          // spuriously large and positive and silently disarm the constraint)
          // and the concavity chord is invalid past the density factor. So
          // apply stage 1 to the evaluation states first, then take the
          // pressure chord on the density-limited state.
          ++m_nbSequential;

          const CFreal theta1 = theta;
          scaleEvalStates(theta1, mode);

          CFreal theta2 = 1.0;
          for (CFuint iPnt = 0; iPnt < m_nbrEvalPnts; ++iPnt)
          {
            CFreal rho, p, B2;
            constraintsAtPoint(m_consEval[iPnt], rho, p, B2);

            CFreal pBase = pMean;
            if (mode == HYDRO)
            {
              buildHydroBaseState(m_consEval[iPnt]);
              CFreal rhoB, pB, B2B;
              constraintsAtPoint(m_consBase, rhoB, pB, B2B);
              pBase = pB;
            }
            theta2 = min(theta2, thetaPressure(pBase, p, epsP));
          }

          theta = theta1*theta2;
        }
        else
        {
          // Single pass. Every constraint is affine or concave in theta and
          // holds at theta = 0, so taking the minimum is rigorous.
          for (CFuint iPnt = 0; iPnt < m_nbrEvalPnts; ++iPnt)
          {
            CFreal rho, p, B2;
            constraintsAtPoint(m_consEval[iPnt], rho, p, B2);

            CFreal pBase = pMean;
            if (mode == HYDRO)
            {
              buildHydroBaseState(m_consEval[iPnt]);
              CFreal rhoB, pB, B2B;
              constraintsAtPoint(m_consBase, rhoB, pB, B2B);
              pBase = pB;
            }
            theta = min(theta, thetaPressure(pBase, p, epsP));
          }
        }
      }

      if (theta < 0.0) theta = 0.0;

      applyThetaFromOrig(theta, mode);
      writeBackStates();

      // Corrector. The closed form theta is exact for interp_Cons, but the
      // scheme samples Cons(interp_Upd), and the two differ whenever the update
      // variables are not the conservative set. Verify against what the scheme
      // sees and bisect down if needed. theta = 0 is always admissible: it makes
      // every nodal value equal, and the Lagrange basis reproduces constants, so
      // the interpolant is then constant and equal to the cell mean. Note
      // admissibility is not monotone in theta, so this is a
      // safeguarded search rather than a bracketing bisection.
      computeTraceStates(m_consEval);

      CFuint nVerify = 0;
      while (!allAdmissible(m_consEval, epsRho, epsP))
      {
        if (nVerify >= m_maxVerifyIters)
        {
          theta = 0.0;
          ++m_nbVerifyToMean;
        }
        else
        {
          theta *= 0.5;
          ++m_nbVerifyBisections;
        }
        ++nVerify;

        applyThetaFromOrig(theta, mode);
        writeBackStates();
        computeTraceStates(m_consEval);

        if (theta == 0.0) break;
      }

      for (CFuint iPnt = 0; iPnt < m_nbrEvalPnts; ++iPnt)
      {
        CFreal rhoO, pO, B2O;
        constraintsAtPoint(m_consEval[iPnt], rhoO, pO, B2O);
        if (rhoO > 0.0) m_minPOut = min(m_minPOut, pO);
      }

      ++m_nbLimited;
      m_minTheta = min(m_minTheta, theta);

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        output[(*m_cellStates)[iSol]->getLocalID()] = 1.0 - theta;
      }

      m_cellBuilder->releaseGE();
    }
  }

  // reductions and reporting
  const std::string nsp = getMethodData().getNamespace();

  CFuint totLimited    = m_nbLimited;
  CFuint totMeanFail   = m_nbMeanFailures;
  CFuint totSequential = m_nbSequential;
  CFuint totHydroFall  = m_nbHydroFallbacks;
  CFuint totBisect     = m_nbVerifyBisections;
  CFuint totToMean     = m_nbVerifyToMean;
  CFreal totMinTheta   = m_minTheta;
  CFreal totMinRho     = m_minRho;
  CFreal totMinP       = m_minP;
  CFreal totMinPOut    = m_minPOut;

#ifdef CF_HAVE_MPI
  MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
  PE::GetPE().setBarrier(nsp);
  MPI_Allreduce(&m_nbLimited,          &totLimited,    1, MPI_UNSIGNED, MPI_SUM, comm);
  MPI_Allreduce(&m_nbMeanFailures,     &totMeanFail,   1, MPI_UNSIGNED, MPI_SUM, comm);
  MPI_Allreduce(&m_nbSequential,       &totSequential, 1, MPI_UNSIGNED, MPI_SUM, comm);
  MPI_Allreduce(&m_nbHydroFallbacks,   &totHydroFall,  1, MPI_UNSIGNED, MPI_SUM, comm);
  MPI_Allreduce(&m_nbVerifyBisections, &totBisect,     1, MPI_UNSIGNED, MPI_SUM, comm);
  MPI_Allreduce(&m_nbVerifyToMean,     &totToMean,     1, MPI_UNSIGNED, MPI_SUM, comm);
  MPI_Allreduce(&m_minTheta,           &totMinTheta,   1, MPI_DOUBLE,   MPI_MIN, comm);
  MPI_Allreduce(&m_minRho,             &totMinRho,     1, MPI_DOUBLE,   MPI_MIN, comm);
  MPI_Allreduce(&m_minP,               &totMinP,       1, MPI_DOUBLE,   MPI_MIN, comm);
  MPI_Allreduce(&m_minPOut,            &totMinPOut,    1, MPI_DOUBLE,   MPI_MIN, comm);
#endif

  if (PE::GetPE().GetRank(nsp) == 0)
  {
    // an inadmissible mean is a scheme failure, so it is always reported
    if (totMeanFail > 0)
    {
      CFLog(NOTICE, "PositivityIDP: " << totMeanFail << " INADMISSIBLE CELL MEANS\n");
    }

    if (iter % m_showrate == 0)
    {
      // Minima are of the state arriving here, which is what says how close the
      // scheme is to trouble. The limited state is bounded below by the floor by
      // construction, so its pressure is worth a line only when it is not.
      CFLog(NOTICE, "PositivityIDP: " << totLimited << " cells"
            << ", theta " << totMinTheta
            << ", min rho " << totMinRho
            << " p " << totMinP);
      if (totSequential > 0) CFLog(NOTICE, ", " << totSequential << " sequential");
      if (totHydroFall  > 0) CFLog(NOTICE, ", " << totHydroFall  << " Hydro->Full");
      if (totBisect     > 0) CFLog(NOTICE, ", " << totBisect     << " corrector bisections");
      if (totToMean     > 0) CFLog(NOTICE, ", " << totToMean     << " fell back to the mean");
      if (totMinPOut <= 0.0) CFLog(NOTICE, ", min p out " << totMinPOut << " WHICH IS A BUG");
      CFLog(NOTICE, "\n");
    }
  }

  PE::GetPE().setBarrier(nsp);

  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
