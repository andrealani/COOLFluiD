// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Common/BadValueException.hh"
#include "Common/StringOps.hh"
#include "Environment/Factory.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/VarSetTransformer.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/ConvRHSFluxReconstructionSubcellBlending.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ConvRHSFluxReconstructionSubcellBlending, FluxReconstructionSolverData, FluxReconstructionModule >
  ConvRHSFluxReconstructionSubcellBlendingProvider("ConvRHSSubcellBlending");

//////////////////////////////////////////////////////////////////////////////

ConvRHSFluxReconstructionSubcellBlending::ConvRHSFluxReconstructionSubcellBlending(const std::string& name) :
  ConvRHSFluxReconstruction(name),
  socket_alpha("alpha"),
  socket_subcellFaceSamples("subcellFaceSamples"),
  m_scData(),
  m_currFaceAlphaF(0.0),
  m_currCellAlpha(0.0),
  m_flxPntFaceConn(CFNULL),
  m_faceRecStates()
{
  addConfigOptionsTo(this);

  m_faceFluxBlending = true;
  setParameter("FaceFluxBlending", &m_faceFluxBlending);

  m_reconstruction = "FirstOrder";
  setParameter("SubcellReconstruction", &m_reconstruction);

  m_limiter = "VanAlbada";
  setParameter("SubcellLimiter", &m_limiter);

  m_limiterEps = 1.0e-3;
  setParameter("SubcellLimiterEps", &m_limiterEps);

  m_reconstructionVar = "";
  setParameter("SubcellReconstructionVar", &m_reconstructionVar);
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("FaceFluxBlending",
    "Blend the face Riemann flux with the first-order flux between the adjacent solution points, "
    "weighted by max(alpha_L, alpha_R). With alpha = 1 the cell then runs a pure subcell P0 scheme.");
  options.addConfigOption< std::string >("SubcellReconstruction",
    "States at the subcell faces: FirstOrder (solution point values, default) or Linear "
    "(limited linear reconstruction to every subcell face, element faces included).");
  options.addConfigOption< std::string >("SubcellLimiter",
    "Slope limiter of the Linear reconstruction: VanAlbada (default, smooth, needed for Newton "
    "convergence), Minmod, or None (average of the two secants, unlimited, for verification only).");
  options.addConfigOption< std::string >("SubcellReconstructionVar",
    "Variables of the Linear reconstruction (a variable set name of the physical model, e.g. Puvt); "
    "default: the update variables. Primitive variables avoid negative pressures of conservative "
    "reconstructions at high Mach.");
  options.addConfigOption< CFreal >("SubcellLimiterEps",
    "VanAlbada smoothing size, relative to the largest value of the stencil (default 1e-3): "
    "differences below it are not limited, so smooth extrema are not clipped.");
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::configure(Config::ConfigArgs& args)
{
  ConvRHSFluxReconstruction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > >
ConvRHSFluxReconstructionSubcellBlending::needsSockets()
{
  std::vector< SafePtr< BaseDataSocketSink > > result = ConvRHSFluxReconstruction::needsSockets();
  result.push_back(&socket_alpha);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSource > >
ConvRHSFluxReconstructionSubcellBlending::providesSockets()
{
  std::vector< SafePtr< BaseDataSocketSource > > result = ConvRHSFluxReconstruction::providesSockets();
  result.push_back(&socket_subcellFaceSamples);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConvRHSFluxReconstructionSubcellBlending::getCellAlpha(const std::vector< State* >& states)
{
  DataHandle< CFreal > alpha = socket_alpha.getDataHandle();
  return alpha[states[0]->getLocalID()];
}

//////////////////////////////////////////////////////////////////////////////

CFreal ConvRHSFluxReconstructionSubcellBlending::computeFaceAlphaF()
{
  if (!m_faceFluxBlending) return 0.0;

  // both neighbours evaluate the same expression, so the face flux stays single-valued
  return std::max(getCellAlpha(*m_states[LEFT]), getCellAlpha(*m_states[RIGHT]));
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::storeFaceSamples()
{
  // the sample of a cell at a flux point is the trace of the other cell there
  DataHandle< CFreal > samples = socket_subcellFaceSamples.getDataHandle();
  const CFuint nbrFlxPnts = m_scData.getNbrFlxPnts();
  for (CFuint side = 0; side < 2; ++side)
  {
    const CFuint other  = 1 - side;
    const CFuint cellID = m_cells[side]->getID();
    for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
    {
      const CFuint flxIdx = (*m_faceFlxPntConnPerOrient)[m_orient][side][iFlxPnt];
      const CFuint start  = m_nbrEqs*(cellID*nbrFlxPnts + flxIdx);
      const State& trace  = *(m_cellStatesFlxPnt[other][iFlxPnt]);
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        samples[start+iEq] = trace[iEq];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

const RealVector& ConvRHSFluxReconstructionSubcellBlending::computeFaceLoFlux(const CFuint iFlxPnt)
{
  // solution points adjacent to this flux point in the left and right cell
  const CFuint flxIdxL = (*m_faceFlxPntConnPerOrient)[m_orient][LEFT ][iFlxPnt];
  const CFuint flxIdxR = (*m_faceFlxPntConnPerOrient)[m_orient][RIGHT][iFlxPnt];

  if (!m_scData.isLinear())
  {
    return m_riemannFluxComputer->computeFlux(*((*m_states[LEFT ])[m_scData.getClosestSol(flxIdxL)]),
                                              *((*m_states[RIGHT])[m_scData.getClosestSol(flxIdxR)]),
                                              m_unitNormalFlxPnts[iFlxPnt]);
  }

  // each side reconstructs its closest solution point to the face, with the trace of
  // the other side as outer sample
  const CFuint flxIdx[2] = {flxIdxL, flxIdxR};
  for (CFuint side = 0; side < 2; ++side)
  {
    m_scData.reconstructAtElementFace(*m_states[side], flxIdx[side], *(m_cellStatesFlxPnt[1-side][iFlxPnt]),
                                      *m_faceRecStates[side]);
  }
  return m_riemannFluxComputer->computeFlux(*m_faceRecStates[LEFT], *m_faceRecStates[RIGHT],
                                            m_unitNormalFlxPnts[iFlxPnt]);
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::computeInterfaceFlxCorrection()
{
  // face Riemann flux from the extrapolated states, unscaled in m_flxPntRiemannFlux
  // and scaled by the face Jacobian in m_cellFlx
  ConvRHSFluxReconstruction::computeInterfaceFlxCorrection();

  // needed by the cell loop also when this face uses no first-order flux
  if (m_scData.isLinear()) storeFaceSamples();

  m_currFaceAlphaF = computeFaceAlphaF();
  if (m_currFaceAlphaF <= 0.0) return;

  const CFreal oneMinusAlphaF = 1.0 - m_currFaceAlphaF;

  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    // first-order flux between the adjacent solution points
    const RealVector& loFlux = computeFaceLoFlux(iFlxPnt);

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_flxPntRiemannFlux[iFlxPnt][iEq] = oneMinusAlphaF*m_flxPntRiemannFlux[iFlxPnt][iEq]
                                        + m_currFaceAlphaF*loFlux[iEq];
    }

    // interface flux in the mapped coordinate frame of each neighbour
    m_cellFlx[LEFT ][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][LEFT ];
    m_cellFlx[RIGHT][iFlxPnt] = (m_flxPntRiemannFlux[iFlxPnt])*m_faceJacobVecSizeFlxPnts[iFlxPnt][RIGHT];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::addSubcellFaceFlux(const CFuint side, const CFreal alpha,
                                                                  vector< RealVector >& corrections)
{
  // each face flux point is the boundary of the subcell of its closest solution point.
  // m_cellFlx holds the face flux times the face Jacobian, oriented along the positive
  // mapped direction of the face normal.
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrFaceFlxPnts; ++iFlxPnt)
  {
    const CFuint flxIdx  = (*m_faceFlxPntConnPerOrient)[m_orient][side][iFlxPnt];
    const CFuint solIdx  = m_scData.getClosestSol(flxIdx);
    const CFint  faceDir = (*m_faceLocalDir)[(*m_flxPntFaceConn)[flxIdx]];

    const CFreal factor = alpha*static_cast<CFreal>(faceDir)/m_scData.getFaceSubcellWidth(flxIdx);

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      corrections[solIdx][iEq] -= factor*m_cellFlx[side][iFlxPnt][iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::computeCorrection(CFuint side, vector< RealVector >& corrections)
{
  // FR correction -(F* divh) of all solution points
  ConvRHSFluxReconstruction::computeCorrection(side, corrections);

  const CFreal alpha = getCellAlpha(*m_states[side]);
  if (alpha <= 0.0) return;

  const CFreal oneMinusAlpha = 1.0 - alpha;
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    corrections[iSolPnt] *= oneMinusAlpha;
  }

  addSubcellFaceFlux(side, alpha, corrections);
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::computeDivDiscontFlx(vector< RealVector >& residuals)
{
  // FR volume term -div(F_D) + div(h F_D) of all solution points
  ConvRHSFluxReconstruction::computeDivDiscontFlx(residuals);

  const CFreal alpha = m_currCellAlpha;
  if (alpha <= 0.0) return;

  const CFreal oneMinusAlpha = 1.0 - alpha;
  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    residuals[iSolPnt] *= oneMinusAlpha;
  }

  m_scData.computeSubcellRes(alpha, *m_cellStates, *m_riemannFluxComputer);
  const vector< RealVector >& subcellRes = m_scData.getSubcellRes();

  for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      residuals[iSolPnt][iEq] += subcellRes[iSolPnt][iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::setCellData()
{
  // flux projection vectors at the solution points
  ConvRHSFluxReconstruction::setCellData();

  m_currCellAlpha = getCellAlpha(*m_cellStates);
  if (m_currCellAlpha <= 0.0) return;

  m_scData.computeCellNormals(m_cell);

  if (m_scData.isLinear())
  {
    DataHandle< CFreal > samples = socket_subcellFaceSamples.getDataHandle();
    const CFuint start = m_nbrEqs*m_cell->getID()*m_scData.getNbrFlxPnts();
    cf_assert(start < samples.size());
    m_scData.computeCellSlopes(*m_cellStates, &samples[start]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::setup()
{
  CFAUTOTRACE;

  ConvRHSFluxReconstruction::setup();

  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() == 1);

  m_scData.setup(frLocalData[0], m_dim, m_nbrEqs);
  // telescoped subcell normals from the FR derivative, correction terms included
  m_scData.setCorrectionFunction(m_corrFctDiv);
  cf_assert(m_scData.getNbrSolPnts1D()*m_scData.getNbrSolPnts1D() == m_nbrSolPnts);

  m_flxPntFaceConn = frLocalData[0]->getFlxPntFaceConn();

  if (m_reconstruction == "Linear")
  {
    m_scData.setReconstruction(m_limiter, m_limiterEps, getMethodData().getUpdateVar());

    const std::string updateVar = getMethodData().getUpdateVarStr();
    if (m_reconstructionVar != "" && m_reconstructionVar != updateVar)
    {
      SafePtr< PhysicalModel > physModel = PhysicalModelStack::getActive();
      const std::string toRecStr =
        VarSetTransformer::getProviderName(physModel->getConvectiveName(), updateVar, m_reconstructionVar);
      const std::string fromRecStr =
        VarSetTransformer::getProviderName(physModel->getConvectiveName(), m_reconstructionVar, updateVar);
      m_toRecTrans.reset(Environment::Factory< VarSetTransformer >::getInstance().getProvider(toRecStr)
                         ->create(physModel->getImplementor()));
      m_fromRecTrans.reset(Environment::Factory< VarSetTransformer >::getInstance().getProvider(fromRecStr)
                           ->create(physModel->getImplementor()));
      m_toRecTrans->setup(1);
      m_fromRecTrans->setup(1);
      m_scData.setReconstructionVars(m_toRecTrans.getPtr(), m_fromRecTrans.getPtr());
    }

    SafePtr< TopologicalRegionSet > cells = MeshDataStack::getActive()->getTrs("InnerCells");
    socket_subcellFaceSamples.getDataHandle().resize
      (cells->getLocalNbGeoEnts()*m_scData.getNbrFlxPnts()*m_nbrEqs);

    RealVector dummyCoord(m_dim);
    dummyCoord = 0.0;
    for (CFuint side = 0; side < 2; ++side)
    {
      m_faceRecStates.push_back(new State());
      m_faceRecStates[side]->setSpaceCoordinates(new Node(dummyCoord, false));
    }
  }
  else if (m_reconstruction != "FirstOrder")
  {
    throw BadValueException(FromHere(),
      "ConvRHSFluxReconstructionSubcellBlending: SubcellReconstruction must be FirstOrder or Linear, not "
      + m_reconstruction);
  }

  CFLog(INFO, "ConvRHSFluxReconstructionSubcellBlending: SubcellReconstruction = " << m_reconstruction
        << (m_reconstruction == "Linear" ? " (" + m_limiter + ", eps " + StringOps::to_str(m_limiterEps)
                                           + ", variables " + (m_reconstructionVar == "" ?
                                             getMethodData().getUpdateVarStr() : m_reconstructionVar) + ")"
                                         : std::string(""))
        << ", FaceFluxBlending = " << m_faceFluxBlending << ", 1D subcell widths =");
  const vector< CFreal >& widths = m_scData.getWidths1D();
  for (CFuint i = 0; i < widths.size(); ++i)
  {
    CFLog(INFO, " " << widths[i]);
  }
  CFLog(INFO, "\n");
}

//////////////////////////////////////////////////////////////////////////////

void ConvRHSFluxReconstructionSubcellBlending::unsetup()
{
  CFAUTOTRACE;

  for (CFuint side = 0; side < m_faceRecStates.size(); ++side)
  {
    // the state owns its coordinates (Node built with isOnMesh = false)
    deletePtr(m_faceRecStates[side]);
  }
  m_faceRecStates.clear();

  ConvRHSFluxReconstruction::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD
