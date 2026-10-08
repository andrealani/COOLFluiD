// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include <limits>
#include "Common/BadValueException.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/BaseOrderBlending.hh"
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

MethodCommandProvider<BaseOrderBlending, FluxReconstructionSolverData, FluxReconstructionModule>
    OrderBlendingFRProvider("OrderBlending");

//////////////////////////////////////////////////////////////////////////////

BaseOrderBlending::BaseOrderBlending(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_alpha("alpha"),
  socket_prevAlpha("prevAlpha"),
  socket_smoothness("smoothness"),
  m_cellBuilder(CFNULL),
  m_cell(CFNULL),
  m_cellStates(CFNULL),
  m_sweepSnapshot(),
  m_obUpdateVarSet(CFNULL),
  m_obPData(),
  m_tempSolPntVec(),
  m_tempSolPntVec2(),
  m_vdmInv(),
  m_maxModalOrder(),
  m_NeighborIDs(),
  m_s(0.0),
  m_alphaRelaxation(1.0),
  m_alphaInitialized(false),
  m_stage(0),
  m_lastIter(std::numeric_limits<CFuint>::max()),
  m_resHistory(),
  m_cflHistory(),
  m_order(0),
  m_nbrSolPnts(0),
  m_nbrEqs(0),
  m_dim(0),
  m_iElemType(0),
  m_elemIdx(0),
  m_bndFaceOrients(),
  m_bndFacesBuilt(false),
  m_faceFlxPntConn(CFNULL),
  m_solPolyValsAtFlxPnts(CFNULL)
{
  addConfigOptionsTo(this);

  // Reference smoothness: s0 = -S0 * log10(N+1). Higher S0 = more lenient.
  m_s0 = 4.0;
  setParameter("S0", &m_s0);

  // Transition half-width for the sinusoidal ramp in log-space.
  m_kappa = 1.5;
  setParameter("Kappa", &m_kappa);

  // Physics-agnostic monitored expression. Physics-specific expressions
  // (e.g. B2 for MHD) are handled by subclasses overriding extractMonitoredField.
  m_modalMonitoredExpression = "rho*p";
  setParameter("ModalMonitoredExpression", &m_modalMonitoredExpression);

  m_alphaMin = 0.01;
  setParameter("AlphaMin", &m_alphaMin);

  m_alphaMax = 1.0;
  setParameter("AlphaMax", &m_alphaMax);

  setParameter("RelaxationFactor", &m_alphaRelaxation);

  // Decay factor for neighbor alpha during Jacobi smoothing.
  m_neighborWeight = 0.5;
  setParameter("NeighborWeight", &m_neighborWeight);

  // Smoothing iterations beyond the initial spread: total passes = m_nbSweeps + 1.
  m_nbSweeps = 0;
  setParameter("NbSweeps", &m_nbSweeps);

  // Iteration to freeze alpha (reuses prevAlpha). Default: never freeze.
  m_freezeAlphaIter = 1000000;
  setParameter("FreezeAlphaIter", &m_freezeAlphaIter);

  // former name of FreezeAlphaIter, still read so that existing cases keep freezing
  m_freezeFilterIter = 1000000;
  setParameter("freezeFilterIter", &m_freezeFilterIter);

  // State variables checked for undershoots below the neighbours. Default: none.
  m_forceAlphaMinVars = std::vector<CFuint>();
  setParameter("ForceAlphaMinVars", &m_forceAlphaMinVars);

  m_forceAlphaMinMargin = std::log(2.0);
  setParameter("ForceAlphaMinMargin", &m_forceAlphaMinMargin);

  // Boundary faces as neighbours in the undershoot test. Default: off.
  m_forceAlphaMinBoundaryGhost = false;
  setParameter("ForceAlphaMinBoundaryGhost", &m_forceAlphaMinBoundaryGhost);

  // Automatic freeze on a stalled residual. Default: off.
  m_freezeStallWindow = 0;
  setParameter("FreezeAlphaStallWindow", &m_freezeStallWindow);
}

//////////////////////////////////////////////////////////////////////////////

BaseOrderBlending::~BaseOrderBlending()
{
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::configure(Config::ConfigArgs& args)
{
  FluxReconstructionSolverCom::configure(args);
  if (m_freezeFilterIter != 1000000)
  {
    CFLog(WARN, "OrderBlending: freezeFilterIter is deprecated, use FreezeAlphaIter\n");
  }
  if (!(m_alphaRelaxation > 0.0 && m_alphaRelaxation <= 1.0))
  {
    throw BadValueException(FromHere(),
      "OrderBlending RelaxationFactor must be in (0, 1].");
  }
  if (m_freezeStallWindow > 0 && getFreezeAlphaIter() < 1000000)
  {
    throw BadValueException(FromHere(),
      "OrderBlending: FreezeAlphaStallWindow and FreezeAlphaIter cannot be used together");
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal, Config::DynamicOption<> >("S0",
    "Reference smoothness: s0 = -S0 * log10(N+1). Higher = more dissipation. Typical: 3.0-5.0.");
  options.addConfigOption< CFreal, Config::DynamicOption<> >("Kappa",
    "Transition half-width for the sinusoidal ramp around s0. Typical: 1.0-2.0.");
  options.addConfigOption< std::string >("ModalMonitoredExpression",
    "Monitored scalar for smoothness detection. Base class accepts: "
    "'rho', 'p', 'rho*p', 'p/rho', 'rho/p', 'velocity_magnitude'. "
    "Physics subclasses may add more (e.g. 'B2' in OrderBlendingMHD).");
  options.addConfigOption< CFreal >("AlphaMin",
    "Lower dead-band threshold: alpha < AlphaMin snaps to 0.");
  options.addConfigOption< CFreal >("AlphaMax",
    "Maximum blending coefficient cap.");
  options.addConfigOption< CFreal >("RelaxationFactor",
    "Fraction of the new sensor field applied after spatial spreading: "
    "alpha = previous + factor*(requested - previous). Default 1.0 "
    "applies the full new field without temporal relaxation. The first evaluation initializes alpha directly.");
  options.addConfigOption< CFreal >("NeighborWeight",
    "Decay factor applied to neighbor alpha during Jacobi smoothing. "
    "0 disables spreading.");
  options.addConfigOption< CFuint >("NbSweeps",
    "Number of smoothing iterations beyond the initial spread. "
    "Total spreading passes = NbSweeps + 1.");
  options.addConfigOption< CFuint, Config::DynamicOption<> >("FreezeAlphaIter",
    "Iteration number at which alpha is frozen (reuses prevAlpha). "
    "The first evaluation always initializes alpha. Very large = never freeze.");
  options.addConfigOption< CFuint, Config::DynamicOption<> >("freezeFilterIter",
    "Deprecated name of FreezeAlphaIter. If both are given, the smaller iteration is used.");
  options.addConfigOption< std::vector<CFuint> >("ForceAlphaMinVars",
    "State variables (e.g. ln rho_i) checked for undershoots: a cell where one of them at a solution point lies "
    "more than ForceAlphaMinMargin below the smallest cell mean of its neighbours gets alpha = 1 for the rest of "
    "the run. Empty (default) disables it.");
  options.addConfigOption< CFreal >("ForceAlphaMinMargin",
    "Margin of the ForceAlphaMinVars undershoot test, in the units of the variables (default ln 2).");
  options.addConfigOption< bool >("ForceAlphaMinBoundaryGhost",
    "Each boundary face of a cell counts as one more neighbour in the undershoot test, with the average of the "
    "cell trace over the face as its mean, so a profile falling into a boundary (a cold wall) is not flagged. "
    "Default false.");
  options.addConfigOption< CFuint >("FreezeAlphaStallWindow",
    "Automatic freeze: once the CFL has not grown for this many iterations and the monitored residual dropped "
    "by less than 0.1 decade over them, alpha is computed from the sensor (without relaxation) on the current "
    "solution and frozen for the rest of the run, together with the subcell reconstruction choices. "
    "The undershoot flags (ForceAlphaMinVars) are kept and new ones still raise alpha to 1. "
    "0 (default): off. Not together with FreezeAlphaIter.");
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSource > >
BaseOrderBlending::providesSockets()
{
  std::vector< SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_alpha);
  result.push_back(&socket_prevAlpha);
  result.push_back(&socket_smoothness);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::execute()
{
  CFTRACEBEGIN;

  SafePtr<vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  DataHandle<CFreal> output = socket_alpha.getDataHandle();
  DataHandle<CFreal> prevAlpha = socket_prevAlpha.getDataHandle();

  const CFuint nbrElemTypes = elemType->size();
  cf_assert(nbrElemTypes == 1);

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // automatic freeze (FreezeAlphaStallWindow), checked once per iteration (this command
  // runs at every nonlinear residual evaluation); the first evaluation only initializes alpha
  bool freezeNow = false;
  if (iter != m_lastIter)
  {
    m_lastIter = iter;
    if (m_freezeStallWindow > 0 && m_alphaInitialized && m_stage == 0) freezeNow = isStalled();
  }

  // FreezeAlphaIter: alpha held from that iteration on
  if (m_alphaInitialized && iter >= getFreezeAlphaIter() && m_stage == 0)
  {
    m_stage = 1;
    getMethodData().setBlendingHoldStage(m_stage);
  }

  // undershoot flags, kept for the whole run (also after the freeze)
  updateForcedCells();

  if (m_stage > 0)
  {
    // held alpha; a newly flagged cell is raised to 1
    for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
    {
      const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
      const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();
      for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
      {
        geoData.idx = elemIdx;
        m_cell = m_cellBuilder->buildGE();
        m_cellStates = m_cell->getStates();
        const CFreal prev = prevAlpha[(*m_cellStates)[0]->getLocalID()];
        const CFreal held = m_forcedCells[elemIdx] ? 1.0 : prev;
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          output[(*m_cellStates)[iSol]->getLocalID()] = held;
          prevAlpha[(*m_cellStates)[iSol]->getLocalID()] = held;
        }
        m_cellBuilder->releaseGE();
      }
    }
  }
  else
  {
    //
    // Phase 1: per-cell physics compute. No neighbor interaction.
    // Writes raw physics-based alpha to socket_alpha. All local cells are
    // processed, including the non-updatable overlap cells of a parallel run,
    // whose states are already synchronised at this point.
    //
    for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
    {
      const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
      const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();

      for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
      {
        geoData.idx = elemIdx;
        m_elemIdx = elemIdx;
        m_cell = m_cellBuilder->buildGE();
        m_cellStates = m_cell->getStates();

        computeSmoothness();

        // a flagged cell gets 1 before the smoothing below, so its neighbours
        // are raised like those of any sensor-detected cell
        const CFreal alpha = m_forcedCells[elemIdx] ? 1.0 : applyAlphaLimits(computeBlendingCoefficient(m_s));
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          output[(*m_cellStates)[iSol]->getLocalID()] = alpha;
        }

        m_cellBuilder->releaseGE();
      }
    }

    //
    // Phase 2: (NbSweeps + 1) Jacobi smoothing iterations.
    // Each iteration snapshots socket_alpha and writes max-pooled values back.
    // The "+1" is the initial neighbor spread; NbSweeps additional passes extend it.
    //
    const CFuint nbIterations = m_nbSweeps + 1;
    for (CFuint it = 0; it < nbIterations; ++it)
    {
      applyJacobiSmoothingPass();
    }

    //
    // Phase 3: relax toward the requested field after all spatial spreading.
    // On the first evaluation and at the automatic freeze the request is used
    // directly. Do not apply the dead-band again after temporal relaxation.
    //
    const bool relax = m_alphaInitialized && !freezeNow && m_alphaRelaxation < 1.0;
    for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
    {
      const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
      const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();
      for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
      {
        geoData.idx = elemIdx;
        m_cell = m_cellBuilder->buildGE();
        m_cellStates = m_cell->getStates();
        const CFuint firstID = (*m_cellStates)[0]->getLocalID();
        CFreal finalAlpha = output[firstID];
        if (m_forcedCells[elemIdx])
        {
          // no relaxation: the flagged cell goes to 1 at once
          finalAlpha = 1.0;
        }
        else if (relax)
        {
          finalAlpha = prevAlpha[firstID] + m_alphaRelaxation * (finalAlpha - prevAlpha[firstID]);
        }
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          output[(*m_cellStates)[iSol]->getLocalID()] = finalAlpha;
          prevAlpha[(*m_cellStates)[iSol]->getLocalID()] = finalAlpha;
        }
        m_cellBuilder->releaseGE();
      }
    }

    m_alphaInitialized = true;

    // this alpha is held from now on; the subcell reconstruction records its
    // choices on the same solution in this iteration
    if (freezeNow)
    {
      m_stage = 1;
      getMethodData().setBlendingHoldStage(m_stage);
    }
  }

  // Say how much blending is actually being applied. Without this there is no
  // way to tell a sensor that never fires from one that is doing its job.
  // Counted over the updatable cells of every rank, so the numbers are global
  // and do not depend on which part of the mesh rank 0 happens to own.
  {
    CFreal aMax = 0., aSum = 0.;
    CFuint nAct = 0, nTot = 0, nForced = 0;
    for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
    {
      const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
      const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();
      for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
      {
        geoData.idx = elemIdx;
        m_cell = m_cellBuilder->buildGE();
        m_cellStates = m_cell->getStates();
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          const CFreal a = output[(*m_cellStates)[0]->getLocalID()];
          aMax = std::max(aMax, a);
          aSum += a;
          if (a > m_alphaMin) ++nAct;
          if (m_forcedCells[elemIdx]) ++nForced;
          ++nTot;
        }
        m_cellBuilder->releaseGE();
      }
    }

    const std::string nsp = getMethodData().getNamespace();
#ifdef CF_HAVE_MPI
    if (PE::GetPE().IsParallel())
    {
      MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
      CFreal aMaxL = aMax, aSumL = aSum;
      CFuint nActL = nAct, nTotL = nTot, nForcedL = nForced;
      MPI_Allreduce(&aMaxL, &aMax, 1, MPI_DOUBLE,   MPI_MAX, comm);
      MPI_Allreduce(&aSumL, &aSum, 1, MPI_DOUBLE,   MPI_SUM, comm);
      MPI_Allreduce(&nActL, &nAct, 1, MPI_UNSIGNED, MPI_SUM, comm);
      MPI_Allreduce(&nTotL, &nTot, 1, MPI_UNSIGNED, MPI_SUM, comm);
      MPI_Allreduce(&nForcedL, &nForced, 1, MPI_UNSIGNED, MPI_SUM, comm);
    }
#endif
    if (nTot > 0 && PE::GetPE().GetRank(nsp) == 0)
    {
      CFLog(INFO, "OrderBlending: alpha max " << aMax
            << ", mean " << aSum/nTot
            << ", " << nAct << "/" << nTot << " cells above AlphaMin");
      if (!m_forceAlphaMinVars.empty())
      {
        CFLog(INFO, ", " << nForced << " forced to 1");
      }
      if (m_stage > 0) CFLog(INFO, ", frozen");
      CFLog(INFO, "\n");
    }
  }

  PE::GetPE().setBarrier(getMethodData().getNamespace());

  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

bool BaseOrderBlending::isStalled()
{
  // monitored residual of the previous iteration and CFL of this one, both global,
  // so every rank takes the same decision
  const CFreal res = SubSystemStatusStack::getActive()->getResidual();
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();

  m_resHistory.push_back(res);
  m_cflHistory.push_back(cfl);
  if (m_resHistory.size() > m_freezeStallWindow + 1)
  {
    m_resHistory.pop_front();
    m_cflHistory.pop_front();
  }
  if (m_resHistory.size() < m_freezeStallWindow + 1) return false;

  // no freeze while the CFL is still growing
  if (m_cflHistory.back() > m_cflHistory.front()) return false;

  // stall: the residual dropped by less than a tenth of a decade over the window
  const CFreal windowMin = *std::min_element(m_resHistory.begin() + 1, m_resHistory.end());
  if (m_resHistory.front() - windowMin >= 0.1) return false;

  const std::string nsp = getMethodData().getNamespace();
  if (PE::GetPE().GetRank(nsp) == 0)
  {
    CFLog(INFO, "OrderBlending: residual stalled over " << m_freezeStallWindow << " iterations at iteration "
          << m_lastIter << " (residual " << res << ", CFL " << cfl
          << "), alpha and the subcell reconstruction choices frozen\n");
  }
  return true;
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::buildBoundaryFaces()
{
  m_bndFacesBuilt = true;

  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  m_bndFaceOrients.assign(cells->getLocalNbGeoEnts(), std::vector<CFuint>());

  // boundary TRSs only: the partition faces have their own start indexes
  std::map< std::string , std::vector< std::vector< CFuint > > >& bndFacesStartIdxsPerTRS =
    getMethodData().getBndFacesStartIdxs();

  SafePtr< GeometricEntityPool< FaceToCellGEBuilder > > faceBuilder = getMethodData().getFaceBuilder();
  FaceToCellGEBuilder::GeoData& faceData = faceBuilder->getDataGE();
  faceData.cellsTRS = cells;
  faceData.isBoundary = true;

  CFuint nbFaces = 0;
  std::map< std::string , std::vector< std::vector< CFuint > > >::const_iterator it;
  for (it = bndFacesStartIdxsPerTRS.begin(); it != bndFacesStartIdxsPerTRS.end(); ++it)
  {
    faceData.facesTRS = MeshDataStack::getActive()->getTrs(it->first);
    const std::vector< std::vector< CFuint > >& startIdxs = it->second;
    for (CFuint iTR = 0; iTR < startIdxs.size(); ++iTR)
    {
      // the orientation of a boundary face is the local index of that face in its cell
      const CFuint nbOrients = startIdxs[iTR].size() - 1;
      for (CFuint orient = 0; orient < nbOrients; ++orient)
      {
        for (CFuint faceID = startIdxs[iTR][orient]; faceID < startIdxs[iTR][orient+1]; ++faceID)
        {
          faceData.idx = faceID;
          GeometricEntity* face = faceBuilder->buildGE();
          m_bndFaceOrients[face->getNeighborGeo(0)->getID()].push_back(orient);
          faceBuilder->releaseGE();
          ++nbFaces;
        }
      }
    }
  }
  CFLog(VERBOSE, "OrderBlending: " << nbFaces << " local boundary faces in the undershoot test\n");
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::computeUndershootReferences()
{
  SafePtr<vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = MeshDataStack::getActive()->getTrs("InnerCells");

  if (m_forceAlphaMinBoundaryGhost && !m_bndFacesBuilt) buildBoundaryFaces();

  const CFuint nbrElemTypes = elemType->size();
  const CFuint nbrVars = m_forceAlphaMinVars.size();

  // cell means of the checked variables, overlap cells included; for a cell with boundary
  // faces also the smallest face mean, where the face mean of a variable is the average of
  // the trace of the cell polynomial over the flux points of that face (the averaging of
  // the cell means, applied to the face)
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();
      m_cellStates = m_cell->getStates();
      for (CFuint iVar = 0; iVar < nbrVars; ++iVar)
      {
        const CFuint var = m_forceAlphaMinVars[iVar];
        CFreal mean = 0.0;
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          mean += (*((*m_cellStates)[iSol]))[var];
        }
        m_minVarsCellMeans[elemIdx][iVar] = mean/m_nbrSolPnts;

        if (m_forceAlphaMinBoundaryGhost)
        {
          CFreal bndMin = MathTools::MathConsts::CFrealMax();
          for (CFuint iFace = 0; iFace < m_bndFaceOrients[elemIdx].size(); ++iFace)
          {
            const std::vector<CFuint>& faceFlxPnts = (*m_faceFlxPntConn)[m_bndFaceOrients[elemIdx][iFace]];
            CFreal faceMean = 0.0;
            for (CFuint iFlx = 0; iFlx < faceFlxPnts.size(); ++iFlx)
            {
              const std::vector<CFreal>& coefs = (*m_solPolyValsAtFlxPnts)[faceFlxPnts[iFlx]];
              for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
              {
                faceMean += coefs[iSol]*(*((*m_cellStates)[iSol]))[var];
              }
            }
            bndMin = std::min(bndMin, faceMean/faceFlxPnts.size());
          }
          m_minVarsBndMeans[elemIdx][iVar] = bndMin;
        }
      }
      m_cellBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::updateForcedCells()
{
  if (m_forceAlphaMinVars.empty())
  {
    return;
  }

  SafePtr<vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();

  computeUndershootReferences();

  // owned unflagged cells, whose neighbour lists are complete
  std::vector<CFuint> newFlags;
  const CFuint nbrElemTypes = elemType->size();
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      if (m_forcedCells[elemIdx] || m_NeighborIDs[elemIdx].empty()) continue;
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();
      m_cellStates = m_cell->getStates();
      if ((*m_cellStates)[0]->isParUpdatable() && hasUndershoot(elemIdx, m_forceAlphaMinMargin))
      {
        m_forcedCells[elemIdx] = true;
        newFlags.push_back((*m_cellStates)[0]->getGlobalID());
      }
      m_cellBuilder->releaseGE();
    }
  }

  // every rank updates its copies of the cells flagged by their owners
  shareForcedCells(newFlags, true);
}

//////////////////////////////////////////////////////////////////////////////

bool BaseOrderBlending::hasUndershoot(const CFuint elemIdx, const CFreal margin)
{
  // a solution point of the current cell below every neighbouring cell mean, and below the
  // face means of the boundary faces with ForceAlphaMinBoundaryGhost, by more than the margin
  const bool ghost = m_forceAlphaMinBoundaryGhost && !m_bndFaceOrients[elemIdx].empty();
  for (CFuint iVar = 0; iVar < m_forceAlphaMinVars.size(); ++iVar)
  {
    const CFuint var = m_forceAlphaMinVars[iVar];
    CFreal neighbourMin = ghost ? m_minVarsBndMeans[elemIdx][iVar] : MathTools::MathConsts::CFrealMax();
    for (CFuint i = 0; i < m_NeighborIDs[elemIdx].size(); ++i)
    {
      neighbourMin = std::min(neighbourMin, m_minVarsCellMeans[m_NeighborIDs[elemIdx][i]][iVar]);
    }
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      if ((*((*m_cellStates)[iSol]))[var] < neighbourMin - margin) return true;
    }
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////

CFuint BaseOrderBlending::shareForcedCells(const std::vector<CFuint>& firstStateGlobalIDs, const bool value)
{
  CFuint nbTotal = firstStateGlobalIDs.size();

#ifdef CF_HAVE_MPI
  if (PE::GetPE().IsParallel())
  {
    const std::string nsp = getMethodData().getNamespace();
    MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
    const int nbRanks = PE::GetPE().GetProcessorCount(nsp);

    int nbLocal = static_cast<int>(firstStateGlobalIDs.size());
    std::vector<int> counts(nbRanks, 0);
    MPI_Allgather(&nbLocal, 1, MPI_INT, &counts[0], 1, MPI_INT, comm);

    std::vector<int> displs(nbRanks, 0);
    for (int r = 1; r < nbRanks; ++r)
    {
      displs[r] = displs[r-1] + counts[r-1];
    }
    nbTotal = displs[nbRanks-1] + counts[nbRanks-1];

    if (nbTotal > 0)
    {
      std::vector<CFuint> allIDs(nbTotal);
      MPI_Allgatherv(firstStateGlobalIDs.empty() ? CFNULL : const_cast<CFuint*>(&firstStateGlobalIDs[0]), nbLocal,
                     MPI_UNSIGNED, &allIDs[0], &counts[0], &displs[0], MPI_UNSIGNED, comm);

      for (CFuint i = 0; i < nbTotal; ++i)
      {
        std::map<CFuint, CFuint>::const_iterator it = m_cellByFirstStateGlobalID.find(allIDs[i]);
        if (it != m_cellByFirstStateGlobalID.end())
        {
          m_forcedCells[it->second] = value;
        }
      }
    }
  }
#endif

  return nbTotal;
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::applyJacobiSmoothingPass()
{
  SafePtr<vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();

  DataHandle<CFreal> output = socket_alpha.getDataHandle();

  // Snapshot current alpha into the scratch buffer before any writes this iteration.
  // This gives true Jacobi updates — neighbor reads are independent of traversal order.
  const CFuint nbStates = output.size();
  cf_assert(m_sweepSnapshot.size() == nbStates);
  for (CFuint i = 0; i < nbStates; ++i)
  {
    m_sweepSnapshot[i] = output[i];
  }

  const CFuint nbrElemTypes = elemType->size();
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();

    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // Scan neighbors first (reading only from the snapshot).
      CFreal alphaNeighborMax = 0.0;
      for (CFuint i = 0; i < m_NeighborIDs[elemIdx].size(); ++i)
      {
        geoData.idx = m_NeighborIDs[elemIdx][i];
        m_cell = m_cellBuilder->buildGE();
        m_cellStates = m_cell->getStates();
        const CFreal alphaN = m_sweepSnapshot[(*m_cellStates)[0]->getLocalID()];
        alphaNeighborMax = std::max(alphaNeighborMax, alphaN);
        m_cellBuilder->releaseGE();
      }

      // Update this cell: origin alpha is preserved (max), neighbor contribution is damped.
      geoData.idx = elemIdx;
      m_elemIdx = elemIdx;
      m_cell = m_cellBuilder->buildGE();
      m_cellStates = m_cell->getStates();

      const CFreal alphaSelf = m_sweepSnapshot[(*m_cellStates)[0]->getLocalID()];
      const CFreal alphaNew  = applyAlphaLimits(
        std::max(alphaSelf, m_neighborWeight * alphaNeighborMax));

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        output[(*m_cellStates)[iSol]->getLocalID()] = alphaNew;
      }

      m_cellBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal BaseOrderBlending::computeBlendingCoefficient(CFreal smoothness) const
{
  // Sinusoidal ramp (Persson-Peraire style) in log-space:
  //   s0 = -S0 * log10(N+1)
  //   alpha = 0                                   if s < s0 - kappa
  //         = 1                                   if s > s0 + kappa
  //         = 0.5 * (1 + sin(pi*(s-s0)/(2*kappa))) otherwise
  const CFreal s0 = -m_s0 * std::log10(static_cast<CFreal>(m_order + 1));

  if (smoothness < s0 - m_kappa)
  {
    return 0.0;
  }
  if (smoothness > s0 + m_kappa)
  {
    return 1.0;
  }
  return 0.5 * (1.0 + std::sin(MathTools::MathConsts::CFrealPi() * (smoothness - s0) / (2.0 * m_kappa)));
}

//////////////////////////////////////////////////////////////////////////////

CFreal BaseOrderBlending::applyAlphaLimits(CFreal alpha) const
{
  if (alpha < m_alphaMin)
  {
    alpha = 0.0;
  }
  // The sinusoidal sensor already reaches one continuously. An early snap
  // to one would introduce a finite jump into the state-dependent residual.
  return std::min(alpha, m_alphaMax);
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::extractMonitoredField()
{
  // Physics-agnostic expressions handled directly. B2 and any other
  // physics-specific expression must be handled by a subclass override.

  if (m_modalMonitoredExpression == "rho")
  {
    // Density: index 0 is universal across variable sets.
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_tempSolPntVec[iSol] = (*((*m_cellStates)[iSol]))[0];
    }
    return;
  }

  if (m_modalMonitoredExpression == "velocity_magnitude")
  {
    // Velocity magnitude from conservative momenta: sqrt((rhoU^2 + rhoV^2 + rhoW^2) / rho^2).
    // Assumes a conservative variable layout with momentum at indices 1..3.
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      const CFreal rho = (*((*m_cellStates)[iSol]))[0];
      const CFreal rhoInv = 1.0 / std::max(std::abs(rho), MathTools::MathConsts::CFrealEps());
      const CFreal u = (*((*m_cellStates)[iSol]))[1] * rhoInv;
      const CFreal v = (*((*m_cellStates)[iSol]))[2] * rhoInv;
      const CFreal w = (m_dim == 3) ? (*((*m_cellStates)[iSol]))[3] * rhoInv : 0.0;
      m_tempSolPntVec[iSol] = std::sqrt(u*u + v*v + w*w);
    }
    return;
  }

  // Pressure-based expressions: use computePhysicalData for physics-agnostic extraction.
  // BaseTerm::P = 1 is the universal pressure slot across all physics models.
  if (m_modalMonitoredExpression == "p" ||
      m_modalMonitoredExpression == "rho*p" ||
      m_modalMonitoredExpression == "p/rho" ||
      m_modalMonitoredExpression == "rho/p")
  {
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_obUpdateVarSet->computePhysicalData(*((*m_cellStates)[iSol]), m_obPData);
      const CFreal p   = m_obPData[1];
      const CFreal rho = (*((*m_cellStates)[iSol]))[0];

      if (m_modalMonitoredExpression == "p")
      {
        m_tempSolPntVec[iSol] = p;
      }
      else if (m_modalMonitoredExpression == "p/rho")
      {
        m_tempSolPntVec[iSol] = p / std::max(std::abs(rho), MathTools::MathConsts::CFrealEps());
      }
      else if (m_modalMonitoredExpression == "rho/p")
      {
        m_tempSolPntVec[iSol] = rho / std::max(std::abs(p), MathTools::MathConsts::CFrealEps());
      }
      else // rho*p
      {
        m_tempSolPntVec[iSol] = rho * p;
      }
    }
    return;
  }

  throw BadValueException(FromHere(),
    "BaseOrderBlending: unknown ModalMonitoredExpression '" + m_modalMonitoredExpression +
    "'. Base class accepts: rho, p, rho*p, p/rho, rho/p, velocity_magnitude. "
    "For B2 or other physics-specific expressions, use a physics-aware subclass "
    "(e.g. OrderBlendingMHD from libFluxReconstructionMHD).");
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::computeSmoothness()
{
  // Step 1: extract the monitored scalar at each solution point (virtual dispatch).
  extractMonitoredField();

  // Step 2: nodal-to-modal transform via inverse Vandermonde.
  m_tempSolPntVec2 = m_vdmInv * m_tempSolPntVec;

  // Step 3: modal energy ratios.
  // E1 = energy in P modes / total energy.
  // E2 = energy in P-1 modes / energy in modes 0..P-1 (guards odd-even decoupling).
  // For P <= 2, E2 would flag physical linear gradients, so only E1 is used.
  const CFreal eps = MathTools::MathConsts::CFrealEps();

  CFreal energyTotal = 0.0;
  CFreal energyP     = 0.0;
  CFreal energyPm1   = 0.0;
  CFreal energyLow   = 0.0;

  for (CFuint j = 0; j < m_nbrSolPnts; ++j)
  {
    const CFreal mj2 = m_tempSolPntVec2[j] * m_tempSolPntVec2[j];
    energyTotal += mj2;

    const CFuint modeOrder = static_cast<CFuint>(m_maxModalOrder[j]);
    if (modeOrder == m_order)
    {
      energyP += mj2;
    }
    else if (modeOrder == m_order - 1)
    {
      energyPm1 += mj2;
    }
    else
    {
      energyLow += mj2;
    }
  }

  const CFreal E1 = energyP / std::max(energyTotal, eps);

  CFreal eVar;
  if (m_order >= 3)
  {
    const CFreal E2 = energyPm1 / std::max(energyLow + energyPm1, eps);
    eVar = std::max(E1, E2);
  }
  else
  {
    eVar = E1;
  }

  // Step 4: log-transform to smoothness indicator.
  m_s = std::log10(std::max(eVar, eps));

  // Store for visualization.
  DataHandle<CFreal> smoothness = socket_smoothness.getDataHandle();
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    smoothness[(*m_cellStates)[iSol]->getLocalID()] = m_s;
  }
}

//////////////////////////////////////////////////////////////////////////////

RealVector BaseOrderBlending::getmaxModalOrder(const CFGeoShape::Type elemShape, const CFuint order)
{
  RealVector maxModalOrder;

  switch (elemShape)
  {
    case CFGeoShape::QUAD:
    {
      // QuadFluxReconstructionElementData stores the modes grouped by shell:
      // the modes of order p occupy the range [p^2, (p+1)^2).
      const CFuint totalModes = (order + 1) * (order + 1);
      maxModalOrder.resize(totalModes);
      CFuint modeIndex = 0;
      for (CFuint p = 0; p <= order; ++p)
      {
        for (CFuint i = p*p; i < (p+1)*(p+1); ++i)
        {
          maxModalOrder[modeIndex++] = p;
        }
      }
      cf_assert(modeIndex == totalModes);
      break;
    }
    case CFGeoShape::TRIAG:
    {
      const CFuint totalModes = ((order + 1) * (order + 2)) / 2;
      maxModalOrder.resize(totalModes);
      CFuint modeIndex = 0;
      for (CFuint totalOrder = 0; totalOrder <= order; ++totalOrder)
      {
        for (CFuint iOrderKsi = totalOrder;; --iOrderKsi)
        {
          const CFuint iOrderEta = totalOrder - iOrderKsi;
          maxModalOrder[modeIndex++] = std::max(iOrderKsi, iOrderEta);
          if (iOrderKsi == 0) break;
        }
      }
      cf_assert(modeIndex == totalModes);
      break;
    }
    case CFGeoShape::TETRA:
    {
      const CFuint totalModes = ((order + 1) * (order + 2) * (order + 3)) / 6;
      maxModalOrder.resize(totalModes);
      CFuint modeIndex = 0;
      for (CFuint totalOrder = 0; totalOrder <= order; ++totalOrder)
      {
        for (CFuint iOrderZta = 0; iOrderZta <= totalOrder; ++iOrderZta)
        {
          for (CFuint iOrderEta = 0; iOrderEta + iOrderZta <= totalOrder; ++iOrderEta)
          {
            maxModalOrder[modeIndex++] = totalOrder;
          }
        }
      }
      cf_assert(modeIndex == totalModes);
      break;
    }
    case CFGeoShape::PRISM:
    {
      const CFuint totalModes = ((order + 1) * (order + 1) * (order + 2)) / 2;
      maxModalOrder.resize(totalModes);
      CFuint modeIndex = 0;
      for (CFuint totalOrderXY = 0; totalOrderXY <= order; ++totalOrderXY)
      {
        for (CFuint iOrderZta = 0; iOrderZta <= order; ++iOrderZta)
        {
          for (CFuint iOrderKsi = totalOrderXY;; --iOrderKsi)
          {
            maxModalOrder[modeIndex++] = std::max(totalOrderXY, iOrderZta);
            if (iOrderKsi == 0) break;
          }
        }
      }
      cf_assert(modeIndex == totalModes);
      break;
    }
    case CFGeoShape::HEXA:
    {
      // HexaFluxReconstructionElementData stores the modes grouped by shell:
      // it walks (iKsi, iEta, iZta) but writes each mode at column
      // max(iKsi,iEta,iZta)^3 + counter, so the modes of order p occupy the
      // contiguous range [p^3, (p+1)^3). Walking the triple loop in its own
      // order instead would mislabel the modes from P2 upward.
      const CFuint totalModes = (order + 1) * (order + 1) * (order + 1);
      maxModalOrder.resize(totalModes);
      CFuint modeIndex = 0;
      for (CFuint p = 0; p <= order; ++p)
      {
        for (CFuint i = p*p*p; i < (p+1)*(p+1)*(p+1); ++i)
        {
          maxModalOrder[modeIndex++] = p;
        }
      }
      cf_assert(modeIndex == totalModes);
      break;
    }
    default:
      throw Common::NotImplementedException(FromHere(),
        "BaseOrderBlending::getmaxModalOrder: unsupported element shape " +
        StringOps::to_str(elemShape));
  }

  return maxModalOrder;
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::setup()
{
  CFAUTOTRACE;

  m_alphaInitialized = false;

  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim    = PhysicalModelStack::getActive()->getDim();

  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

  vector<FluxReconstructionElementData*>& frLocalData = getMethodData().getFRLocalData();
  const CFuint nbrElemTypes = frLocalData.size();
  cf_assert(nbrElemTypes > 0);

  m_order      = static_cast<CFuint>(frLocalData[0]->getPolyOrder());
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();
  const CFGeoShape::Type elemShape = frLocalData[0]->getShape();

  m_vdmInv = *(frLocalData[0]->getVandermondeMatrixInv());

  m_tempSolPntVec.resize(m_nbrSolPnts);
  m_tempSolPntVec2.resize(m_nbrSolPnts);
  m_maxModalOrder.resize(m_nbrSolPnts);
  m_maxModalOrder = getmaxModalOrder(elemShape, m_order);

  // Physical data vector for pressure-based expressions.
  m_obUpdateVarSet = getMethodData().getUpdateVar();
  SafePtr<BaseTerm> convTerm = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm();
  convTerm->resizePhysicalData(m_obPData);

  // Allocate per-state sockets.
  SafePtr<vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbrCells = (*elemType)[0].getEndIdx();
  const CFuint nbStates = nbrCells * m_nbrSolPnts;

  socket_alpha.getDataHandle().resize(nbStates);
  socket_prevAlpha.getDataHandle().resize(nbStates);
  socket_smoothness.getDataHandle().resize(nbStates);

  // undershoot flags (ForceAlphaMinVars) and the means they are tested against
  m_forcedCells.assign(nbrCells, false);
  m_bndFacesBuilt = false;
  m_faceFlxPntConn = frLocalData[0]->getFaceFlxPntConn();
  m_solPolyValsAtFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();

  // alpha follows the sensor at the start of every run
  m_stage = 0;
  m_lastIter = std::numeric_limits<CFuint>::max();
  m_resHistory.clear();
  m_cflHistory.clear();
  getMethodData().setBlendingHoldStage(0);

  for (CFuint i = 0; i < m_forceAlphaMinVars.size(); ++i)
  {
    if (m_forceAlphaMinVars[i] >= m_nbrEqs)
    {
      throw BadValueException(FromHere(), "OrderBlending: ForceAlphaMinVars index outside the state");
    }
  }
  if (!(m_forceAlphaMinMargin >= 0.0))
  {
    throw BadValueException(FromHere(), "OrderBlending: ForceAlphaMinMargin must be >= 0");
  }
  m_minVarsCellMeans.assign(nbrCells, std::vector<CFreal>(m_forceAlphaMinVars.size(), 0.0));
  m_minVarsBndMeans.assign(nbrCells, std::vector<CFreal>(m_forceAlphaMinVars.size(), 0.0));

  m_cellByFirstStateGlobalID.clear();
  if (!m_forceAlphaMinVars.empty())
  {
    SafePtr<TopologicalRegionSet> innerCells = MeshDataStack::getActive()->getTrs("InnerCells");
    StdTrsGeoBuilder::GeoData& geoDataFlags = m_cellBuilder->getDataGE();
    geoDataFlags.trs = innerCells;
    for (CFuint elemIdx = 0; elemIdx < nbrCells; ++elemIdx)
    {
      geoDataFlags.idx = elemIdx;
      GeometricEntity* cell = m_cellBuilder->buildGE();
      m_cellByFirstStateGlobalID[(*cell->getStates())[0]->getGlobalID()] = elemIdx;
      m_cellBuilder->releaseGE();
    }
  }

  // Scratch buffer for Jacobi smoothing passes.
  m_sweepSnapshot.assign(nbStates, 0.0);

  // Pre-compute node-sharing neighbor IDs for each cell.
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  m_NeighborIDs.resize(nbrCells);

  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();

    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      std::set<CFuint> currNodeIDs;
      const CFuint nbNodesCurr = cells->getNbNodesInGeo(elemIdx);
      for (CFuint iNode = 0; iNode < nbNodesCurr; ++iNode)
      {
        currNodeIDs.insert(cells->getNodeID(elemIdx, iNode));
      }

      std::set<CFuint> neighborSet;
      for (CFuint nIdx = startIdx; nIdx < endIdx; ++nIdx)
      {
        if (nIdx == elemIdx) continue;

        const CFuint nbNodesNeighbor = cells->getNbNodesInGeo(nIdx);
        for (CFuint jNode = 0; jNode < nbNodesNeighbor; ++jNode)
        {
          if (currNodeIDs.find(cells->getNodeID(nIdx, jNode)) != currNodeIDs.end())
          {
            neighborSet.insert(nIdx);
            break;
          }
        }
      }

      m_NeighborIDs[elemIdx].assign(neighborSet.begin(), neighborSet.end());
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
