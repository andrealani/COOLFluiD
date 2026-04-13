// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"
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
  m_order(0),
  m_nbrSolPnts(0),
  m_nbrEqs(0),
  m_dim(0),
  m_iElemType(0),
  m_elemIdx(0)
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

  // Decay factor for neighbor alpha during Jacobi smoothing.
  m_neighborWeight = 0.5;
  setParameter("NeighborWeight", &m_neighborWeight);

  // Smoothing iterations beyond the initial spread: total passes = m_nbSweeps + 1.
  m_nbSweeps = 0;
  setParameter("NbSweeps", &m_nbSweeps);

  // Iteration to freeze alpha (reuses prevAlpha). Default: never freeze.
  m_freezeFilterIter = 1000000;
  setParameter("freezeFilterIter", &m_freezeFilterIter);
}

//////////////////////////////////////////////////////////////////////////////

BaseOrderBlending::~BaseOrderBlending()
{
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::configure(Config::ConfigArgs& args)
{
  FluxReconstructionSolverCom::configure(args);
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
    "Dead-band threshold: alpha < AlphaMin snaps to 0, alpha > 1-AlphaMin snaps to 1.");
  options.addConfigOption< CFreal >("AlphaMax",
    "Maximum blending coefficient cap.");
  options.addConfigOption< CFreal >("NeighborWeight",
    "Decay factor applied to neighbor alpha during Jacobi smoothing. "
    "0 disables spreading.");
  options.addConfigOption< CFuint >("NbSweeps",
    "Number of smoothing iterations beyond the initial spread. "
    "Total spreading passes = NbSweeps + 1.");
  options.addConfigOption< CFuint, Config::DynamicOption<> >("freezeFilterIter",
    "Iteration number at which alpha is frozen (reuses prevAlpha). "
    "Very large = never freeze.");
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

  // Freeze path: reuse prevAlpha without recomputing anything.
  if (iter >= m_freezeFilterIter)
  {
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
          const CFreal frozen = prevAlpha[(*m_cellStates)[0]->getLocalID()];
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            output[(*m_cellStates)[iSol]->getLocalID()] = frozen;
          }
        }
        m_cellBuilder->releaseGE();
      }
    }
    PE::GetPE().setBarrier(getMethodData().getNamespace());
    CFTRACEEND;
    return;
  }

  //
  // Phase 1: per-cell physics compute. No neighbor interaction.
  // Writes raw physics-based alpha to socket_alpha.
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

      if ((*m_cellStates)[0]->isParUpdatable())
      {
        computeSmoothness();

        const CFreal alpha = applyAlphaLimits(computeBlendingCoefficient(m_s));
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          output[(*m_cellStates)[iSol]->getLocalID()] = alpha;
        }
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
  // Phase 3: snapshot final alpha into prevAlpha for the freeze mechanism.
  //
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
        const CFreal finalAlpha = output[(*m_cellStates)[0]->getLocalID()];
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          prevAlpha[(*m_cellStates)[iSol]->getLocalID()] = finalAlpha;
        }
      }
      m_cellBuilder->releaseGE();
    }
  }

  PE::GetPE().setBarrier(getMethodData().getNamespace());

  CFTRACEEND;
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
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          const CFreal alphaN = m_sweepSnapshot[(*m_cellStates)[0]->getLocalID()];
          alphaNeighborMax = std::max(alphaNeighborMax, alphaN);
        }
        m_cellBuilder->releaseGE();
      }

      // Update this cell: origin alpha is preserved (max), neighbor contribution is damped.
      geoData.idx = elemIdx;
      m_elemIdx = elemIdx;
      m_cell = m_cellBuilder->buildGE();
      m_cellStates = m_cell->getStates();

      if ((*m_cellStates)[0]->isParUpdatable())
      {
        const CFreal alphaSelf = m_sweepSnapshot[(*m_cellStates)[0]->getLocalID()];
        const CFreal alphaNew  = applyAlphaLimits(
          std::max(alphaSelf, m_neighborWeight * alphaNeighborMax));

        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          output[(*m_cellStates)[iSol]->getLocalID()] = alphaNew;
        }
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
  else if (alpha > 1.0 - m_alphaMin)
  {
    alpha = 1.0;
  }
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
      const CFuint totalModes = (order + 1) * (order + 1) * (order + 1);
      maxModalOrder.resize(totalModes);
      CFuint modeIndex = 0;
      for (CFuint iOrderKsi = 0; iOrderKsi <= order; ++iOrderKsi)
      {
        for (CFuint iOrderEta = 0; iOrderEta <= order; ++iOrderEta)
        {
          for (CFuint iOrderZta = 0; iOrderZta <= order; ++iOrderZta)
          {
            maxModalOrder[modeIndex++] = std::max({iOrderKsi, iOrderEta, iOrderZta});
          }
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
