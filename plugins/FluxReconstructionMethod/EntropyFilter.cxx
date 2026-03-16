#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/BaseTerm.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/EntropyFilter.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/BaseOrderBlending.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<EntropyFilter,
                      FluxReconstructionSolverData,
                      FluxReconstructionModule>
    entropyFilterProvider("EntropyFilter");

//////////////////////////////////////////////////////////////////////////////

void EntropyFilter::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Gamma","Heat capacity ratio (default 1.4).");
  options.addConfigOption< CFreal >("MinDensity","Minimum density positivity floor (default 1e-6).");
  options.addConfigOption< CFreal >("MinPressure","Minimum pressure positivity floor (default 1e-6).");
  options.addConfigOption< CFreal >("ZetaMax","Maximum filter strength (default 20.0).");
  options.addConfigOption< CFuint >("NBisect","Maximum root-bracketing iterations (default 20).");
  options.addConfigOption< CFreal >("Tolerance","Root-bracketing stopping tolerance on f (default 1e-4).");
  options.addConfigOption< CFreal >("EntropyTol","Entropy constraint tolerance: s >= sMin - eTol (default 1e-6).");
  options.addConfigOption< CFuint >("ShowRate","Print frequency for filter statistics (default 1).");
}

//////////////////////////////////////////////////////////////////////////////

EntropyFilter::EntropyFilter(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_filterStrength("filterStrength"),
  m_cellBuilder(CFNULL),
  m_cell(),
  m_cellStates(),
  m_nbrEqs(),
  m_nbrSolPnts(),
  m_order(),
  m_nbrModeOrders(),
  m_cellAvgSolCoefs(CFNULL),
  m_tmpState(CFNULL),
  m_tmpCoord(CFNULL),
  m_updateVarSet(CFNULL),
  m_nbrFlxPnts(0),
  m_nbrEvalPnts(0),
  m_solPolyInFlxPnts(CFNULL),
  m_faceBuilder(CFNULL),
  m_totalCells(0),
  m_nbFiltered(0),
  m_totalNbFiltered(0)
{
  addConfigOptionsTo(this);

  m_gamma = 1.4;
  setParameter("Gamma", &m_gamma);

  m_minDensity = 1e-6;
  setParameter("MinDensity", &m_minDensity);

  m_minPressure = 1e-6;
  setParameter("MinPressure", &m_minPressure);

  m_zetaMax = 20.0;
  setParameter("ZetaMax", &m_zetaMax);

  m_nBisect = 20;
  setParameter("NBisect", &m_nBisect);

  m_tolerance = 1e-4;
  setParameter("Tolerance", &m_tolerance);

  m_eTol = 1e-6;
  setParameter("EntropyTol", &m_eTol);

  m_showrate = 1;
  setParameter("ShowRate", &m_showrate);

}

//////////////////////////////////////////////////////////////////////////////

EntropyFilter::~EntropyFilter()
{
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilter::configure(Config::ConfigArgs& args)
{
  FluxReconstructionSolverCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
EntropyFilter::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_filterStrength);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilter::setup()
{
  CFAUTOTRACE;

  // get number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();

  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);

  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  m_order = static_cast<CFuint>(order);

  // get nbr of sol pnts
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();

  // get nbr of flux pnts
  m_nbrFlxPnts = frLocalData[0]->getNbrOfFlxPnts();

  // total evaluation points = sol + flx
  m_nbrEvalPnts = m_nbrSolPnts + m_nbrFlxPnts;

  // get Vandermonde matrices (SafePtr — owned by frLocalData, persistent)
  m_vdm = frLocalData[0]->getVandermondeMatrix();
  m_vdmInv = frLocalData[0]->getVandermondeMatrixInv();

  // get cell-averaged solution coefficients
  m_cellAvgSolCoefs = frLocalData[0]->getCellAvgSolCoefs();

  // get interpolation from sol points to flux points
  m_solPolyInFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();

  // get modal orders using BaseOrderBlending's static method
  const CFGeoShape::Type elemShape = frLocalData[0]->getShape();
  const RealVector tmpModalOrder = BaseOrderBlending::getmaxModalOrder(elemShape, m_order);
  m_maxModalOrder.resize(tmpModalOrder.size());
  for (CFuint i = 0; i < tmpModalOrder.size(); ++i)
  {
    m_maxModalOrder[i] = static_cast<CFuint>(tmpModalOrder[i]);
  }

  // determine number of distinct mode orders
  m_nbrModeOrders = 0;
  for (CFuint i = 0; i < m_maxModalOrder.size(); ++i)
  {
    if (m_maxModalOrder[i] + 1 > m_nbrModeOrders) m_nbrModeOrders = m_maxModalOrder[i] + 1;
  }

  // resize mode contributions for sol points only (used by computeFilteredStates)
  m_modeContrib.resize(m_nbrModeOrders);
  for (CFuint k = 0; k < m_nbrModeOrders; ++k)
  {
    m_modeContrib[k].resize(m_nbrSolPnts);
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_modeContrib[k][iSol].resize(m_nbrEqs);
    }
  }

  // resize mode contributions for ALL evaluation points (sol + flx)
  // m_modeContribPt[k][iPt] for iPt = 0..nbrEvalPnts-1
  m_modeContribPt.resize(m_nbrModeOrders);
  for (CFuint k = 0; k < m_nbrModeOrders; ++k)
  {
    m_modeContribPt[k].resize(m_nbrEvalPnts);
    for (CFuint iPt = 0; iPt < m_nbrEvalPnts; ++iPt)
    {
      m_modeContribPt[k][iPt].resize(m_nbrEqs);
    }
  }

  // resize filtered states (sol points only — for final write-back)
  m_filteredStates.resize(m_nbrSolPnts);
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_filteredStates[iSol].resize(m_nbrEqs);
  }

  // resize single-point filtered state
  m_filteredPtState.resize(m_nbrEqs);

  // resize cell-averaged state
  m_cellAvgState.resize(m_nbrEqs);

  // allocate temporary state for physics queries
  m_tmpState = new State();
  m_tmpState->resize(m_nbrEqs);

  // allocate dummy coordinates for m_tmpState (some variable sets
  // call getCoordinates() inside computePhysicalData)
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  m_tmpCoord = new Node(false);
  m_tmpCoord->resize(dim);
  *m_tmpCoord = 0.0;
  m_tmpState->setSpaceCoordinates(m_tmpCoord);

  // resize physical data
  PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->resizePhysicalData(m_physData);

  // get update variable set
  m_updateVarSet = getMethodData().getUpdateVar();

  // resize modal coefficients working vectors
  m_modalCoeffs.resize(m_nbrSolPnts);
  m_uHat.resize(m_nbrSolPnts);

  // get face builder for neighbor entropy exchange
  m_faceBuilder = getMethodData().getFaceBuilder();

  // resize filter strength socket — count ALL cells across all element types
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  m_totalCells = 0;
  for (CFuint iType = 0; iType < elemType->size(); ++iType)
  {
    m_totalCells += (*elemType)[iType].getNbElems();
  }
  DataHandle< CFreal > filterStrength = socket_filterStrength.getDataHandle();
  filterStrength.resize(m_totalCells * m_nbrSolPnts);
  filterStrength = 0.0;

  // resize per-cell entropy arrays for neighbor exchange
  m_localMinEntropy.resize(m_totalCells, MathTools::MathConsts::CFrealMax());
  m_neighborMinEntropy.resize(m_totalCells, MathTools::MathConsts::CFrealMax());
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilter::unsetup()
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

bool EntropyFilter::computePhysics(const RealVector& state, CFreal& rho, CFreal& p, CFreal& s)
{
  m_tmpState->copyData(state);
  m_updateVarSet->computePhysicalData(*m_tmpState, m_physData);
  rho = m_physData[0];
  p   = m_physData[1];

  if (rho > 0.0 && p > 0.0)
  {
    s = p / std::pow(rho, m_gamma);
    return true;
  }
  else
  {
    s = MathTools::MathConsts::CFrealMax();
    return false;
  }
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilter::execute()
{
  CFTRACEBEGIN;

  // P0 has no higher modes to filter
  if (m_order == 0)
  {
    CFTRACEEND;
    return;
  }

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  const CFuint nbrElemTypes = elemType->size();
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // get filter strength data handle for output
  DataHandle< CFreal > filterStrength = socket_filterStrength.getDataHandle();

  // reset filter count
  m_nbFiltered = 0;
  m_totalNbFiltered = 0;

  // --- Neighbor entropy exchange ---
  // m_localMinEntropy persists across iterations: it contains the post-filter
  // entropy from the PREVIOUS call (initialized to CFrealMax in setup).
  // This is critical: the current solution may have new oscillations from
  // the time stepper that violate the previous entropy bound.
  // On the first call, CFrealMax disables the entropy constraint (positivity only).
  exchangeNeighborEntropy();

  // --- Pass 2: Filter cells where constraints are violated ---
  // loop over all element types
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // Entropy bound from neighbor exchange: minimum over cell AND face-neighbors
      const CFreal sMin = m_neighborMinEntropy[elemIdx];

      // Compute minimum density, pressure, and entropy across all evaluation points
      CFreal dMin = MathTools::MathConsts::CFrealMax();
      CFreal pMin = MathTools::MathConsts::CFrealMax();
      CFreal eMin = MathTools::MathConsts::CFrealMax();

      // Check at solution points
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        CFreal rho, p, s;
        computePhysics(*((*m_cellStates)[iSol]), rho, p, s);
        dMin = std::min(dMin, rho);
        pMin = std::min(pMin, p);
        eMin = std::min(eMin, s);
      }

      // Check at flux points (interpolate from sol points)
      for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
      {
        m_filteredPtState = 0.0;
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          const CFreal coef = (*m_solPolyInFlxPnts)[iFlx][iSol];
          for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
          {
            m_filteredPtState[iEq] += coef * (*((*m_cellStates)[iSol]))[iEq];
          }
        }

        CFreal rho, p, s;
        computePhysics(m_filteredPtState, rho, p, s);
        dMin = std::min(dMin, rho);
        pMin = std::min(pMin, p);
        eMin = std::min(eMin, s);
      }

      // Check if entropy bound is finite for the entropy trigger
      const bool finiteEntBound = (sMin < MathTools::MathConsts::CFrealMax() * 0.5);

      // Trigger filtering if ANY constraint is violated at ANY evaluation point:
      // (1) density below floor, (2) pressure below floor, or
      // (3) entropy below cell-average entropy (minus tolerance)
      const bool entropyViolated = finiteEntBound && (eMin < sMin - m_eTol);
      if (dMin >= m_minDensity && pMin >= m_minPressure && !entropyViolated)
      {
        // No filtering needed — update entropy for next iteration
        m_localMinEntropy[elemIdx] = eMin;
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          filterStrength[((*m_cellStates)[iSol])->getLocalID()] = 0.0;
        }
        m_cellBuilder->releaseGE();
        continue;
      }

      // --- Step 2: Pre-decompose solution by mode order ---
      // For each equation, compute modal coefficients and group by order.
      // Build mode contributions at ALL evaluation points (sol + flx).

      // Initialize all mode contributions to zero
      for (CFuint k = 0; k < m_nbrModeOrders; ++k)
      {
        for (CFuint iPt = 0; iPt < m_nbrEvalPnts; ++iPt)
        {
          m_modeContribPt[k][iPt] = 0.0;
        }
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          m_modeContrib[k][iSol] = 0.0;
        }
      }

      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        // Extract nodal values for this equation
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          m_modalCoeffs[iSol] = (*((*m_cellStates)[iSol]))[iEq];
        }

        // Transform to modal space: u_hat = V^{-1} * u
        m_uHat = 0.0;
        for (CFuint i = 0; i < m_nbrSolPnts; ++i)
        {
          for (CFuint j = 0; j < m_nbrSolPnts; ++j)
          {
            m_uHat[i] += (*m_vdmInv)(i, j) * m_modalCoeffs[j];
          }
        }

        // For each mode, accumulate its contribution grouped by mode order k
        for (CFuint iMode = 0; iMode < m_nbrSolPnts; ++iMode)
        {
          const CFuint k = m_maxModalOrder[iMode];

          // Contribution at solution points: V(iSol, iMode) * uHat[iMode]
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            const CFreal contrib = (*m_vdm)(iSol, iMode) * m_uHat[iMode];
            m_modeContrib[k][iSol][iEq] += contrib;
            m_modeContribPt[k][iSol][iEq] += contrib;
          }

          // Contribution at flux points: interpolated from sol points
          // V_flx(iFlx, iMode) = sum_j m_solPolyInFlxPnts[iFlx][j] * V(j, iMode)
          for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
          {
            CFreal vFlx = 0.0;
            for (CFuint j = 0; j < m_nbrSolPnts; ++j)
            {
              vFlx += (*m_solPolyInFlxPnts)[iFlx][j] * (*m_vdm)(j, iMode);
            }
            m_modeContribPt[k][m_nbrSolPnts + iFlx][iEq] += vFlx * m_uHat[iMode];
          }
        }
      }

      // --- Step 3: Per-point progressive tightening ---
      // Start with f = 1 (no filter). For each evaluation point, check if
      // constraints are satisfied with the current f. If not, bisect between
      // [0, f] to find a smaller f that satisfies constraints at this point.
      // After processing all points, f is the global minimum.

      CFreal f = 1.0;

      for (CFuint iPt = 0; iPt < m_nbrEvalPnts; ++iPt)
      {
        // Compute filtered state at this point with current f
        computeFilteredStateAtPoint(iPt, f, m_filteredPtState);

        CFreal rho, p, s;
        computePhysics(m_filteredPtState, rho, p, s);

        // Check if constraints are violated at this point
        const bool ptEntropyBad = finiteEntBound && (s < sMin - m_eTol);
        if (rho < m_minDensity || p < m_minPressure || ptEntropyBad)
        {
          // Bisect between [0, f] to find tighter bound
          CFreal fHi = f;
          CFreal fLo = 0.0;

          for (CFuint bisectIter = 0; bisectIter < m_nBisect; ++bisectIter)
          {
            if (fHi - fLo < m_tolerance) break;

            const CFreal fMid = 0.5 * (fLo + fHi);

            computeFilteredStateAtPoint(iPt, fMid, m_filteredPtState);
            computePhysics(m_filteredPtState, rho, p, s);

            const bool midEntropyBad = finiteEntBound && (s < sMin - m_eTol);
            if (rho < m_minDensity || p < m_minPressure || midEntropyBad)
            {
              fHi = fMid;  // still violating — need more filtering
            }
            else
            {
              fLo = fMid;  // constraints satisfied — try less filtering
            }
          }

          // Set current minimum f as the bounds-preserving value
          f = fLo;
        }
      }

      // --- Step 4: Apply filter to all solution points with the found f ---
      computeFilteredStates(f);

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        (*((*m_cellStates)[iSol])) = m_filteredStates[iSol];
        // store filter strength: 0 = no filter (f=1), 1 = max filter (f=0)
        filterStrength[((*m_cellStates)[iSol])->getLocalID()] = 1.0 - f;
      }

      if (f < 1.0 && (*m_cellStates)[0]->isParUpdatable()) ++m_nbFiltered;

      // --- Step 5: Update local min entropy from filtered solution ---
      // Recompute eMin from the filtered states so that entropy bounds
      // reflect the post-filter state (needed for multi-stage time integrators)
      {
        CFreal eMinFiltered = MathTools::MathConsts::CFrealMax();
        // Check at sol points (already updated in-place)
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          CFreal rho, p, s;
          computePhysics(*((*m_cellStates)[iSol]), rho, p, s);
          eMinFiltered = std::min(eMinFiltered, s);
        }
        // Check at flux points (interpolated from filtered sol points)
        for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
        {
          m_filteredPtState = 0.0;
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            const CFreal coef = (*m_solPolyInFlxPnts)[iFlx][iSol];
            for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
            {
              m_filteredPtState[iEq] += coef * (*((*m_cellStates)[iSol]))[iEq];
            }
          }
          CFreal rho, p, s;
          computePhysics(m_filteredPtState, rho, p, s);
          eMinFiltered = std::min(eMinFiltered, s);
        }
        m_localMinEntropy[elemIdx] = eMinFiltered;
      }

      // release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }

  // MPI reduction and logging
  const std::string nsp = this->getMethodData().getNamespace();

#ifdef CF_HAVE_MPI
  MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
  PE::GetPE().setBarrier(nsp);
  const CFuint count = 1;
  MPI_Allreduce(&m_nbFiltered, &m_totalNbFiltered, count, MPI_UNSIGNED, MPI_SUM, comm);
#else
  m_totalNbFiltered = m_nbFiltered;
#endif

  if (PE::GetPE().GetRank(nsp) == 0 && iter % m_showrate == 0)
  {
    CFLog(NOTICE, "EntropyFilter: " << m_totalNbFiltered << " cells filtered\n");
  }

  PE::GetPE().setBarrier(nsp);

  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilter::computeFilteredStates(const CFreal f)
{
  // Filtered state at sol points only:
  // u_tilde_j[iEq] = sum_{k=0}^{P} f^(k^2) * m_modeContrib[k][j][iEq]
  // Note: f^(0^2) = f^0 = 1 always — cell average is preserved exactly
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_filteredStates[iSol] = m_modeContrib[0][iSol]; // k=0: f^0 = 1

    for (CFuint k = 1; k < m_nbrModeOrders; ++k)
    {
      const CFreal fPow = std::pow(f, static_cast<CFreal>(k * k));
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        m_filteredStates[iSol][iEq] += fPow * m_modeContrib[k][iSol][iEq];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilter::computeFilteredStateAtPoint(const CFuint iPt, const CFreal f, RealVector& filteredState)
{
  // Filtered state at a single evaluation point:
  // u_tilde[iEq] = sum_{k=0}^{P} f^(k^2) * m_modeContribPt[k][iPt][iEq]
  filteredState = m_modeContribPt[0][iPt]; // k=0: always preserved

  for (CFuint k = 1; k < m_nbrModeOrders; ++k)
  {
    const CFreal fPow = std::pow(f, static_cast<CFreal>(k * k));
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      filteredState[iEq] += fPow * m_modeContribPt[k][iPt][iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilter::computeLocalMinEntropy()
{
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  const CFuint nbrElemTypes = elemType->size();

  // loop over all element types
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();

    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();
      m_cellStates = m_cell->getStates();

      CFreal eMin = MathTools::MathConsts::CFrealMax();

      // Entropy at solution points
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        CFreal rho, p, s;
        computePhysics(*((*m_cellStates)[iSol]), rho, p, s);
        eMin = std::min(eMin, s);
      }

      // Entropy at flux points (interpolated from sol points)
      for (CFuint iFlx = 0; iFlx < m_nbrFlxPnts; ++iFlx)
      {
        m_filteredPtState = 0.0;
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          const CFreal coef = (*m_solPolyInFlxPnts)[iFlx][iSol];
          for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
          {
            m_filteredPtState[iEq] += coef * (*((*m_cellStates)[iSol]))[iEq];
          }
        }

        CFreal rho, p, s;
        computePhysics(m_filteredPtState, rho, p, s);
        eMin = std::min(eMin, s);
      }

      m_localMinEntropy[elemIdx] = eMin;

      m_cellBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void EntropyFilter::exchangeNeighborEntropy()
{
  // Initialize neighbor entropy from local entropy
  m_neighborMinEntropy = m_localMinEntropy;

  // get InnerCells and InnerFaces TRS
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  // get the face start indexes
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();

  // get number of face orientations
  const CFuint nbrFaceOrients = innerFacesStartIdxs.size() - 1;

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoDataFace = m_faceBuilder->getDataGE();
  geoDataFace.cellsTRS = cells;
  geoDataFace.facesTRS = faces;
  geoDataFace.isBoundary = false;

  // loop over different orientations
  for (CFuint orient = 0; orient < nbrFaceOrients; ++orient)
  {
    const CFuint faceStartIdx = innerFacesStartIdxs[orient];
    const CFuint faceStopIdx  = innerFacesStartIdxs[orient + 1];

    // loop over faces with this orientation
    for (CFuint faceID = faceStartIdx; faceID < faceStopIdx; ++faceID)
    {
      // build the face GeometricEntity
      geoDataFace.idx = faceID;
      Framework::GeometricEntity* face = m_faceBuilder->buildGE();

      // get the neighbouring cell IDs
      const CFuint leftIdx  = face->getNeighborGeo(LEFT)->getID();
      const CFuint rightIdx = face->getNeighborGeo(RIGHT)->getID();

      // exchange: each side gets min of both local entropies
      // Read from m_localMinEntropy (immutable) to avoid transitive propagation
      m_neighborMinEntropy[leftIdx]  = std::min(m_neighborMinEntropy[leftIdx],  m_localMinEntropy[rightIdx]);
      m_neighborMinEntropy[rightIdx] = std::min(m_neighborMinEntropy[rightIdx], m_localMinEntropy[leftIdx]);

      m_faceBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
