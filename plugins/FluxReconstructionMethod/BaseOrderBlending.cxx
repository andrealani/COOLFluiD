// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/BaseOrderBlending.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<BaseOrderBlending, FluxReconstructionSolverData, FluxReconstructionModule>
    OrderBlendingFRProvider("OrderBlending");
//////////////////////////////////////////////////////////////////////////////

BaseOrderBlending::BaseOrderBlending(const std::string& name) :
  FluxReconstructionSolverCom(name),
  m_showrate(),
  m_freezeFilterIter(),
  m_cellBuilder(CFNULL),
  m_iElemType(),
  m_cell(),
  m_cellStates(),
  m_nbrEqs(),
  m_dim(),
  m_nbrSolPnts(),
  m_maxNbrFlxPnts(),
  m_order(),
  m_vdm(),
  m_vdmInv(),
  m_solPolyValsAtFlxPnts(CFNULL),
  m_cellStatesFlxPnt(),
  m_nbLimits(),
  m_totalNbLimits(),
  m_NeighborIDs(),
  m_tempSolPntVec(),
  m_tempSolPntVec2(),
  socket_prevAlpha("prevAlpha"),
  socket_alpha("alpha"),
  socket_smoothness("smoothness")
{
  addConfigOptionsTo(this);
  
  //==========================================================================
  // MAIN PARAMETERS (Log-Scale Method - RECOMMENDED)
  //==========================================================================
  
  // Use log-scale method (recommended) or legacy sigmoid
  m_useLogScale = true;
  setParameter( "UseLogScale", &m_useLogScale );
  
  // Reference smoothness: s0 = -S0 * log10(N+1)
  // Higher S0 = more dissipation. Typical values: 3.0 - 5.0
  m_s0 = 4.0;
  setParameter( "S0", &m_s0);
  
  // Transition half-width in log-space
  // Blending ramps over [s0 - Kappa, s0 + Kappa]
  // Typical values: 1.0 - 2.0
  m_kappa = 1.5;
  setParameter( "Kappa", &m_kappa);
  
  // Monitored expression for smoothness detection
  // Options: 'rho', 'p', 'rho*p', 'p/rho', 'rho/p', 'velocity_magnitude'
  m_modalMonitoredExpression = "rho*p";
  setParameter( "ModalMonitoredExpression", &m_modalMonitoredExpression );
  
  //==========================================================================
  // BLENDING LIMITS
  //==========================================================================
  
  // Minimum alpha (below this, alpha = 0)
  m_alphaMin = 0.01;
  setParameter( "AlphaMin", &m_alphaMin );
  
  // Maximum alpha cap
  m_alphaMax = 1.0;
  setParameter( "AlphaMax", &m_alphaMax );
  
  //==========================================================================
  // NEIGHBOR AVERAGING & SWEEPS
  //==========================================================================
  
  // Weight for neighbor influence (0 = disabled)
  m_neighborWeight = 0.5;
  setParameter( "NeighborWeight", &m_neighborWeight );
  
  // Number of smoothing sweeps (0 = disabled)
  m_nbSweeps = 0;
  setParameter( "NbSweeps", &m_nbSweeps );
  
  // Sweep parameters
  m_sweepWeight = 0.5;
  setParameter( "SweepWeight", &m_sweepWeight );
  
  m_sweepDamping = 0.7;
  setParameter( "SweepDamping", &m_sweepDamping );
  
  //==========================================================================
  // GENERAL SETTINGS
  //==========================================================================
  
  // Iteration to freeze blending coefficient (very large = never freeze)
  m_freezeFilterIter = 1000000;
  setParameter( "freezeFilterIter", &m_freezeFilterIter );
  
  // Show rate for console output
  m_showrate = 1;
  setParameter( "ShowRate", &m_showrate );
  
  // Smoothness computation method (Modal recommended)
  m_smoothnessMethod = "Modal";
  setParameter( "SmoothnessMethod", &m_smoothnessMethod );
  
  // Monitored variable index (unused when ModalMonitoredExpression is set)
  m_monitoredVar = 0;
  setParameter( "MonitoredVar", &m_monitoredVar);
  
  //==========================================================================
  // DEPRECATED PARAMETERS (Legacy Sigmoid - UseLogScale = false)
  //==========================================================================
  
  // Threshold: T = a * 10^(-c * (N+1)^0.25)
  m_thresholdA = 0.5;
  setParameter( "ThresholdA", &m_thresholdA );
  
  m_thresholdC = 1.8;
  setParameter( "ThresholdC", &m_thresholdC );
  
  // Sigmoid steepness: alpha = 1/(1 + exp(-s/T * (S - T)))
  m_sigmoidS = 9.21024;
  setParameter( "SigmoidS", &m_sigmoidS );

}

//////////////////////////////////////////////////////////////////////////////

BaseOrderBlending::~BaseOrderBlending()
{
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
BaseOrderBlending::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_alpha);
  result.push_back(&socket_prevAlpha);
  result.push_back(&socket_smoothness);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::defineConfigOptions(Config::OptionList& options)
{
  //==========================================================================
  // MAIN PARAMETERS (Log-Scale Method - RECOMMENDED)
  //==========================================================================
  options.addConfigOption< bool >("UseLogScale",
    "Use log-scale smoothness indicator (true, RECOMMENDED) or legacy sigmoid (false). Default: true.");
  options.addConfigOption< CFreal,Config::DynamicOption<> >("S0",
    "Reference smoothness: s0 = -S0 * log10(N+1). Higher = more dissipation. Typical: 3.0-5.0. Default: 4.0.");
  options.addConfigOption< CFreal,Config::DynamicOption<> >("Kappa",
    "Transition half-width in log-space. Blending ramps over [s0-Kappa, s0+Kappa]. Typical: 1.0-2.0. Default: 1.5.");
  options.addConfigOption< std::string >("ModalMonitoredExpression",
    "Expression for smoothness detection: 'rho', 'p', 'rho*p' (default), 'p/rho', 'rho/p', 'velocity_magnitude'.");
  
  //==========================================================================
  // BLENDING LIMITS
  //==========================================================================
  options.addConfigOption< CFreal >("AlphaMin",
    "Minimum blending coefficient. Below this, alpha = 0. Default: 0.01.");
  options.addConfigOption< CFreal >("AlphaMax",
    "Maximum blending coefficient cap. Default: 1.0.");
  
  //==========================================================================
  // NEIGHBOR AVERAGING & SWEEPS
  //==========================================================================
  options.addConfigOption< CFreal >("NeighborWeight",
    "Weight for neighbor influence (0 = disabled). Default: 0.5.");
  options.addConfigOption< CFuint >("NbSweeps",
    "Number of smoothing sweeps (0 = disabled). Default: 0.");
  options.addConfigOption< CFreal >("SweepWeight",
    "Weight factor for neighbor influence during sweeps. Default: 0.5.");
  options.addConfigOption< CFreal >("SweepDamping",
    "Damping factor applied after each sweep. Default: 0.7.");
  
  //==========================================================================
  // GENERAL SETTINGS
  //==========================================================================
  options.addConfigOption< CFuint,Config::DynamicOption<> >("freezeFilterIter",
    "Iteration to freeze blending coefficient. Very large = never. Default: 1000000.");
  options.addConfigOption< CFuint >("ShowRate",
    "Show rate for console output. Default: 1.");
  options.addConfigOption< std::string >("SmoothnessMethod",
    "Method for smoothness: 'Modal' (recommended) or 'Projection'. Default: Modal.");
  options.addConfigOption< CFuint,Config::DynamicOption<> >("MonitoredVar",
    "Index of monitored variable (unused when ModalMonitoredExpression is set). Default: 0.");
  
  //==========================================================================
  // DEPRECATED PARAMETERS (Legacy Sigmoid - only used if UseLogScale = false)
  //==========================================================================
  options.addConfigOption< CFreal >("ThresholdA",
    "[DEPRECATED] Parameter 'a' in threshold T = a * 10^(-c * (N+1)^0.25). Default: 0.5.");
  options.addConfigOption< CFreal >("ThresholdC",
    "[DEPRECATED] Parameter 'c' in threshold T = a * 10^(-c * (N+1)^0.25). Default: 1.8.");
  options.addConfigOption< CFreal >("SigmoidS",
    "[DEPRECATED] Sigmoid steepness in alpha = 1/(1 + exp(-s/T * (S - T))). Default: 9.21024.");
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::execute()
{
  CFTRACEBEGIN;

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  DataHandle< CFreal > output = socket_alpha.getDataHandle();
  DataHandle< CFreal > prevf = socket_prevAlpha.getDataHandle();

  // number of element types, should be 1 in non-mixed grids
  const CFuint nbrElemTypes = elemType->size();
  
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  // variable to store the number of limits done
  m_nbLimits = 0;
  m_totalNbLimits = 0;
  
  // [NOTE] Threshold is now computed inside computeBlendingCoefficient() 
  // based on m_useLogScale flag

  // loop over element types, for the moment there should only be one
  if (iter<m_freezeFilterIter)
  {
  cf_assert(nbrElemTypes == 1);
  for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();


    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      CFreal alpha_neighbor = 0.0;
      CFreal alpha_neighbor_i = 0.0;

      // Loop over the neighbor IDs for the current element
      for (CFuint i = 0; i < m_NeighborIDs[elemIdx].size(); ++i) 
      {
        // Build the GeometricEntity for each neighbor
        geoData.idx = m_NeighborIDs[elemIdx][i];
        m_cell = m_cellBuilder->buildGE();

        // Get the states in this neighbor cell
        m_cellStates = m_cell->getStates();
        
        // Check if neighbor is updatable before computing its influence
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          computeSmoothness();

          alpha_neighbor_i = computeBlendingCoefficient(m_s);
          alpha_neighbor_i = applyAlphaLimits(alpha_neighbor_i);
          alpha_neighbor = std::max(alpha_neighbor, alpha_neighbor_i);
        }

        // Release the GeometricEntity of the neighbor
        m_cellBuilder->releaseGE();
      }

      // build the GeometricEntity for current cell
      geoData.idx = elemIdx;
      m_elemIdx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // Only compute blending for parallel updatable cells
      if ((*m_cellStates)[0]->isParUpdatable())
      {
        // compute the smoothness
        computeSmoothness();

        // Compute Filter Strength
        CFreal alpha = computeBlendingCoefficient(m_s);
        alpha = applyAlphaLimits(alpha);

        m_filterStrength = std::max(alpha, m_neighborWeight * alpha_neighbor);

        // Apply alpha limits after neighbor influence
        m_filterStrength = applyAlphaLimits(m_filterStrength);

        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt) 
        {
          output[((*m_cellStates)[iSolPnt])->getLocalID()] = m_filterStrength;
        }
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }

  /// Sweep to smooth out the blending coefficient strength/value 
  for (CFuint sweep = 0; sweep < m_nbSweeps; ++sweep)
  {
    for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
    {
      // get start and end indexes for this type of element
      const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
      const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();

      // loop over cells
      for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
      {
        CFreal alpha_neighbor = 0.0;
        CFreal alpha_neighbor_i = 0.0;

        // Loop over the neighbor IDs for the current element
        for (CFuint i = 0; i < m_NeighborIDs[elemIdx].size(); ++i) 
        {
          // Build the GeometricEntity for each neighbor
          geoData.idx = m_NeighborIDs[elemIdx][i];
          m_cell = m_cellBuilder->buildGE();

          // Get the states in this neighbor cell
          m_cellStates = m_cell->getStates();

          // Only read alpha from updatable neighbors (avoid ghost cells in parallel)
          if ((*m_cellStates)[0]->isParUpdatable())
          {
            alpha_neighbor_i = output[((*m_cellStates)[0])->getLocalID()];
            alpha_neighbor = std::max(alpha_neighbor, alpha_neighbor_i);
          }

          // Release the GeometricEntity of the neighbor
          m_cellBuilder->releaseGE();
        }

        // build the GeometricEntity for current cell
        geoData.idx = elemIdx;
        m_elemIdx = elemIdx;
        m_cell = m_cellBuilder->buildGE();

        // get the states in this cell
        m_cellStates = m_cell->getStates();

        // Only apply sweeps to parallel updatable cells
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          // get the cell blending coefficient
          m_filterStrength = output[((*m_cellStates)[0])->getLocalID()];

          m_filterStrength = std::max(m_filterStrength, m_sweepWeight * alpha_neighbor);

          // Apply alpha limits after sweeping operations
          CFreal finalAlpha = applyAlphaLimits(m_sweepDamping * m_filterStrength);

          for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt) 
          {
            output[((*m_cellStates)[iSolPnt])->getLocalID()] = finalAlpha;
          }
        }

        //release the GeometricEntity
        m_cellBuilder->releaseGE();
      }
    }
  }

  // After sweeps complete, store final alpha values to prevf for proper freezing
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
        CFreal finalAlpha = output[((*m_cellStates)[0])->getLocalID()];
        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt) 
        {
          prevf[((*m_cellStates)[iSolPnt])->getLocalID()] = finalAlpha;
        }
      }
      m_cellBuilder->releaseGE();
    }
  }
}
  else
  {

    for (m_iElemType = 0; m_iElemType < nbrElemTypes; ++m_iElemType)
    {
      // get start and end indexes for this type of element
      const CFuint startIdx = (*elemType)[m_iElemType].getStartIdx();
      const CFuint endIdx   = (*elemType)[m_iElemType].getEndIdx();
  
      // loop over cells
      for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
      {
        // build the GeometricEntity
        geoData.idx = elemIdx;
        m_elemIdx = elemIdx;
        m_cell = m_cellBuilder->buildGE();

        // get the states in this cell
        m_cellStates = m_cell->getStates();

        // Only update values for parallel updatable cells
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          m_filterStrength = prevf[((*m_cellStates)[0])->getLocalID()];

          for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt) 
          {
            output[((*m_cellStates)[iSolPnt])->getLocalID()] = m_filterStrength;
          }
        }

        //release the GeometricEntity
        m_cellBuilder->releaseGE();
      }
    }
  }


  const std::string nsp = this->getMethodData().getNamespace();
  
#ifdef CF_HAVE_MPI
    MPI_Comm comm = PE::GetPE().GetCommunicator(nsp);
    PE::GetPE().setBarrier(nsp);
    const CFuint count = 1;
    MPI_Allreduce(&m_nbLimits, &m_totalNbLimits, count, MPI_UNSIGNED, MPI_SUM, comm);
    //MPI_Allreduce(&m_nbAvLimits, &m_totalNbAvLimits, count, MPI_UNSIGNED, MPI_SUM, comm);
#endif
    
  if (PE::GetPE().GetRank(nsp) == 0 && iter%m_showrate == 0) 
  {
    // print number of limits
    //CFLog(NOTICE, "Number of times Flux Filtering enforced: " << m_totalNbLimits << "\n");
  }

  PE::GetPE().setBarrier(nsp);
  
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////
CFreal BaseOrderBlending::computeBlendingCoefficient(CFreal smoothness)
{
  if (m_useLogScale)
  {
    // Log-scale approach with sinusoidal ramp (Persson-Peraire style)
    // Smoothness s is already in log10 scale from computeSmoothnessModal()
    // Reference threshold: s0 = -S0 * log10(N+1)
    const CFreal s0 = computeReferenceThreshold();
    
    if (smoothness < s0 - m_kappa)
    {
      // Smooth region: no blending needed
      return 0.0;
    }
    else if (smoothness > s0 + m_kappa)
    {
      // Rough region: full blending
      return 1.0;
    }
    else
    {
      // Transition region: smooth sinusoidal ramp
      // alpha = 0.5 * (1 + sin(π*(s-s0)/(2*κ)))
      return 0.5 * (1.0 + std::sin(MathTools::MathConsts::CFrealPi() * (smoothness - s0) / (2.0 * m_kappa)));
    }
  }
  else
  {
    // [DEPRECATED] Legacy sigmoid approach
    const CFreal threshold = computeThreshold_Legacy(m_order, m_thresholdA, m_thresholdC);
    return 1.0 / (1.0 + std::exp(-m_sigmoidS / threshold * (smoothness - threshold)));
  }
}

//////////////////////////////////////////////////////////////////////////////
CFreal BaseOrderBlending::computeReferenceThreshold() const
{
  // Reference threshold for log-scale method: s0 = -S0 * log10(N+1)
  // Higher S0 means more dissipation (stricter smoothness requirement)
  return -m_s0 * std::log10(static_cast<CFreal>(m_order + 1));
}

//////////////////////////////////////////////////////////////////////////////
CFreal BaseOrderBlending::computeThreshold_Legacy(CFuint N, CFreal a, CFreal c)
{
  // [DEPRECATED] Legacy threshold function for sigmoid blending
  // T = a * 10^(-c * (N+1)^0.25)
  return a * std::pow(10.0, -c * std::pow(N + 1, 0.25));
}

//////////////////////////////////////////////////////////////////////////////
CFreal BaseOrderBlending::applyAlphaLimits(CFreal alpha)
{
  if (alpha < m_alphaMin)
  {
    alpha= 0.0;
  }
  else if (alpha > 1.0 - m_alphaMin)
  {
    alpha= 1.0;
  }

  return std::min(alpha, m_alphaMax);
}

//////////////////////////////////////////////////////////////////////////////
void BaseOrderBlending::computeSmoothness()
{
  if (m_smoothnessMethod == "Projection")
  {
    computeSmoothnessProjection();
  }
  else // Default to Modal method
  {
    computeSmoothnessModal();
  }
}

//////////////////////////////////////////////////////////////////////////////
void BaseOrderBlending::computeSmoothnessProjection()
{ 
  CFreal sNum = 0.0;
  CFreal sDenom = 0.0;
  
  // Compute projected states to P-1
  computeProjStates(m_statesPMinOne);
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    CFreal stateP = (*((*m_cellStates)[iSol]))[m_monitoredVar];
    CFreal diffStatesPPMinOne = stateP - m_statesPMinOne[iSol][m_monitoredVar];
    sNum += diffStatesPPMinOne*diffStatesPPMinOne;
    sDenom += stateP*stateP;
  }
  
  m_s = sNum / std::max(sDenom, MathTools::MathConsts::CFrealEps());
  
  // Store the smoothness value to all solution points of this element
  DataHandle< CFreal > smoothness = socket_smoothness.getDataHandle();
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    smoothness[(((*m_cellStates)[iSol]))->getLocalID()] = m_s;
  }
  
  if (m_s > m_Smax)
  {
      m_Smax = m_s;
  }
}

//////////////////////////////////////////////////////////////////////////////
void BaseOrderBlending::computeSmoothnessModal()
{
  std::vector<CFreal> modalCoeffs(m_nbrSolPnts, 0.0);
  CFreal max_indicator = 0.0;

  // Step 1: Load sol pnts values for monitored variable based on expression (optimized)
  if (m_modalMonitoredExpression == "rho")
  {
    // Extract density values
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_tempSolPntVec[iSol] = (*((*m_cellStates)[iSol]))[0];
    }
  }
  else if (m_modalMonitoredExpression == "p")
  {
    // Extract pressure values
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_tempSolPntVec[iSol] = (*((*m_cellStates)[iSol]))[7];
    }
  }
  else
  {
    // For more complex expressions, we need the loop but it's more efficient
    // to check the expression once outside the loop
    const bool isRhoP = (m_modalMonitoredExpression == "rho*p");
    const bool isPRho = (m_modalMonitoredExpression == "p/rho");
    const bool isRhoOverP = (m_modalMonitoredExpression == "rho/p");
    const bool isVelMag = (m_modalMonitoredExpression == "velocity_magnitude");
    
    if (!isRhoP && !isPRho && !isRhoOverP && !isVelMag)
    {
      CFLog(WARN, "Unknown ModalMonitoredExpression '" << m_modalMonitoredExpression 
                  << "', defaulting to 'rho*p'\n");
    }
    
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      if (isRhoP || (!isPRho && !isRhoOverP && !isVelMag)) // Default case
      {
        m_tempSolPntVec[iSol] = (*((*m_cellStates)[iSol]))[0] * (*((*m_cellStates)[iSol]))[7];
      }
      else if (isPRho)
      {
        CFreal rho = (*((*m_cellStates)[iSol]))[0];
        m_tempSolPntVec[iSol] = (*((*m_cellStates)[iSol]))[7] / std::max(rho, MathTools::MathConsts::CFrealEps());
      }
      else if (isRhoOverP)
      {
        CFreal p = (*((*m_cellStates)[iSol]))[7];
        m_tempSolPntVec[iSol] = (*((*m_cellStates)[iSol]))[0] / std::max(p, MathTools::MathConsts::CFrealEps());
      }
      else if (isVelMag)
      {
        CFreal vx = (*((*m_cellStates)[iSol]))[1];
        CFreal vy = (*((*m_cellStates)[iSol]))[2];
        CFreal vz = (*((*m_cellStates)[iSol]))[3];
        m_tempSolPntVec[iSol] = std::sqrt(vx*vx + vy*vy + vz*vz);
      }
    }
  }

  // Step 2: Transform to modal space: modalCoeffs = V^{-1} * solPntsValues
  m_tempSolPntVec2 = m_vdmInv * m_tempSolPntVec;

  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    modalCoeffs[iSol] = m_tempSolPntVec2[iSol];
  }

  // Step 3: Compute modal energy ratios
  // E1 = energy in P modes (highest) / total energy
  // E2 = energy in P-1 modes / energy in modes 0 to P-1 (to catch odd-even decoupling)
  // 
  // For P<=2, E2 is problematic because P-1 modes represent physical linear gradients
  // So we only use E2 for P>=3
  
  const CFreal eps = MathTools::MathConsts::CFrealEps();
  
  CFreal energy_total = 0.0;   // Total energy (all modes)
  CFreal energy_P = 0.0;       // Energy in P modes (highest order)
  CFreal energy_Pm1 = 0.0;     // Energy in P-1 modes (second highest)
  CFreal energy_low = 0.0;     // Energy in modes 0 to P-2 (low modes)

  for (CFuint j = 0; j < m_nbrSolPnts; ++j)
  {
    const CFreal mj2 = modalCoeffs[j] * modalCoeffs[j];
    energy_total += mj2;
    
    if (j >= m_nbrSolPntsMinOne)
    {
      // P modes (highest)
      energy_P += mj2;
    }
    else if (j >= m_nbrSolPntsMinTwo)
    {
      // P-1 modes
      energy_Pm1 += mj2;
    }
    else
    {
      // Modes 0 to P-2 (low modes)
      energy_low += mj2;
    }
  }

  // E1: Highest mode energy ratio (always computed)
  const CFreal E1 = energy_P / std::max(energy_total, eps);
  
  // Compute final energy ratio based on polynomial order
  CFreal E_var;
  if (m_order >= 3)
  {
    // For P>=3: Use both E1 and E2 to catch odd-even decoupling
    // E2 = energy in P-1 modes / energy in modes 0 to P-1
    const CFreal E2 = energy_Pm1 / std::max(energy_low + energy_Pm1, eps);
    E_var = std::max(E1, E2);
  }
  else
  {
    // For P<=2: Only use E1 (E2 would flag physical linear gradients)
    E_var = E1;
  }

  // Apply log transform if using log-scale method
  if (m_useLogScale)
  {
    // Log-scale smoothness indicator: s = log10(E_var)
    // Smooth regions: E_var ~ 10^-10 to 10^-6 => s ~ -10 to -6
    // Rough regions:  E_var ~ 10^-2 to 10^0  => s ~ -2 to 0
    m_s = std::log10(std::max(E_var, eps));
  }
  else
  {
    // [DEPRECATED] Legacy linear smoothness indicator
    m_s = E_var;
  }

  // Store the smoothness value to all solution points of this element
  DataHandle< CFreal > smoothness = socket_smoothness.getDataHandle();
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    smoothness[(((*m_cellStates)[iSol]))->getLocalID()] = m_s;
  }

  if (m_s > m_Smax)
  {
    m_Smax = m_s;
  }
}

//////////////////////////////////////////////////////////////////////////////
void BaseOrderBlending::computeProjStates(std::vector< RealVector >& projStates)
{
  cf_assert(m_nbrSolPnts == projStates.size());
  
  if (m_order != 1)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        m_tempSolPntVec[iSol] = (*((*m_cellStates)[iSol]))[iEq];
      }

      m_tempSolPntVec2 = m_transformationMatrix*m_tempSolPntVec;

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        projStates[iSol][iEq] = m_tempSolPntVec2[iSol];
      }
    }
  }
  else
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      CFreal stateSum = 0.0;
      
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        stateSum += (*((*m_cellStates)[iSol]))[iEq];
      }

      stateSum /= m_nbrSolPnts;

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        projStates[iSol][iEq] = stateSum;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////
RealVector BaseOrderBlending::getmaxModalOrder(const CFGeoShape::Type elemShape, const CFuint m_order) 
{
  RealVector maxModalOrder; // This will be filled based on element shape and order

  // Compute the correct function based on element shape
  switch(elemShape) {
      case CFGeoShape::QUAD:
      {
        // Calculate the number of modes for quadrilateral elements
        CFuint totalModes = (m_order + 1) * (m_order + 1);
        maxModalOrder.resize(totalModes);

        CFuint modeIndex = 0;
        for (CFuint p = 0; p <= m_order; ++p) {
            for (CFuint i = (p)*(p); i < (p+1)*(p+1); ++i) {
                maxModalOrder[modeIndex] = p;
                modeIndex+=1;
            }
        }
        cf_assert(modeIndex == totalModes);
        break;
      }
      case CFGeoShape::TRIAG:
      {
        CFuint totalModes = ((m_order + 1)*(m_order + 2))/2.;
        maxModalOrder.resize(totalModes);

        CFuint modeIndex = 0;
        for (CFuint totalOrder = 0;  totalOrder <= m_order; ++totalOrder)
        {
        // Within each total order, iterate over powers of ksi and eta
        for (CFuint iOrderKsi = totalOrder; iOrderKsi >= 0; --iOrderKsi)
          {
            CFuint iOrderEta = totalOrder - iOrderKsi;
            maxModalOrder[modeIndex] = std::max(iOrderKsi, iOrderEta);
            modeIndex+=1;
            if (iOrderKsi == 0) break;
          }
        }
        cf_assert(modeIndex == totalModes);
        break;
      }
      case CFGeoShape::PRISM:
      {
        CFuint totalModes = ((m_order + 1)*(m_order + 1)*(m_order + 2))/2.;
        maxModalOrder.resize(totalModes);

        CFuint modeIndex = 0;
        // Loop over the total order of polynomial terms for ξ and η
        for (CFuint totalOrderXY = 0; totalOrderXY <= m_order; ++totalOrderXY)
        {
          // Loop for the third dimension (Zta)
          for (CFuint iOrderZta = 0; iOrderZta <= m_order; ++iOrderZta)
          {
            // Within each total order, iterate over powers of ksi and eta in reverse for ksi
            for (CFuint iOrderKsi = totalOrderXY; iOrderKsi >= 0; --iOrderKsi)
            {
              CFuint iOrderEta = totalOrderXY - iOrderKsi;
              maxModalOrder[modeIndex] = totalOrderXY;
              modeIndex+=1;
              if (iOrderKsi == 0) break;
            }
          }
        }
        cf_assert(modeIndex == totalModes);
        break;
      }      
      case CFGeoShape::HEXA:
      {
        // Calculate the number of modes for hexahedral elements
        CFuint totalModes = (m_order + 1) * (m_order + 1) * (m_order + 1);
        maxModalOrder.resize(totalModes);

        CFuint modeIndex = 0;
        for (CFuint iOrderKsi = 0; iOrderKsi <= m_order; ++iOrderKsi) {
            for (CFuint iOrderEta = 0; iOrderEta <= m_order; ++iOrderEta) {
                for (CFuint iOrderZta = 0; iOrderZta <= m_order; ++iOrderZta) {
                    CFuint maxOrder = std::max({iOrderKsi, iOrderEta, iOrderZta});
                    maxModalOrder[modeIndex] = maxOrder;
                    modeIndex += 1;
                }
            }
        }
        cf_assert(modeIndex == totalModes);
        break;
      }
      default:
      {
          throw Common::NotImplementedException(FromHere(), "Maximal modal order calculation not implemented for elements of type " + StringOps::to_str(elemShape) + ".");
      }
  }

  return maxModalOrder;
}

//////////////////////////////////////////////////////////////////////////////

void BaseOrderBlending::setup()
{
  CFAUTOTRACE;

  // get number of equations and dimensions
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim    = PhysicalModelStack::getActive()->getDim ();

  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  const CFuint nbrElemTypes = frLocalData.size();
  cf_assert(nbrElemTypes > 0);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  const CFGeoShape::Type elemShape = frLocalData[0]->getShape();
  
  m_order = static_cast<CFuint>(order);
  
  // get nbr of sol pnts
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();

  // get the maximum number of flux points
  m_maxNbrFlxPnts = 0;
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  
  // this is only necessary for mixed grids and is a bit superfluous for now
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbrFlxPnts = frLocalData[iElemType]->getNbrOfFlxPnts();
    m_maxNbrFlxPnts = m_maxNbrFlxPnts > nbrFlxPnts ? m_maxNbrFlxPnts : nbrFlxPnts;
  }
  
  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();
  
  // initialize extrapolated states
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  {
    RealVector temp(m_nbrEqs);
    m_cellStatesFlxPnt.push_back(temp);
  }
  
  // get the number of cells in the mesh
  const CFuint nbrCells = (*elemType)[0].getEndIdx();
  
  // get datahandle
  DataHandle< CFreal > prevF = socket_prevAlpha.getDataHandle();
  DataHandle< CFreal > output = socket_alpha.getDataHandle();
  DataHandle< CFreal > smoothness = socket_smoothness.getDataHandle();
  
  const CFuint nbStates = nbrCells*m_nbrSolPnts;

  // resize socket
  output.resize(nbStates);
  prevF.resize(nbStates);
  smoothness.resize(nbStates);

  // initialize extrapolated states
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    RealVector temp(m_nbrEqs);
    m_tempStates.push_back(temp);
  }
  m_tempSolPntVec.resize(m_nbrSolPnts);
  m_tempSolPntVec2.resize(m_nbrSolPnts);

  m_vdm = *(frLocalData[0]->getVandermondeMatrix());
  
  m_vdmInv = *(frLocalData[0]->getVandermondeMatrixInv());
  
  m_NeighborIDs.resize(nbrCells);

  m_maxModalOrder.resize(m_nbrSolPnts);
  m_maxModalOrder = getmaxModalOrder(elemShape, m_order);

  m_Smax = m_s0 + m_kappa;


  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    RealVector temp(m_nbrEqs);
    temp = 0.0;
    m_statesPMinOne.push_back(temp);
  }

  RealMatrix temp(m_nbrSolPnts,m_nbrSolPnts);
  temp = 0.0;
  if (m_dim == 2)
  {
    if (elemShape == CFGeoShape::TRIAG){  
      for (CFuint idx = 0; idx < (m_order)*(m_order+1)/2; ++idx)
      {
        temp(idx,idx) = 1.0;
      }
    }
    else{
      for (CFuint idx = 0; idx < (m_order)*(m_order); ++idx)
      {
        temp(idx,idx) = 1.0;
      }
    }  
  }
  else if (m_dim == 3)
  {
    if (elemShape == CFGeoShape::TETRA){ 
      for (CFuint idx = 0; idx < (m_order)*(m_order+1)*(m_order+2)/6; ++idx)
      {
        temp(idx,idx) = 1.0;
      }
    }
    else if (elemShape == CFGeoShape::PRISM){ 
      for (CFuint idx = 0; idx < (m_order)*(m_order)*(m_order+1)/2; ++idx)
      {
        temp(idx,idx) = 1.0;
      }
    }
    else{
        for (CFuint idx = 0; idx < (m_order)*(m_order)*(m_order); ++idx)
      {
        temp(idx,idx) = 1.0;
      }
    }  
  }


  m_transformationMatrix.resize(m_nbrSolPnts,m_nbrSolPnts);

  m_transformationMatrix = (m_vdm)*temp*(m_vdmInv);


  if (m_dim == 2)
  {
    if (elemShape == CFGeoShape::TRIAG){  
      m_nbrSolPntsMinOne = (m_order)*(m_order+1)/2; 
      m_nbrSolPntsMinTwo = (m_order > 1) ? ((m_order-1)*(m_order)/2) : 1 ; 
    }
    else{ //QUADS
      m_nbrSolPntsMinOne = m_order*m_order; 
      m_nbrSolPntsMinTwo = ((m_order-1)*(m_order-1)); 
    }  
  }
  else if (m_dim == 3)
  {
    if (elemShape == CFGeoShape::TETRA){ 
      m_nbrSolPntsMinOne = (m_order)*(m_order+1)*(m_order+2)/6; 
      m_nbrSolPntsMinTwo = (m_order > 1) ? ((m_order-1)*(m_order)*(m_order+1)/6) : 1 ; 
    }
    else if (elemShape == CFGeoShape::PRISM){ 
      for (CFuint idx = 0; idx < (m_order)*(m_order)*(m_order+1)/2; ++idx)
      {
        temp(idx,idx) = 1.0;
      }
      m_nbrSolPntsMinOne = (m_order)*(m_order)*(m_order+1)/2; 
      m_nbrSolPntsMinTwo = (m_order > 1) ? ((m_order-1)*(m_order-1)*(m_order)/2) : 1 ; 
    }
    else{
        for (CFuint idx = 0; idx < (m_order)*(m_order)*(m_order); ++idx)
      {
        temp(idx,idx) = 1.0;
      }
      m_nbrSolPntsMinOne = (m_order)*(m_order)*(m_order); 
      m_nbrSolPntsMinTwo = (m_order > 1) ? ((m_order-1)*(m_order-1)*(m_order-1)) : 1 ; 
    }  
  }

  ////////
  m_geoBuilder = getMethodData().getCellBuilder();
  
  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  CellToFaceGEBuilder::GeoData& geoData = m_geoBuilder->getDataGE();
  geoData.trs = cells;

  // Loop on elements to find neighbors via shared nodes
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType) 
  {
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx();
  
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      std::set<CFuint> currNodeIDs;
  
      // Step 1: Gather node IDs for the current element
      const CFuint nbNodesCurr = cells->getNbNodesInGeo(elemIdx);
      for (CFuint iNode = 0; iNode < nbNodesCurr; ++iNode)
      {
        CFuint nodeID = cells->getNodeID(elemIdx, iNode);
        currNodeIDs.insert(nodeID);
      }
  
      std::set<CFuint> neighborSet; // avoid duplicates
  
      // Step 2: Loop over all other elements in this element type
      for (CFuint nIdx = startIdx; nIdx < endIdx; ++nIdx)
      {
        if (nIdx == elemIdx) continue;
  
        const CFuint nbNodesNeighbor = cells->getNbNodesInGeo(nIdx);
        for (CFuint jNode = 0; jNode < nbNodesNeighbor; ++jNode)
        {
          CFuint neighborNodeID = cells->getNodeID(nIdx, jNode);
  
          if (currNodeIDs.find(neighborNodeID) != currNodeIDs.end())
          {
            neighborSet.insert(nIdx); // shared node found
            break; // one common node is enough
          }
        }
      }
  
      // Step 3: Store unique neighbor IDs
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