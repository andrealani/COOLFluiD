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
  
  m_showrate = 1;
  setParameter( "ShowRate", &m_showrate );

  m_freezeFilterIter = 1000000;
  setParameter( "freezeFilterIter", &m_freezeFilterIter );

  m_kappa = 1.0;
  setParameter( "Kappa", &m_kappa);

  m_s0 = 0.0;
  setParameter( "S0", &m_s0);

  m_monitoredVar = 0;
  setParameter( "MonitoredVar", &m_monitoredVar);

  m_nbSweeps = 1;
  setParameter( "NbSweeps", &m_nbSweeps );

  m_thresholdA = 0.5;
  setParameter( "ThresholdA", &m_thresholdA );

  m_thresholdC = 1.8;
  setParameter( "ThresholdC", &m_thresholdC );

  m_sigmoidS = 9.21024;
  setParameter( "SigmoidS", &m_sigmoidS );

  m_alphaMin = 0.01;
  setParameter( "AlphaMin", &m_alphaMin );

  m_alphaMax = 1.0;
  setParameter( "AlphaMax", &m_alphaMax );

  m_neighborWeight = 0.7;
  setParameter( "NeighborWeight", &m_neighborWeight );

  m_sweepWeight = 0.5;
  setParameter( "SweepWeight", &m_sweepWeight );

  m_sweepDamping = 0.7;
  setParameter( "SweepDamping", &m_sweepDamping );

  m_smoothnessMethod = "Modal";
  setParameter( "SmoothnessMethod", &m_smoothnessMethod );

  m_modalMonitoredExpression = "rho*p";
  setParameter( "ModalMonitoredExpression", &m_modalMonitoredExpression );

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
  options.addConfigOption< CFuint >("ShowRate","Show rate of the filtering information.");
  options.addConfigOption< CFuint,Config::DynamicOption<> >("freezeFilterIter","Iteration at which the filter strength will be frozen.");
  options.addConfigOption< CFreal,Config::DynamicOption<> >("Kappa","Kappa factor of artificial viscosity.");
  options.addConfigOption< CFreal,Config::DynamicOption<> >("S0","Reference smoothness factor, will be multiplied by -log(P).");
  options.addConfigOption< CFuint,Config::DynamicOption<> >("MonitoredVar","Index of the monitored var for positivity preservation.");
  
  // Order blending specific parameters
  options.addConfigOption< CFuint >("NbSweeps","Number of smoothing sweeps for blending coefficients.");
  options.addConfigOption< CFreal >("ThresholdA","Parameter 'a' in threshold function T = a * 10^(-c * (N+1)^0.25).");
  options.addConfigOption< CFreal >("ThresholdC","Parameter 'c' in threshold function T = a * 10^(-c * (N+1)^0.25).");
  options.addConfigOption< CFreal >("SigmoidS","Parameter 's' in sigmoid function alpha = 1/(1 + exp(-s/T * (S - T))).");
  options.addConfigOption< CFreal >("AlphaMin","Minimum blending coefficient (below this value, alpha = 0).");
  options.addConfigOption< CFreal >("AlphaMax","Maximum blending coefficient (above 1-AlphaMin, alpha = 1).");
  options.addConfigOption< CFreal >("NeighborWeight","Weight factor for neighbor influence in blending computation.");
  options.addConfigOption< CFreal >("SweepWeight","Weight factor for neighbor influence during sweeps.");
  options.addConfigOption< CFreal >("SweepDamping","Damping factor applied after each sweep iteration.");
  options.addConfigOption< std::string >("SmoothnessMethod","Method for computing smoothness: 'Modal' or 'Projection'.");
  options.addConfigOption< std::string >("ModalMonitoredExpression","Expression for monitored variable in modal method: 'rho', 'p', 'rho*p', 'p/rho', 'rho/p', or 'velocity_magnitude'.");
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
  
  // Compute threshold using configurable parameters
  CFreal T = computeThreshold(m_order, m_thresholdA, m_thresholdC);

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

          alpha_neighbor_i = computeBlendingCoefficient(m_s, T);
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
        CFreal alpha = computeBlendingCoefficient(m_s, T);
        alpha = applyAlphaLimits(alpha);

        m_filterStrength = std::max(alpha, m_neighborWeight * alpha_neighbor);

        // Apply alpha limits after neighbor influence
        m_filterStrength = applyAlphaLimits(m_filterStrength);

        for (CFuint iSolPnt = 0; iSolPnt < m_nbrSolPnts; ++iSolPnt) 
        {
          output[((*m_cellStates)[iSolPnt])->getLocalID()] = m_filterStrength;
          prevf[((*m_cellStates)[iSolPnt])->getLocalID()] = m_filterStrength;
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

          alpha_neighbor_i = output[((*m_cellStates)[0])->getLocalID()];
          alpha_neighbor = std::max(alpha_neighbor, alpha_neighbor_i);

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
CFreal BaseOrderBlending::computeBlendingCoefficient(CFreal smoothness, CFreal threshold)
{
  return 1.0 / (1.0 + std::exp(-m_sigmoidS / threshold * (smoothness - threshold)));
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
CFreal BaseOrderBlending::computeThreshold(CFuint N, CFreal a, CFreal c)
{
  return a * std::pow(10.0, -c * std::pow(N + 1, 0.25));
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

  // Step 3: Compute modal energy and high-mode content
  CFreal energy_total = 0.0;
  CFreal energy_pMin1 = 0.0;
  CFreal mN2 = 0.0;
  CFreal mNminOne2 = 0.0;

  for (CFuint j = 0; j < m_nbrSolPnts; ++j)
  {
    CFreal mj2 = modalCoeffs[j] * modalCoeffs[j];
    energy_total += mj2;
    if (j > m_nbrSolPntsMinOne - 1) mN2 += mj2;
  }

  for (CFuint j = m_nbrSolPntsMinTwo ; j < m_nbrSolPntsMinOne; ++j)
  {
    CFreal mj2 = modalCoeffs[j] * modalCoeffs[j];
    mNminOne2 += mj2;
  }

  for (CFuint j = 0; j < m_nbrSolPntsMinOne; ++j)
  {
    CFreal mj2 = modalCoeffs[j] * modalCoeffs[j];
    energy_pMin1 += mj2;
  }

  CFreal E1 = mN2 / std::max(energy_total, MathTools::MathConsts::CFrealEps());
  CFreal E2 = mNminOne2 / std::max(energy_pMin1, MathTools::MathConsts::CFrealEps());
  
  if (m_order==2) {E2 = 0.5*E2;}
 
  CFreal E_var = std::max(E1, E2);

  if (m_order==1) {E_var = E1;}

  max_indicator = std::max(max_indicator, E_var);

  m_s = max_indicator;

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