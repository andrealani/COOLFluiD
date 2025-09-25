#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/MLPLimiter.hh"
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

MethodCommandProvider<MLPLimiter, FluxReconstructionSolverData, FluxReconstructionModule> MLPLimiterFluxReconstructionProvider("MLPLimiter");

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal,Config::DynamicOption<> >("LimFactor","K factor in the MLP-u2 limiting part.");
  
  options.addConfigOption< CFreal,Config::DynamicOption<> >("FreezeLimiterRes","Residual after which to freeze the residual.");
  
  options.addConfigOption< CFuint,Config::DynamicOption<> >("FreezeLimiterIter","Iteration after which to freeze the residual.");

  options.addConfigOption< CFuint >("ShowRate","Showrate of the limiter information.");
}

//////////////////////////////////////////////////////////////////////////////

MLPLimiter::MLPLimiter(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_nodeNghbCellMinAvgStates("nodeNghbCellMinAvgStates"),
  socket_nodeNghbCellMaxAvgStates("nodeNghbCellMaxAvgStates"),
  socket_outputLimiting("outputLimiting"),
  m_limiterValues(),
  m_cellBuilder(CFNULL),
  m_cell(),
  m_cellStates(),
  m_cellStatesBackup(),
  m_cellNodes(),
  m_cellAvgSolCoefs(),
  m_cellCenterDerivCoefs(),
  m_flxPntsRecCoefs(),
  m_allFlxPntIdxs(),
  m_solPntsLocalCoords(),
  m_cellAvgState(),
  m_cellCenterDerivVar(),
  m_minAvgState(),
  m_maxAvgState(),
  m_minAvgStateAll(),
  m_maxAvgStateAll(),
  m_nbrEqs(),
  m_dim(),
  m_nbrFlxPnts(),
  m_nbrSolPnts(),
  m_nbrCornerNodes(),
  m_applyLimiter(),
  m_applyLimiterToPhysVar(),
  m_tvbLimitFactor(),
  m_freezeLimiterRes(),
  m_freezeLimiterIter(),
  m_lengthScaleExp(),
  m_cellStatesNodes(),
  m_cellStatesNodesP1(),
  m_maxNbrFlxPnts(),
  m_solPolyValsAtNodes(CFNULL),
  m_transformationMatrices(),
  m_statesP1(),
  m_states2(),
  m_nbrNodesElem(),
  m_solPolyDerivAtNodes(CFNULL),
  m_order(),
  m_solPolyValsAtFlxPnts(CFNULL),
  m_cellStatesFlxPnt(),
  m_prevStates()
{
  addConfigOptionsTo(this);

  m_tvbLimitFactor = 0.0;
  setParameter( "LimFactor", &m_tvbLimitFactor);
  
  m_freezeLimiterRes = -20.0;
  setParameter( "FreezeLimiterRes", &m_freezeLimiterRes);
  
  m_freezeLimiterIter = MathTools::MathConsts::CFuintMax();
  setParameter( "FreezeLimiterIter", &m_freezeLimiterIter);

  m_showrate= 1;
  setParameter( "ShowRate", &m_showrate );
}

//////////////////////////////////////////////////////////////////////////////

MLPLimiter::~MLPLimiter()
{
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::execute()
{
  CFTRACEBEGIN;

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;

  // RESET THE MINIMUM AND MAXIMUM CELL AVERAGED SOLUTIONS (STORED IN THE NODES!!!)
  resetNodeNghbrCellAvgStates();
  
  // get datahandle for limiting output and initialize to 0.0 (no limiting)
  DataHandle< CFreal > outputLimiting = socket_outputLimiting.getDataHandle();
  outputLimiting = 0.0;

  const CFreal residual = SubSystemStatusStack::getActive()->getResidual();
  
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  
  const bool recomputeLimiter = true; //residual > m_freezeLimiterRes && iter < m_freezeLimiterIter;
  
  const bool useMin = residual < m_freezeLimiterRes || iter > m_freezeLimiterIter;
  
  RealVector duL2Norm(m_nbrEqs);
  duL2Norm = 0.0;
  
  const std::string nsp = this->getMethodData().getNamespace();

  const CFuint nbrElemTypes = elemType->size();
  
  if (recomputeLimiter)
  {
  // LOOP OVER ELEMENT TYPES, TO SET THE MINIMUM AND MAXIMUM CELL AVERAGED STATES IN THE NODES
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx  ();

    // set reconstruction data for averaged solutions
    setAvgReconstructionData(iElemType);

    // loop over cells, to put the minimum and maximum neighbouring cell averaged states in the nodes
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // get the nodes in this cell
      m_cellNodes  = m_cell->getNodes();
      
      // reconstruct cell averaged state
      reconstructCellAveragedState();

      // set the min and max node neighbour cell averaged states
      setMinMaxNodeNghbrCellAvgStates();

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }

  CFuint limitCounter = 0;
  CFuint orderLimit = 0;

  // LOOP OVER ELEMENT TYPES, TO LIMIT THE SOLUTIONS
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx  ();

    // set reconstruction data
    setAllReconstructionData(iElemType);

    // loop over cells
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // get the nodes in this cell
      m_cellNodes  = m_cell->getNodes();

      // set the min and max neighbouring cell averaged states
      setMinMaxNghbrCellAvgStates();
      
      // compute P1 projection of states
      computeProjStates(m_statesP1,1);

      // compute the P1 states in the nodes
      computeNodeStates(m_statesP1,m_cellStatesNodesP1);

      // reconstruct cell averaged state
      reconstructCellAveragedState();

      if(m_cell->getID() == 36)
      {
	CFLog(VERBOSE, "states: \n");
	for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
	{
	  CFLog(VERBOSE, (*((*m_cellStates)[iSol])) << "\n");
	}
	CFLog(VERBOSE, "av state: " << m_cellAvgState << "\n");
	CFLog(VERBOSE, "P1 states: \n");
	for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
	{
	  CFLog(VERBOSE, m_statesP1[iSol] << "\n");
	}
	CFLog(VERBOSE, "P1 node states: \n");
	for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
        {
	  CFLog(VERBOSE, m_cellStatesNodesP1[iNode][0] << "\n");
	}
	CFLog(VERBOSE, "min state: " << m_minAvgState[0] << ", max state: " << m_maxAvgState[0] << "\n");
      }

      if (checkSpecialLimConditions())
      {

      // Initiate hierarchical procedure at order P and break from loop if no limiting is necessary
      for (CFuint iOrder = m_order; iOrder > 0; --iOrder)
      {
	// copy the states in the correct format
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          m_states2[iSol] = (*((*m_cellStates)[iSol]));
        }

        // compute the states in the nodes
        computeNodeStates(m_states2,m_cellStatesNodes);
	
	if(m_cell->getID() == 36)
        {
	  CFLog(VERBOSE, "P2NodeStates: \n");
	  for (CFuint iSol = 0; iSol < m_nbrNodesElem; ++iSol)
	  {
	    CFLog(VERBOSE, m_cellStatesNodes[iSol] << "\n");
	  }
	}

	bool limitFlag = false;
	
        // check if P1u is within vertex limiting averages
        for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
        {
	  m_applyLimiter[iNode] = false;

	  // check if we are in constant area of sol
	  if (!( fabs(m_cellStatesNodes[iNode][0] - m_cellAvgState[0]) < max(1e-3 * fabs(m_cellAvgState[0]),m_cell->computeVolume())))
	  {   
            m_applyLimiter[iNode] = m_cellStatesNodesP1[iNode][0] < 1.0*m_minAvgState[0] || m_cellStatesNodesP1[iNode][0] > 1.0*m_maxAvgState[0];
	    
	    if(m_cell->getID() == 36)
            {
	      bool tempBool = m_cellStatesNodesP1[iNode][0] <  m_minAvgState[0];
	      CFLog(VERBOSE, "min? " << tempBool << "\n");
	      tempBool = m_cellStatesNodesP1[iNode][0] >  m_maxAvgState[0];
	      CFLog(VERBOSE, "max? " << tempBool << "\n");
	      CFLog(VERBOSE, "prelim? " << m_applyLimiter[iNode] << "\n");
	    }
	  }

	  if (m_applyLimiter[iNode])
	  {
	    limitFlag = true;
	    CFLog(VERBOSE, "prelimiting cell " << m_cell->getID() << "\n");
	  }
        }
      
        //const bool unphysical = !checkPhysicality();
        const bool unphysical = false; 
        
        if (limitFlag || unphysical)
        {
	
	  if(m_cell->getID() == 36)
          {
	    CFLog(VERBOSE, "node states: \n");
	    for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
            {
	      CFLog(VERBOSE, m_cellStatesNodes[iNode][0] << "\n");
	    }
	  }

	  bool stillLimit = false;
	
	  // get the data handles for the minimum and maximum nodal states
          DataHandle< RealVector > nodeNghbCellMinAvgStates = socket_nodeNghbCellMinAvgStates.getDataHandle();
          DataHandle< RealVector > nodeNghbCellMaxAvgStates = socket_nodeNghbCellMaxAvgStates.getDataHandle();
	
	  if (iOrder != 1)
	  {
	    // check for smooth extrema
            for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
            {
	      if (m_applyLimiter[iNode])
	      {
	        // get node ID
                const CFuint nodeID = (*m_cellNodes)[iNode]->getLocalID();

	        // get min and max states
                RealVector minAvgState = nodeNghbCellMinAvgStates[nodeID];
                RealVector maxAvgState = nodeNghbCellMaxAvgStates[nodeID];
		
		// Check if the vertex is on the boundary, if so ignore the vertex
		///@todo this won't work for simplex elements
		if(m_applyLimiter[iNode] && maxAvgState[m_nbrEqs] < m_nbrCornerNodes)
                {
		  //m_applyLimiter[iNode] = false;
		}
	    
	        if (m_applyLimiter[iNode])
                {
                  // check if we are in a smooth local maximum
	          bool cond1 = m_cellStatesNodesP1[iNode][0] - 1.0*m_cellAvgState[0] > 0.0;
	          bool cond2 = 1.0*m_cellStatesNodes[iNode][0] - m_cellStatesNodesP1[iNode][0] < 0.0;
	          bool cond3 = m_cellStatesNodes[iNode][0] > 1.0*minAvgState[0];
                  m_applyLimiter[iNode] = !( cond1 && cond2 && cond3 );	
	          CFLog(VERBOSE, "max cond: " << cond1 << ", " << cond2 << ", " << cond3 << "\n");
                }
       
                if(m_applyLimiter[iNode])
                {
	          // check if we are in a smooth local minimum
	          bool cond1 = m_cellStatesNodesP1[iNode][0] - 1.0*m_cellAvgState[0] < 0.0;
	          bool cond2 = 1.0*m_cellStatesNodes[iNode][0] - m_cellStatesNodesP1[iNode][0] > 0.0;
	          bool cond3 = m_cellStatesNodes[iNode][0] < 1.0*maxAvgState[0];
                  m_applyLimiter[iNode] = !( cond1 && cond2 && cond3 );	
	          CFLog(VERBOSE, "min cond: " << cond1 << ", " << cond2 << ", " << cond3 << "\n");
                }
              }
              if (m_applyLimiter[iNode])
	      {
	        stillLimit = true;
	      }
            }
	 }
	 else
	 {
	   stillLimit = true;
	 }

         if (stillLimit || unphysical)
	 {
	   CFLog(VERBOSE, "limiting cell " << m_cell->getID() << "\n");
	   if (iOrder > 1)
	   {
             // compute Pn-1 projection of states
             computeProjStates(m_states2,iOrder-1);
	    
	     for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
             {
	       for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	       {
	         (*((*m_cellStates)[iSol]))[iEq] = m_states2[iSol][iEq];
	       }
             }
             
             ++limitCounter;
             ++orderLimit;
             
             // Update limiting output socket - order reduction limiting (1.0)
             updateLimitingOutput(1.0);
	   }
	   else
	   {
	     executeSlopeLimiter(elemIdx, useMin);
	     ++limitCounter;
	     
	     // Update limiting output socket - slope limiting (2.0)
	     updateLimitingOutput(2.0);
	     
// 	     // only for purposes of printing the shock detector!!
// 	     for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
// 	     {
// 	       (*((*m_cellStates)[iSol]))[0] = 15.0;
//              }
	   }
	 }
	 else
	 {
	   // store the order to which was limited
	   const CFuint prevOrder = static_cast<CFuint>(m_limiterValues[elemIdx][0]+0.5);
	   
	   if (useMin)
	   {
	     m_limiterValues[elemIdx][0] = min(iOrder,prevOrder);
	   }
	   else
	   {
	     m_limiterValues[elemIdx][0] = iOrder;
	   }
	   
	   if (!useMin || prevOrder == iOrder)
	   {
	     if (m_order != iOrder)
	     {
	       ++limitCounter;
	       ++orderLimit;
	       applyChecks(1.0);
	       
	       // Update limiting output socket - order reduction limiting (1.0)
	       updateLimitingOutput(1.0);
// 	       // only for purposes of printing the shock detector!!
// 	       for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
// 	       {
// 	         (*((*m_cellStates)[iSol]))[0] = 15.0;
// 	       }
	     }
	     break;
	   }
	   else if (useMin && prevOrder == 0)
	   {
	     executeSlopeLimiter(elemIdx,useMin);
	     ++limitCounter;
	     
	     // Update limiting output socket - slope limiting (2.0)
	     updateLimitingOutput(2.0);
	     break;
	   }
	 }
       }
       else
       {
	 // store the order to which was limited
	 const CFuint prevOrder = static_cast<CFuint>(m_limiterValues[elemIdx][0]+0.5);
	 
	 if (useMin)
	 {
	   m_limiterValues[elemIdx][0] = min(iOrder,prevOrder);
	 }
	 else
	 {
	   m_limiterValues[elemIdx][0] = iOrder;
	 }
	 
	 if (!useMin || prevOrder == iOrder)
	 {
	   if (iOrder != m_order) 
	   {
	     ++limitCounter;
	     ++orderLimit;
	     
	     // Update limiting output socket - order reduction limiting (1.0)
	     updateLimitingOutput(1.0);
// 	     // only for purposes of printing the shock detector!!
// 	     for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
// 	     {
// 	       (*((*m_cellStates)[iSol]))[0] = 15.0;
// 	     }
	   }
	   break;
	 }
	 else if (useMin && prevOrder == 0)
	 {
	   executeSlopeLimiter(elemIdx,useMin);
	   ++limitCounter;
	   
	   // Update limiting output socket - slope limiting (2.0)
	   updateLimitingOutput(2.0);
	   break;
	 }
       }
     }
      }
//       else
//       {
// 	///@warning for plotting limitFac
// 	for(CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
// 	{
// 	  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
// 	  {
// 	    (*((*m_cellStates)[iSol]))[iEq] = 0.0;
// 	  }
// 	}
//       }

      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
	CFuint stateID = (*m_cellStates)[iSol]->getLocalID();

	RealVector diff = (*((*m_cellStates)[iSol])-m_prevStates[stateID])*(*((*m_cellStates)[iSol])-m_prevStates[stateID]);

	duL2Norm += diff;
	m_prevStates[stateID] = *((*m_cellStates)[iSol]);
	//CFLog(NOTICE, "diff: " << diff << "\n");
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  
  duL2Norm = 0.5*log10(duL2Norm);
  
  if (PE::GetPE().GetRank(nsp) == 0 && iter%m_showrate == 0) 
  {
    CFLog(INFO,"Limits done: " << limitCounter << " of which " << orderLimit << " partial limits\n");
    CFLog(INFO,"L2 norm of the state difference: " << duL2Norm << "\n");
  }

  }
  else
  {      
    applyPrevLimiter();
  }
  
//   CFreal sumLim = 0.0;
//   
//   for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
//   {
//     // get start and end indexes for this type of element
//     const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
//     const CFuint endIdx   = (*elemType)[iElemType].getEndIdx  ();
// 
//     // loop over cells
//     for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
//     {
//       for (CFuint iEq = 0; iEq < m_nbrEqs+1; ++iEq)
//       {
//         sumLim += m_limiterValues[elemIdx][iEq];
//       }
//     }
//   }
//   CFLog(NOTICE,"Sum limiters: " << sumLim << "\n");
  CFTRACEEND;
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::computeProjStates(std::vector< RealVector >& projStates, CFuint order)
{
  if (m_order != 1)
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      RealVector temp(projStates.size());
      for (CFuint iSol = 0; iSol < projStates.size(); ++iSol)
      {
        temp[iSol] = (*((*m_cellStates)[iSol]))[iEq];
      }
      RealVector tempProj(projStates.size());
      tempProj = m_transformationMatrices[order]*temp;
      for (CFuint iSol = 0; iSol < projStates.size(); ++iSol)
      {
        projStates[iSol][iEq] = tempProj[iSol];
      }
    }
  }
  else
  {
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      RealVector temp(projStates.size());
      for (CFuint iSol = 0; iSol < projStates.size(); ++iSol)
      {
        projStates[iSol][iEq] = (*((*m_cellStates)[iSol]))[iEq];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::applyPrevLimiter()
{
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  
  CFuint limitCounter = 0;
  CFuint orderLimit = 0;
  
  // LOOP OVER ELEMENT TYPES, TO SET THE MINIMUM AND MAXIMUM CELL AVERAGED STATES IN THE NODES
  const CFuint nbrElemTypes = elemType->size();
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx  ();

    // set reconstruction data for averaged solutions
    setAvgReconstructionData(iElemType);

    // loop over cells to apply previous limiter
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      // get the nodes in this cell
      m_cellNodes  = m_cell->getNodes();
      
      const RealVector limiterValues = m_limiterValues[elemIdx];
      
      const CFuint limitOrder = static_cast<CFuint>(limiterValues[0]+0.5);
      
      // check if limiting is needed
      if (limitOrder != m_order)
      {
	// check if order reduction is needed or P1 slope limiting
        if (limitOrder != 0)
        {
	  // compute Pn-1 projection of states
          computeProjStates(m_states2,limitOrder);
	    
	  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
	    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	    {
	      (*((*m_cellStates)[iSol]))[iEq] = m_states2[iSol][iEq];
	    }
          }
          
          ++limitCounter;
          ++orderLimit;
          
          applyChecks(1.0);
          
          // Update limiting output socket - order reduction limiting (1.0)
          updateLimitingOutput(1.0);
        }
        else
        {
	  // reconstruct cell averaged state
          reconstructCellAveragedState();
	  
	  // compute P1 projection of states
          computeProjStates(m_statesP1,1);
	  
	  CFreal phiMin = 1.0;
	  
	  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	  {
	  
	    const CFreal phi = limiterValues[iEq+1];
	    
	    if (phi < phiMin)
	    {
	      phiMin = phi;
	    }
	  
	    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
	    {
  	      (*((*m_cellStates)[iSol]))[iEq] = (1.0-phi)*m_cellAvgState[iEq] + phi*m_statesP1[iSol][iEq];
	    }
	  }
	  
	  applyChecks(phiMin);
	  
	  ++limitCounter;
	  
	  // Update limiting output socket - slope limiting (2.0)
	  updateLimitingOutput(2.0);
        }
      }

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  
  // Print limiting statistics for frozen limiter mode
  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  const std::string nsp = this->getMethodData().getNamespace();
  
  if (PE::GetPE().GetRank(nsp) == 0 && iter%m_showrate == 0) 
  {
    CFLog(INFO,"Limits done (frozen): " << limitCounter << " of which " << orderLimit << " partial limits\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::resetNodeNghbrCellAvgStates()
{
  // get the data handles for the minimum and maximum nodal states
  DataHandle< RealVector > nodeNghbCellMinAvgStates = socket_nodeNghbCellMinAvgStates.getDataHandle();
  DataHandle< RealVector > nodeNghbCellMaxAvgStates = socket_nodeNghbCellMaxAvgStates.getDataHandle();

  // get the number of nodes
  const CFuint nbrNodes = nodeNghbCellMinAvgStates.size();
  cf_assert(nbrNodes == nodeNghbCellMaxAvgStates.size());

  // reset the cell averaged states
  for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
  {
    nodeNghbCellMinAvgStates[iNode] = +MathTools::MathConsts::CFrealMax();
    nodeNghbCellMaxAvgStates[iNode] = -MathTools::MathConsts::CFrealMax();
    nodeNghbCellMaxAvgStates[iNode][m_nbrEqs] = 0;
  }
  m_minAvgStateAll = +MathTools::MathConsts::CFrealMax();
  m_maxAvgStateAll = -MathTools::MathConsts::CFrealMax();
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::setAvgReconstructionData(CFuint iElemType)
{
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

  // number of cell corner nodes
  /// @note in the future, hanging nodes should be taken into account here
  m_nbrCornerNodes = frLocalData[iElemType]->getNbrCornerNodes();

    // get solution point local coordinates
  m_cellAvgSolCoefs = frLocalData[iElemType]->getCellAvgSolCoefs();

    // number of solution points
  m_nbrSolPnts = frLocalData[iElemType]->getNbrOfSolPnts();
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::setAllReconstructionData(CFuint iElemType)
{
  // set reconstruction data for averaged solutions
  setAvgReconstructionData(iElemType);

  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

  // get coefficients for cell center derivatives
  m_cellCenterDerivCoefs = frLocalData[iElemType]->getCellCenterDerivCoefs();

  // get the mapped coordinates of the solution points
  m_solPntsLocalCoords = frLocalData[iElemType]->getSolPntsLocalCoords();
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::reconstructCellAveragedState()
{
  //m_cellAvgState = (*m_cellAvgSolCoefs)[0]*(*(*m_cellStates)[0]);
  //for (CFuint iSol = 1; iSol < m_nbrSolPnts; ++iSol)
  //{
    //m_cellAvgState += (*m_cellAvgSolCoefs)[iSol]*(*(*m_cellStates)[iSol]);
  //}
  computeProjStates(m_states2,0);
  m_cellAvgState = 0.0;
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    if(m_cell->getID() == 36)
      {
	CFLog(VERBOSE, "state: " << m_states2[iSol] << "\n");
      }
    m_cellAvgState += m_states2[iSol];
  }
  m_cellAvgState /= m_nbrSolPnts;
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::reconstructCellAveragedVariable(const CFuint iEq)
{
  m_cellAvgState[iEq] = (*m_cellAvgSolCoefs)[0]*(*(*m_cellStates)[0])[iEq];
  for (CFuint iSol = 1; iSol < m_nbrSolPnts; ++iSol)
  {
    m_cellAvgState[iEq] += (*m_cellAvgSolCoefs)[iSol]*(*(*m_cellStates)[iSol])[iEq];
  }
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::computeCellCenterDerivVariable(const CFuint iEq)
{
  for (CFuint iDeriv = 0; iDeriv < m_dim; ++iDeriv)
  {
    m_cellCenterDerivVar[iDeriv] = (*m_cellCenterDerivCoefs)[iDeriv][0]*(*(*m_cellStates)[0])[iEq];
    for (CFuint iSol = 1; iSol < m_nbrSolPnts; ++iSol)
    {
      m_cellCenterDerivVar[iDeriv] += (*m_cellCenterDerivCoefs)[iDeriv][iSol]*(*(*m_cellStates)[iSol])[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::setMinMaxNodeNghbrCellAvgStates()
{
  // get the data handles for the minimum and maximum nodal states
  DataHandle< RealVector > nodeNghbCellMinAvgStates = socket_nodeNghbCellMinAvgStates.getDataHandle();
  DataHandle< RealVector > nodeNghbCellMaxAvgStates = socket_nodeNghbCellMaxAvgStates.getDataHandle();

  // loop over corner nodes, to set the neighbouring min and max states
  for (CFuint iNode = 0; iNode < m_nbrCornerNodes; ++iNode)
  {
    // get node ID
    const CFuint nodeID = (*m_cellNodes)[iNode]->getLocalID();

    // get min and max states
    RealVector& minAvgState = nodeNghbCellMinAvgStates[nodeID];
    RealVector& maxAvgState = nodeNghbCellMaxAvgStates[nodeID];

    // loop over physical variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      // minimum average state
      minAvgState[iEq] = m_cellAvgState[iEq] < minAvgState[iEq] ? m_cellAvgState[iEq] : minAvgState[iEq];

      // maximum average state
      maxAvgState[iEq] = m_cellAvgState[iEq] > maxAvgState[iEq] ? m_cellAvgState[iEq] : maxAvgState[iEq];
    }
    
    // indicate that the node has been checked once
    maxAvgState[m_nbrEqs] += 1;
  }

  // set min and max of all states on the mesh
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
      // minimum average state
    m_minAvgStateAll[iEq] = m_cellAvgState[iEq] < m_minAvgStateAll[iEq] ? m_cellAvgState[iEq] : m_minAvgStateAll[iEq];

      // maximum average state
    m_maxAvgStateAll[iEq] = m_cellAvgState[iEq] > m_maxAvgStateAll[iEq] ? m_cellAvgState[iEq] : m_maxAvgStateAll[iEq];
  }
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::setMinMaxNghbrCellAvgStates()
{
  // get the data handles for the minimum and maximum nodal states
  DataHandle< RealVector > nodeNghbCellMinAvgStates = socket_nodeNghbCellMinAvgStates.getDataHandle();
  DataHandle< RealVector > nodeNghbCellMaxAvgStates = socket_nodeNghbCellMaxAvgStates.getDataHandle();

  // set the minimum and maximum average states in the node neighbours
  const CFuint nodeID = (*m_cellNodes)[0]->getLocalID();
  m_minAvgState = nodeNghbCellMinAvgStates[nodeID];
  m_maxAvgState = nodeNghbCellMaxAvgStates[nodeID];
  for (CFuint iNode = 1; iNode < m_nbrCornerNodes; ++iNode)
  {
    // get node ID
    const CFuint nodeID = (*m_cellNodes)[iNode]->getLocalID();

    // get min and max states in nodes
    const RealVector& noEminAvgState = nodeNghbCellMinAvgStates[nodeID];
    const RealVector& nodeMaxAvgState = nodeNghbCellMaxAvgStates[nodeID];

    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      m_minAvgState[iEq] = m_minAvgState[iEq] < noEminAvgState[iEq] ? m_minAvgState[iEq] : noEminAvgState[iEq];
      m_maxAvgState[iEq] = m_maxAvgState[iEq] > nodeMaxAvgState[iEq] ? m_maxAvgState[iEq] : nodeMaxAvgState[iEq];
    }
  }
  //CFLog(VERBOSE, "min: " << m_minAvgState[0] << "\n");
  //CFLog(VERBOSE, "max: " << m_maxAvgState[0] << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::limitPerturbedStates(CFuint elemIdx, CFuint pertSol, CFuint pertVar)
{
  const RealVector limiterValues = m_limiterValues[elemIdx];
      
  const CFuint limitOrder = static_cast<CFuint>(limiterValues[0]+0.5);
  
  if (limitOrder != m_order)
  {
    // get InnerCells TopologicalRegionSet
    SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

    // get the geodata of the geometric entity builder and set the TRS
    StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
    geoData.trs = cells;
    
    // build the GeometricEntity
    geoData.idx = elemIdx;
    m_cell = m_cellBuilder->buildGE();

    // get the states in this cell
    m_cellStates = m_cell->getStates();
    
    // copy the states in the correct format
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      m_cellStatesBackup[iSol] = (*((*m_cellStates)[iSol]));
    }
  
    if (limitOrder == 1)
    {
      // reconstruct cell averaged state
      reconstructCellAveragedState();
	  
      // compute P1 projection of states
      computeProjStates(m_statesP1,1);
	  
      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      {
        const CFreal phi = limiterValues[iEq+1];
	  
        for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
        {
          (*((*m_cellStates)[iSol]))[iEq] = (1.0-phi)*m_cellAvgState[iEq] + phi*m_statesP1[iSol][iEq];
        }
      }
    }
    else
    {   
      // compute Pn-1 projection of states
      computeProjStates(m_states2,limitOrder);
	    
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
        {
          (*((*m_cellStates)[iSol]))[iEq] = m_states2[iSol][iEq];
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::restoreStates(CFuint elemIdx)
{
  const RealVector limiterValues = m_limiterValues[elemIdx];
      
  const CFuint limitOrder = static_cast<CFuint>(limiterValues[0]+0.5);
  
  if (limitOrder != m_order)
  {
    // get InnerCells TopologicalRegionSet
    SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

    // get the geodata of the geometric entity builder and set the TRS
    StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
    geoData.trs = cells;
    
    // build the GeometricEntity
    geoData.idx = elemIdx;
    m_cell = m_cellBuilder->buildGE();

    // get the states in this cell
    m_cellStates = m_cell->getStates();
    
    // copy the states in the correct format
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      (*((*m_cellStates)[iSol])) = m_cellStatesBackup[iSol];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::computeNodeStates(std::vector< RealVector > states, vector< RealVector >& statesNodes)
{
  // loop over nodes to extrapolate the states to the flux points and get the 
  // discontinuous normal flux in the flux points and the Riemann flux
  for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
  {        
    statesNodes[iNode] = 0.0;

    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      statesNodes[iNode] += (*m_solPolyValsAtNodes)[iNode][iSol]*states[iSol];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::computeFlxPntStates(std::vector< RealVector > states, std::vector< RealVector >& statesFlxPnt)
{
  for (CFuint iFlxPnt = 0; iFlxPnt < m_maxNbrFlxPnts; ++iFlxPnt)
  {        
    statesFlxPnt[iFlxPnt] = 0.0;

    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      statesFlxPnt[iFlxPnt] += (*m_solPolyValsAtFlxPnts)[iFlxPnt][iSol]*states[iSol];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::executeSlopeLimiter(const CFuint elemIdx, const bool useMin)
{
  // get the data handles for the minimum and maximum nodal states
  DataHandle< RealVector > nodeNghbCellMinAvgStates = socket_nodeNghbCellMinAvgStates.getDataHandle();
  DataHandle< RealVector > nodeNghbCellMaxAvgStates = socket_nodeNghbCellMaxAvgStates.getDataHandle();

  CFreal phiMin = 1.0;
	     
  bool comparePhi = false;
	     
  // store the order to which was limited
  if (useMin)
  {
    comparePhi = static_cast<CFuint>(m_limiterValues[elemIdx][0]+0.5) == 0;
  }
  m_limiterValues[elemIdx][0] = 0.0;
	   
  for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
  {
    CFreal phi = 2.0;
	        
    for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
    {
      // get node ID
      const CFuint nodeID = (*m_cellNodes)[iNode]->getLocalID();
	      
      // get min and max states
      RealVector minAvgState = nodeNghbCellMinAvgStates[nodeID];
      RealVector maxAvgState = nodeNghbCellMaxAvgStates[nodeID];
		
      if (fabs(m_cellStatesNodesP1[iNode][iEq] - m_cellAvgState[iEq]) > MathConsts::CFrealEps())
      {
        CFreal r1 = (minAvgState[iEq]-m_cellAvgState[iEq])/(m_cellStatesNodesP1[iNode][iEq] - m_cellAvgState[iEq]);
        CFreal r2 = (maxAvgState[iEq]-m_cellAvgState[iEq])/(m_cellStatesNodesP1[iNode][iEq] - m_cellAvgState[iEq]);
        CFreal r = max(r1,r2);
	if(r < 0.0)
	{
	  CFLog(NOTICE, "r < 0!!\n");
	  CFLog(VERBOSE, "r: " << r << ", r1: " << r1 << ", r2: " << r2 << "\n");
	  CFLog(VERBOSE, "min: " << minAvgState[iEq] << ", max: " << maxAvgState[iEq] << ", av: " << m_cellAvgState[iEq] << ", P1: " << m_cellStatesNodesP1[iNode][iEq] << "\n");
	  r = 0.0;
	}
	CFreal Dm = m_cellStatesNodesP1[iNode][iEq] - m_cellAvgState[iEq];
	CFreal Dp = Dm*r;

	//CFreal epsilon2 = m_tvbLimitFactor*m_cell->computeVolume();
	//CFreal epsilon2 = pow(m_tvbLimitFactor*(m_maxAvgStateAll[iEq]-m_minAvgStateAll[iEq]),2);
	CFreal epsilon2 = m_tvbLimitFactor*pow((maxAvgState[iEq]-minAvgState[iEq]),2.0)/(1.0+(maxAvgState[iEq]-minAvgState[iEq])/(m_tvbLimitFactor*pow(m_cell->computeVolume(),0.75)));
	r = ((Dp*Dp+epsilon2)*Dm + 2.0*Dm*Dm*Dp)/(Dp*Dp + 2.0*Dm*Dm + Dm*Dp + epsilon2)*1.0/Dm;
	phi = min(phi,r);
      }
      else
      {
	phi = min(phi,1.0);
      }
    }
       
    if (useMin && comparePhi)
    {
      phi = min(phi,m_limiterValues[elemIdx][iEq+1]);
    }

    CFLog(VERBOSE, "phi: " << phi << "\n");
    cf_assert(phi > -0.01 && phi < 2.01);
    if (phi < phiMin)
    {
      phiMin = phi;
    }
	      
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      //for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
      //{
      (*((*m_cellStates)[iSol]))[iEq] = (1.0-phi)*m_cellAvgState[iEq] + phi*m_statesP1[iSol][iEq];
      //} 
    }
    if(m_cell->getID() == 36)
    {
      CFLog(VERBOSE, "new states: \n");
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
	CFLog(VERBOSE, (*((*m_cellStates)[iSol])) << "\n");
      }
    }
	       
    // store the limiting value
    m_limiterValues[elemIdx][iEq+1] = phi;
  }
  applyChecks(phiMin);
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::updateLimitingOutput(CFreal limitingType)
{
  DataHandle< CFreal > outputLimiting = socket_outputLimiting.getDataHandle();
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();
    // Use maximum value if multiple types of limiting occur
    outputLimiting[stateID] = max(outputLimiting[stateID], limitingType);
  }
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::setup()
{
  CFAUTOTRACE;
  
  FluxReconstructionSolverCom::setup();

  // get number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim    = PhysicalModelStack::getActive()->getDim ();

  // resize variables
  m_applyLimiterToPhysVar.resize(m_nbrEqs);
  m_cellAvgState.resize(m_nbrEqs);
  m_cellCenterDerivVar.resize(m_dim);
  m_minAvgState.resize(m_nbrEqs);
  m_maxAvgState.resize(m_nbrEqs+1);
  m_minAvgStateAll.resize(m_nbrEqs);
  m_maxAvgStateAll.resize(m_nbrEqs);

  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  const CFuint nbrElemTypes = frLocalData.size();
  cf_assert(nbrElemTypes > 0);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
  m_order = static_cast<CFuint>(order);

  
  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtNodes = frLocalData[0]->getCoefSolPolyInNodes();
  
  m_nbrNodesElem = m_solPolyValsAtNodes->size();
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();

  // get the maximum number of flux points
  m_maxNbrFlxPnts = 0;
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbrFlxPnts = frLocalData[iElemType]->getNbrOfFlxPnts();
    m_maxNbrFlxPnts = m_maxNbrFlxPnts > nbrFlxPnts ? m_maxNbrFlxPnts : nbrFlxPnts;
    
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx  ();

    // loop over cells to apply previous limiter
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      RealVector temp(m_nbrEqs+1);
      temp = 1.0;
      temp[0] = m_order;
      
      m_limiterValues.push_back(temp);
      
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        RealVector temp2(m_nbrEqs);
        temp2 = 0.0;
        m_prevStates.push_back(temp2);
      }
    }
  }

  // get the number of nodes in the mesh
  const CFuint nbrNodes = MeshDataStack::getActive()->getNbNodes();

  // resize the sockets for the minimum and maximum states in the node neighbouring cells
  /// @todo this socket is too big, only the first order nodes (in the cell corners, and later on, hanging nodes)
  /// have to be taken into account
  DataHandle< RealVector > nodeNghbCellMinAvgStates = socket_nodeNghbCellMinAvgStates.getDataHandle();
  DataHandle< RealVector > nodeNghbCellMaxAvgStates = socket_nodeNghbCellMaxAvgStates.getDataHandle();
  nodeNghbCellMinAvgStates.resize(nbrNodes);
  nodeNghbCellMaxAvgStates.resize(nbrNodes);
  for (CFuint iNode = 0; iNode < nbrNodes; ++iNode)
  {
    nodeNghbCellMinAvgStates[iNode].resize(m_nbrEqs);
    nodeNghbCellMaxAvgStates[iNode].resize(m_nbrEqs+1);
  }
  
  // get the number of cells in the mesh to compute total number of solution points
  const CFuint nbrCells = (*elemType)[0].getEndIdx();
  
  // resize socket for limiting output (one value per solution point)
  DataHandle< CFreal > outputLimiting = socket_outputLimiting.getDataHandle();
  const CFuint nbStates = nbrCells * m_nbrSolPnts;
  outputLimiting.resize(nbStates);

  m_lengthScaleExp = 2.0/static_cast<CFreal>(m_dim);
  
  for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
  {
    RealVector temp(m_nbrEqs);
    RealVector temp2(m_nbrEqs);
    m_cellStatesNodes.push_back(temp);
    m_cellStatesNodesP1.push_back(temp2);
  }
  
  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();
  
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  {
    RealVector temp(m_nbrEqs);
    m_cellStatesFlxPnt.push_back(temp);
  }
  
  SafePtr<RealMatrix> vdm = frLocalData[0]->getVandermondeMatrix();
  
  SafePtr<RealMatrix> vdmInv = frLocalData[0]->getVandermondeMatrixInv();
  
  const CFuint nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();
  
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    RealVector temp(m_nbrEqs);
    temp = 0.0;
    m_statesP1.push_back(temp);
    m_states2.push_back(temp);
    m_cellStatesBackup.push_back(temp);
  }
  
  m_applyLimiter.resize(m_nbrNodesElem);
  
  // Create transformation matrices for order reduction
  // These matrices project the solution from higher order to lower order polynomial spaces
  // For order reduction from P to P-k, we keep only the first (P-k+1) basis functions
  for (CFuint iOrder = 0; iOrder < m_order; ++iOrder)
  {
    RealMatrix temp(nbrSolPnts,nbrSolPnts);
    temp = 0.0;
    
    // Get the number of solution points for the reduced order (iOrder+1)
    // - Triangular: (P+1)*(P+2)/2
    // - Quadrilateral: (P+1)^2  
    // - Tetrahedral: (P+1)*(P+2)*(P+3)/6
    // - Prism: (P+1)*(P+1)*(P+2)/2
    // - Hexahedral: (P+1)^3
    CFuint nbrSolPntsReduced;
    if (m_dim == DIM_2D && frLocalData[0]->getShape() == CFGeoShape::TRIAG) {
      // For triangular elements: (P+1)*(P+2)/2
      nbrSolPntsReduced = (iOrder+1)*(iOrder+2)/2;
    } else if (m_dim == DIM_2D && frLocalData[0]->getShape() == CFGeoShape::QUAD) {
      // For quadrilateral elements: (P+1)^2
      nbrSolPntsReduced = (iOrder+1)*(iOrder+1);
    } else if (m_dim == DIM_3D && frLocalData[0]->getShape() == CFGeoShape::TETRA) {
      // For tetrahedral elements: (P+1)*(P+2)*(P+3)/6
      nbrSolPntsReduced = (iOrder+1)*(iOrder+2)*(iOrder+3)/6;
    } else if (m_dim == DIM_3D && frLocalData[0]->getShape() == CFGeoShape::PRISM) {
      // For prism elements: (P+1)*(P+1)*(P+2)/2
      nbrSolPntsReduced = (iOrder+1)*(iOrder+1)*(iOrder+2)/2;
    } else if (m_dim == DIM_3D && frLocalData[0]->getShape() == CFGeoShape::HEXA) {
      // For hexahedral elements: (P+1)^3
      nbrSolPntsReduced = (iOrder+1)*(iOrder+1)*(iOrder+1);
    } else {
      throw Common::NotImplementedException(FromHere(), "MLPLimiter::setup() - unsupported element type for order reduction");
    }
    
    // Create identity matrix for the first nbrSolPntsReduced modes (keeps lower-order modes)
    // Higher-order modes are set to zero, effectively reducing the polynomial order
    for (CFuint idx = 0; idx < nbrSolPntsReduced; ++idx)
    {
      temp(idx,idx) = 1.0;
    }
    m_transformationMatrices.push_back((*vdm)*temp*(*vdmInv));
  }
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::unsetup()
{
  CFAUTOTRACE;

  DataHandle< RealVector > nodeNghbCellMinAvgStates = socket_nodeNghbCellMinAvgStates.getDataHandle();
  DataHandle< RealVector > nodeNghbCellMaxAvgStates = socket_nodeNghbCellMaxAvgStates.getDataHandle();
  nodeNghbCellMinAvgStates.resize(0);
  nodeNghbCellMaxAvgStates.resize(0);
  
  FluxReconstructionSolverCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSource> >
MLPLimiter::providesSockets()
{
  vector<SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_nodeNghbCellMinAvgStates);
  result.push_back(&socket_nodeNghbCellMaxAvgStates);
  result.push_back(&socket_outputLimiting);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
