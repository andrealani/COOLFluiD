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
  options.addConfigOption< CFreal >("LimFactor","K factor in the MLP-u2 limiting part.");
}

//////////////////////////////////////////////////////////////////////////////

MLPLimiter::MLPLimiter(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_nodeNghbCellMinAvgStates("nodeNghbCellMinAvgStates"),
  socket_nodeNghbCellMaxAvgStates("nodeNghbCellMaxAvgStates"),
  m_cellBuilder(CFNULL),
  m_solPolyValsAtNodes(CFNULL),
  m_solPolyDerivAtNodes(CFNULL),
  m_cell(),
  m_cellStates(),
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
  m_nbrFlxPnts(),
  m_nbrSolPnts(),
  m_nbrCornerNodes(),
  m_applyLimiter(),
  m_applyLimiterToPhysVar(),
  m_tvbLimitFactor(),
  m_lengthScaleExp(),
  m_cellStatesNodes(),
  m_cellStatesNodesP1(),
  m_transformationMatrices(),
  m_statesP1(),
  m_states2(),
  m_nbrNodesElem(),
  m_order(),
  m_solPolyValsAtFlxPnts(CFNULL),
  m_cellStatesFlxPnt()
{
  addConfigOptionsTo(this);

  m_tvbLimitFactor = 0.0;
  setParameter( "LimFactor", &m_tvbLimitFactor );
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

  // LOOP OVER ELEMENT TYPES, TO SET THE MINIMUM AND MAXIMUM CELL AVERAGED STATES IN THE NODES
  const CFuint nbrElemTypes = elemType->size();
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
  CFuint MLPuLimit = 0;
  CFuint beforeExt = 0;

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
      m_cellNodes  = m_cell->getNodes ();

      // set the min and max neighbouring cell averaged states
      setMinMaxNghbrCellAvgStates();
      
      // compute P1 projection of states
      computeProjStates(m_statesP1,1);

      // compute the P1 states in the nodes
      computeNodeStates(m_statesP1,m_cellStatesNodesP1);

      // reconstruct cell averaged state
      reconstructCellAveragedState();
      
      if(m_cell->getID() == 1876)
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

      bool limitFlag = false;

      // check if P1u is within vertex limiting averages
      for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
      {
        m_applyLimiter[iNode] = m_cellStatesNodesP1[iNode][0] <  1.01*m_minAvgState[0] ||
                                m_cellStatesNodesP1[iNode][0] >  0.99*m_maxAvgState[0];
	if(m_cell->getID() == 1876)
        {
	  bool tempBool = m_cellStatesNodesP1[iNode][0] <  m_minAvgState[0];
	  CFLog(VERBOSE, "min? " << tempBool << "\n");
	  tempBool = m_cellStatesNodesP1[iNode][0] >  m_maxAvgState[0];
	  CFLog(VERBOSE, "max? " << tempBool << "\n");
	}
	if (m_applyLimiter[iNode])
	{
	  limitFlag = true;
	  CFLog(VERBOSE, "prelimiting cell " << m_cell->getID() << "\n");
	  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
	  {
	    //(*((*m_cellStates)[iSol]))[1] = 5.0;
	  }
	}
      }
      
      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
      {
        m_states2[iSol] = (*((*m_cellStates)[iSol]));
      }
        
      // compute the states in the nodes
      computeNodeStates(m_states2,m_cellStatesNodes);
      
      bool unphysical = !checkPhysicality();
      
      if (limitFlag || unphysical)
      {
	++MLPuLimit;
      for (CFuint iOrder = m_order; iOrder > 0; --iOrder)
      {
	
	if(m_cell->getID() == 1876)
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
 
	// check for cte area/smooth extrema
        for (CFuint iNode = 0; iNode < m_nbrNodesElem; ++iNode)
        {
	  if (m_applyLimiter[iNode])
	  {
	    // check if we are in constant area of sol
	    m_applyLimiter[iNode] = !( fabs(m_cellStatesNodes[iNode][0] - m_cellAvgState[0]) < 1e-3 * m_cellAvgState[0]);// ,m_cell->computeVolume()) );
	    
	    // get node ID
            const CFuint nodeID = (*m_cellNodes)[iNode]->getLocalID();
	      
	    // get min and max states
            RealVector minAvgState = nodeNghbCellMinAvgStates[nodeID];
            RealVector maxAvgState = nodeNghbCellMaxAvgStates[nodeID];
	    
	    if (m_applyLimiter[iNode])
            {
	      ++beforeExt;
              // check if we are in a smooth local maximum
	      bool cond1 = m_cellStatesNodesP1[iNode][0] - m_cellAvgState[0] > 0.0;
	      bool cond2 = m_cellStatesNodes[iNode][0] - m_cellStatesNodesP1[iNode][0] < 0.0;
	      bool cond3 = m_cellStatesNodes[iNode][0] > minAvgState[0];
              m_applyLimiter[iNode] = !( cond1 && cond2 && cond3 );	
	      CFLog(VERBOSE, "max cond: " << cond1 << ", " << cond2 << ", " << cond3 << "\n");
            }
      
            if(m_applyLimiter[iNode])
            {
	      // check if we are in a smooth local minimum
	      bool cond1 = m_cellStatesNodesP1[iNode][0] - m_cellAvgState[0] < 0.0;
	      bool cond2 = m_cellStatesNodes[iNode][0] - m_cellStatesNodesP1[iNode][0] > 0.0;
	      bool cond3 = m_cellStatesNodes[iNode][0] < maxAvgState[0];
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
	  if (iOrder > 2)
	  {
            // compute P1 projection of states
            computeProjStates(m_states2,iOrder-1);
	    
	    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
            {
	      for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	      {
	        (*((*m_cellStates)[iSol]))[iEq] = m_states2[iSol][iEq];
	      }
            }
	  }
	  else
	  {
	    CFreal phi = 1.0;
	    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	    {
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
		    CFLog(VERBOSE, "r: " << r << ", r1: " << r1 << ", r2: " << r2 << "\n");
		    CFLog(VERBOSE, "min: " << minAvgState[iEq] << ", max: " << maxAvgState[iEq] << ", av: " << m_cellAvgState[iEq] << ", P1: " << m_cellStatesNodesP1[iNode][iEq] << "\n");
		    r = 0.0;
		  }
		  CFreal Dm = m_cellStatesNodesP1[iNode][iEq] - m_cellAvgState[iEq];
		  CFreal Dp = 0.0;
		  if (r1 > r2)
		  {
		    Dp = minAvgState[iEq]-m_cellAvgState[iEq];
		  }
		  else
		  {
		    Dp = maxAvgState[iEq]-m_cellAvgState[iEq];
		  }
		  //CFreal epsilon2 = m_tvbLimitFactor*m_cell->computeVolume();
		  //CFreal epsilon2 = pow(m_tvbLimitFactor*(m_maxAvgStateAll[iEq]-m_minAvgStateAll[iEq]),2);
		  CFreal epsilon2 = m_tvbLimitFactor*pow((maxAvgState[iEq]-minAvgState[iEq]),2)/(1.0+(maxAvgState[iEq]-minAvgState[iEq])/(m_tvbLimitFactor*pow(m_cell->computeVolume(),0.5)));
		  //r = ((Dp*Dp+epsilon2)*Dm + 2.0*Dm*Dm*Dp)/(Dp*Dp + 2.0*Dm*Dm + Dm*Dp + epsilon2)*1.0/Dm;
		  phi = min(phi,r);
		}
	      }
	    
	      CFLog(VERBOSE, "phi: " << phi << "\n");
	      cf_assert(phi > -0.01 && phi < 1.01);
	      for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
	      {
	        for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
	        {
  	          (*((*m_cellStates)[iSol]))[iEq] = (1.0-phi)*m_cellAvgState[iEq] + phi*m_statesP1[iSol][iEq];
	        }
	      }
	    }
	    applyChecks(phi);
	    
	  }
	  ++limitCounter;
	  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
	  {
	    //(*((*m_cellStates)[iSol]))[0] = 5.0;
	  }
	}
	else
	{
	  if (m_order != iOrder)
	  {
	    ++limitCounter;
	    applyChecks(1.0);
	  }
	  break;
	}
        

      // apply limiter if necessary
     // if (m_applyLimiter)
      //{
	//CFLog(VERBOSE,"Limits done: " << limitCounter << "\n");
	
//	for (CFuint iState = 0; iState < m_cellStates->size(); ++iState)
//	{
//	  CFLog(VERBOSE, "state " << iState << " : " << *((*m_cellStates)[iState]) << "\n");
//	}
//	CFLog(VERBOSE, "average state : " << m_cellAvgState << "\n");
//
  //      limitStates();
	//++limitCounter;
      //}
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

      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  
  CFLog(NOTICE,"MLPu detected: " << MLPuLimit << "\n");
  CFLog(NOTICE,"before extrema detector: " << beforeExt << "\n");
  CFLog(NOTICE,"Limits done: " << limitCounter << "\n");

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
    if(m_cell->getID() == 1876)
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

void MLPLimiter::setLimitBooleans()
{
//   // compute length scale factor
//   const CFreal lengthScaleFactor = m_tvbLimitFactor*pow(m_cell->computeVolume(),m_lengthScaleExp);
// 
//   if (m_applyLimiter)
//   {
//     // check if we are in a constant part of u
//     for (CFuint iNode = 0; iNode < m_nbrNodesElem && m_applyLimiter; ++iNode)
//     {
//       m_applyLimiter = !( fabs(m_cellStatesNodes[iNode][0] - m_cellAvgState[0]) < 1e-3 * m_cellAvgState[0] );	     
//     }
//     if (m_applyLimiter)
//     {
//       // check if we are in a smooth local maximum
//       for (CFuint iNode = 0; iNode < m_nbrNodesElem && m_applyLimiter; ++iNode)
//       {
// 	bool cond1 = m_cellStatesNodesP1[iNode][0] - m_cellAvgState[0] > 0.0;
// 	bool cond2 = m_cellStatesNodes[iNode][0] - m_cellStatesNodesP1[iNode][0] < 0.0;
// 	bool cond3 = m_cellStatesNodes[iNode][0] > m_minAvgStateAll[0];
//         m_applyLimiter = !( cond1 && cond2 && cond3 );	     
//       }
//       
//       if(m_applyLimiter)
//       {
// 	// check if we are in a smooth local minimum
//         for (CFuint iNode = 0; iNode < m_nbrNodesElem && m_applyLimiter; ++iNode)
//         {
// 	  bool cond1 = m_cellStatesNodesP1[iNode][0] - m_cellAvgState[0] < 0.0;
// 	  bool cond2 = m_cellStatesNodes[iNode][0] - m_cellStatesNodesP1[iNode][0] > 0.0;
// 	  bool cond3 = m_cellStatesNodes[iNode][0] < m_maxAvgStateAll[0];
//           m_applyLimiter = !( cond1 && cond2 && cond3 );	     
//         }
//       }
//     }
//   }
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::limitStates()
{
//   for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
//   {
//     for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
//     {
//       
//     }
//   }
}

//////////////////////////////////////////////////////////////////////////////

void MLPLimiter::computeNodeStates(std::vector< RealVector > states, vector< RealVector >& statesNodes)
{
  // loop over nodes to extrapolate the states to the flux points and get the 
  // discontinuous normal flux in the flux points and the Riemann flux
  for (CFuint iFlxPnt = 0; iFlxPnt < m_nbrNodesElem; ++iFlxPnt)
  {        
    statesNodes[iFlxPnt] = 0.0;

    // extrapolate the left and right states to the flx pnts
    for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
    {
      statesNodes[iFlxPnt] += (*m_solPolyValsAtNodes)[iFlxPnt][iSol]*states[iSol];
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

void MLPLimiter::setup()
{
  CFAUTOTRACE;

  // get number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim    = PhysicalModelStack::getActive()->getDim ();

  // resize variables
  m_applyLimiterToPhysVar.resize(m_nbrEqs);
  m_cellAvgState.resize(m_nbrEqs);
  m_cellCenterDerivVar.resize(m_dim);
  m_minAvgState.resize(m_nbrEqs);
  m_maxAvgState.resize(m_nbrEqs);
  m_minAvgStateAll.resize(m_nbrEqs);
  m_maxAvgStateAll.resize(m_nbrEqs);

  // get cell builder
  m_cellBuilder = getMethodData().getStdTrsGeoBuilder();

  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  const CFuint nbrElemTypes = frLocalData.size();
  cf_assert(nbrElemTypes > 0);
  
  
  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtNodes = frLocalData[0]->getCoefSolPolyInNodes();
  
  m_nbrNodesElem = m_solPolyValsAtNodes->size();
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();

  // get the maximum number of flux points
  m_maxNbrFlxPnts = 0;
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbrFlxPnts = frLocalData[iElemType]->getNbrOfFlxPnts();
    m_maxNbrFlxPnts = m_maxNbrFlxPnts > nbrFlxPnts ? m_maxNbrFlxPnts : nbrFlxPnts;
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
    nodeNghbCellMaxAvgStates[iNode].resize(m_nbrEqs);
  }

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
  }
  
  m_applyLimiter.resize(m_nbrNodesElem);
  
  const CFPolyOrder::Type order = frLocalData[0]->getPolyOrder();
  
  m_order = static_cast<CFuint>(order);
  
  for (CFuint iOrder = 0; iOrder < order; ++iOrder)
  {
    RealMatrix temp(nbrSolPnts,nbrSolPnts);
    temp = 0.0;
    for (CFuint idx = 0; idx < (iOrder+1)*(iOrder+1); ++idx)
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
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSource> >
MLPLimiter::providesSockets()
{
  vector<SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_nodeNghbCellMinAvgStates           );
  result.push_back(&socket_nodeNghbCellMaxAvgStates           );

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
