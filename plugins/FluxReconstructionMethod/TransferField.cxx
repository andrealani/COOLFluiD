#include "Framework/CFSide.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/TransferField.hh"
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

MethodCommandProvider<TransferField, FluxReconstructionSolverData, FluxReconstructionModule> TransferFieldFluxReconstructionProvider("TransferField");

//////////////////////////////////////////////////////////////////////////////

void TransferField::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFuint> > ("VariableIDs", "IDs of the variable to be assigned to the newly computed field.");
  
  options.addConfigOption< string > ("OtherNamespace", "Name of the other namespace (providing the potential).");
  
  options.addConfigOption< CFreal > ("InterRadius", "Radius corresponding to the internal boundary between donor and current grids (<= 0 assumes one mesh).");
  
  options.addConfigOption< CFreal > ("DeltaSelection", "Distance within which points in the smaller mesh are selected.");
  
  options.addConfigOption< bool > ("UsePFSSBInit", "Flag to use the PFSS B solution as initialization for B (default true).");
}

//////////////////////////////////////////////////////////////////////////////

TransferField::TransferField(const std::string& name) :
  FluxReconstructionSolverCom(name),
  m_cellBuilder(CFNULL),
  m_cell(),
  m_cellStates(),
  m_flxPntsRecCoefs(),
  m_allFlxPntIdxs(),
  m_solPntsLocalCoords(),
  m_nbrEqs(),
  m_dim(),
  m_nbrFlxPnts(),
  m_nbrSolPnts(),
  m_maxNbrFlxPnts(),
  m_solPolyValsAtNodes(CFNULL),
  m_solPolyDerivAtNodes(CFNULL),
  m_order(),
  m_solPolyValsAtFlxPnts(CFNULL),
  m_cellStatesFlxPnt(),
  socket_otherStates("states"),
  socket_states("states"),
  socket_pastStates("pastStates"),
  m_applyProcessing(true)
{
  addConfigOptionsTo(this);

  m_variableIDs = vector<CFuint>();
  setParameter("VariableIDs",&m_variableIDs);
  
  m_otherNamespace = "";
  setParameter("OtherNamespace", &m_otherNamespace);
  
  m_interRadius = -1.0;
  setParameter("InterRadius", &m_interRadius);

  m_deltaSelection = 0.0;
  setParameter("DeltaSelection", &m_deltaSelection);
  
  m_usePFSSBInit = true;
  setParameter("UsePFSSBInit", &m_usePFSSBInit);
}

//////////////////////////////////////////////////////////////////////////////

TransferField::~TransferField()
{
}

//////////////////////////////////////////////////////////////////////////////

void TransferField::configure ( Config::ConfigArgs& args )
{
  FluxReconstructionSolverCom::configure(args);
  
  cf_assert(m_otherNamespace != "");
  CFLog(VERBOSE, "TransferField::configure() => m_otherNamespace = " << m_otherNamespace << "\n");
  socket_otherStates.setDataSocketNamespace(m_otherNamespace);
}

//////////////////////////////////////////////////////////////////////////////

void TransferField::execute()
{  
  if (m_applyProcessing)
  {
      
  CFLog(INFO, "TransferField::execute() => START\n");

  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get the geodata of the geometric entity builder and set the TRS
  StdTrsGeoBuilder::GeoData& geoData = m_cellBuilder->getDataGE();
  geoData.trs = cells;
  
  const CFuint nbrElemTypes = elemType->size();
  
  // get the state socket
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  
  // LOOP OVER ELEMENT TYPES
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    // get start and end indexes for this type of element
    const CFuint startIdx = (*elemType)[iElemType].getStartIdx();
    const CFuint endIdx   = (*elemType)[iElemType].getEndIdx  ();
    
    // get the local FR data
    vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();

    // number of solution points
    m_nbrSolPnts = frLocalData[iElemType]->getNbrOfSolPnts();

    // get the mapped coordinates of the solution points
    m_solPntsLocalCoords = frLocalData[iElemType]->getSolPntsLocalCoords();

    // loop over cells to set the varIDs of the cells' states to the ones from otherNameSpace
    for (CFuint elemIdx = startIdx; elemIdx < endIdx; ++elemIdx)
    {
      // build the GeometricEntity
      geoData.idx = elemIdx;
      m_cell = m_cellBuilder->buildGE();

      // get the states in this cell
      m_cellStates = m_cell->getStates();

      if (SubSystemStatusStack::getActive()->getNbIter() >= 1) 
      {
        Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance(SubSystemStatusStack::getCurrentName()).getNamespace(m_otherNamespace);
        Common::SafePtr<SubSystemStatus> otherSubSystemStatus = SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);

        DataHandle<State*, GLOBAL> otherStates = socket_otherStates.getDataHandle();
//    // States on the extended corona MHD mesh
//    DataHandle<State*, GLOBAL> states = socket_states.getDataHandle(); // LARGER mesh
//    
//    const CFuint nbStates = states.size();
        cf_assert(m_dim >= 2);
        cf_assert(m_variableIDs.size() >= 2);
        
        const CFuint xVar = m_variableIDs[0];
        cf_assert(xVar < m_nbrEqs);
        cf_assert(xVar == 4);
        
        const CFuint yVar = m_variableIDs[1];
        cf_assert(yVar < m_nbrEqs);
        cf_assert(yVar == 5);
        
        const CFuint zVar = (m_dim == DIM_3D) ? m_variableIDs[2] : 0;
        cf_assert(zVar < m_nbrEqs);
        cf_assert(zVar == 6);
    
        //Stopwatch<WallTime> stp;
        //stp.start();
    
        if (m_interRadius <= 0.0) 
        {
          //CFLog(INFO, "TransferField::execute() => transferring field\n");
          
          for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
          {
            const CFuint stateID = (*m_cellStates)[iSol]->getLocalID();
              
            (*(*m_cellStates)[iSol])[xVar] = (*otherStates[stateID])[1]; 
            (*(*m_cellStates)[iSol])[yVar] = (*otherStates[stateID])[2];
            if (m_dim == 3) (*(*m_cellStates)[iSol])[zVar] = (*otherStates[stateID])[3];    
            
            (*states[stateID])[xVar] = (*otherStates[stateID])[1];
            (*states[stateID])[yVar] = (*otherStates[stateID])[2];
            if (m_dim == 3) (*states[stateID])[zVar] = (*otherStates[stateID])[3]; 
            
            if (m_usePFSSBInit)
            {
              (*pastStates[stateID])[xVar] = (*otherStates[stateID])[1];
              (*pastStates[stateID])[yVar] = (*otherStates[stateID])[2];
              if (m_dim == 3) (*pastStates[stateID])[zVar] = (*otherStates[stateID])[3]; 
            }
          }
        }
        else 
        {
          CFLog(INFO, "TransferField::execute() => Extrapolation to different radius mesh (InterRadius > 0) not implemented for FR!!\n");
          cf_assert(false);
//      DataHandle<State*, GLOBAL> otherStates = socket_otherStates.getDataHandle(); // SMALLER mesh
//      
//      // loop over internal faces
//      //   if face has all its nodes radius < (m_interRadius + eps) and > (m_interRadius - eps)
//      //     store the face states (or stateIDs) whose radius < m_interRadius as interStates
//      //     those states will be those within which looking for the matching "fronteer" otherStates 
//      vector<State*> interStates;
//      interStates.reserve(otherStates.size()/4); // rough estimation
//      
//      SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");
//      SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> > faceBuilder = 
//	this->getMethodData().getFaceTrsGeoBuilder();
//      FaceTrsGeoBuilder::GeoData& geoData = faceBuilder->getDataGE();
//      SafePtr<FaceTrsGeoBuilder> faceBuilderPtr = faceBuilder->getGeoBuilder();
//      faceBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
//      geoData.trs = faces;
//      geoData.isBFace = false;
//
//      cf_assert(m_deltaSelection > 0.);
//      const CFreal rMax = m_interRadius + m_deltaSelection;
//      const CFreal rMin = m_interRadius - m_deltaSelection;
//      
//      const CFuint nbFaces = faces->getLocalNbGeoEnts();
//      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
//	geoData.idx = iFace;
//	const GeometricEntity *const face = faceBuilder->buildGE();
//	const vector<Node*>& nodesInFace = face->getNodes();
//	const CFuint nbNodesInFace = nodesInFace.size();
//	cf_assert(nbNodesInFace > 1);
//	CFuint countInterNodes = 0;
//	CFreal radius = 0.;
//	for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
//	  radius = nodesInFace[iNode]->norm2();
//	  if (radius < rMax && radius > rMin) {
//	    countInterNodes++;
//	  }
//	}
//	if (countInterNodes == nbNodesInFace) {
//	  // store state on the inner side of the internal surface
//	  const CFreal radius0 = face->getState(0)->getCoordinates().norm2();
//	  const CFreal radius1 = face->getState(1)->getCoordinates().norm2();
//	  (radius0 < radius1) ? interStates.push_back(face->getState(0)) : interStates.push_back(face->getState(1)); 
//	  CFLog(DEBUG_MIN, "TransferField::execute() => #" << interStates.size() <<  " face with radius [" << radius << "] detected\n");
//	}
//	faceBuilder->releaseGE();
//      }
//
//      const CFreal maxFloat = std::numeric_limits<double>::max();
//      
//      // find the closest otherStates for each of the selected interStates
//      const CFuint nbOtherStates = otherStates.size();
//      const CFuint nbInterStates = interStates.size();
//      vector<State*> interOtherStates(nbInterStates); // internal states (smaller mesh) attached to m_interRadius boundary
//      vector<CFreal> distMin(nbInterStates, maxFloat);
//      for (CFuint iState = 0; iState < nbInterStates; ++iState) {
//	const Node& coord = interStates[iState]->getCoordinates();
//	for (CFuint jState = 0; jState < nbOtherStates; ++jState) {
//	  const CFreal dist2 = MathFunctions::getSquaredDistance(coord, otherStates[jState]->getCoordinates());
//	  if (dist2 < distMin[iState]) {
//	    distMin[iState] = dist2;
//	    interOtherStates[iState] = otherStates[jState];
//	  }
//	}
//	
//	CFLog(DEBUG_MIN, "TransferField::execute() => distMin[" << iState << "] = " << std::sqrt(distMin[iState]) << "\n");
//      }
//      
//      // count the external states and store the corresponding state IDs in a list
//      vector<State*> externalStates;
//      externalStates.reserve(nbStates); // size overestimated to avoid reallocation
//      vector<State*> internalStates;
//      internalStates.reserve(nbStates); // size overestimated to avoid reallocation
//      
//      const CFreal radiusPFSS = m_interRadius + m_deltaSelection;
//      for (CFuint iState = 0; iState < nbStates; ++iState) {
//	if (states[iState]->getCoordinates().norm2() > radiusPFSS) {
//	  cf_assert(iState == states[iState]->getLocalID());
//	  externalStates.push_back(states[iState]);
//	}
//	else {
//	  internalStates.push_back(states[iState]);
//	}
//      }
//      cf_assert(externalStates.size() + internalStates.size() == nbStates);
//      
//      CFLog(INFO, "TransferField::execute() => detected [" << externalStates.size() << "] external states\n");
//      
//      // loop over all otherStates
//      //   if state radius > m_interRadius
//      //     loop over all interStates to find the closest state and store the minimal distance
//      //   else increment internal state counter (which will have to match otherStates.size() at the end)
//      const CFuint nbExternalStates = externalStates.size();
//
//      vector<State*> closestStates(nbExternalStates);
//      RealVector closestStateDist(maxFloat, nbExternalStates);
//      
//      for (CFuint iState = 0; iState < nbExternalStates; ++iState) {
//	const Node& coord = externalStates[iState]->getCoordinates();
//	for (CFuint jState = 0; jState < nbInterStates; ++jState) {
//	  State* const intState = interOtherStates[jState];
//	  // use square distance for comparison since it is faster to compute 
//	  const CFreal dist2 =
//	    MathFunctions::getSquaredDistance(coord, intState->getCoordinates());
//	  if (dist2 < closestStateDist[iState]) {
//	    closestStateDist[iState] = dist2;
//	    closestStates[iState] = intState;
//	  }
//	}
//      }            
//      closestStateDist = sqrt(closestStateDist);
//      
//      const CFuint nbInternalStates = internalStates.size();
//      vector<State*> closestInternalStates(nbInternalStates);
//
//      for (CFuint iState = 0; iState < nbInternalStates; ++iState) {
//	const Node& coord = internalStates[iState]->getCoordinates();
//	CFreal distMin = maxFloat;
//	for (CFuint jState = 0; jState < nbOtherStates; ++jState) {
//	  State* const intState = otherStates[jState];
//	  // use square distance for comparison since it is faster to compute 
//	  const CFreal dist2 =
//	    MathFunctions::getSquaredDistance(coord, intState->getCoordinates());
//	  if (dist2 < distMin) {
//	    distMin = dist2;
//	    closestInternalStates[iState] = intState;
//	  }
//	}
//      } 
//      
//      // from here on:
//      // 1) externalStates[k]   gives the State from LARGER mesh with radius > m_interRadius 
//      // 2) closestStates[k]    gives the State from SMALLER mesh to be used for extrapolation
//      // 3) closestStateDist[k] gives the distance between 1) and 2)
//      // 4) closestInternalStates[k]->getLocalID() gives the state ID to be used to fetch the gradient in uX, uY, uZ
//      CFLog(VERBOSE, "size internalStates => " << internalStates[0]->size() << "\n");
//      CFLog(VERBOSE, "size closestInternalStates  => " << closestInternalStates[0]->size() << "\n");
//
//      // array to keep track of the state IDs that correspond to internal and external states
//      vector<bool> flag(nbStates, false);
//      
//      // for (CFuint iState = 0; iState < nbExternalStates; ++iState) {
//      //	const CFuint estateID = externalStates[iState]->getLocalID();
//      //	const CFuint ostateID = closestStates[iState]->getLocalID();
//      /// use ux_i[ostateID] to compute Br, Btheta, Bphi
//      //	flag[estateID] = true;
//      //	(*states[estateID])[xVar] = ux[ostateID];
//      //	(*states[estateID])[yVar] = uy[ostateID];
//      //	if (dim == DIM_3D) {
//      //	  (*states[estateID])[zVar] = uz[ostateID];
//      //	}
//      // }
//
//      // from here on:
//      // 1) internalStates[k]   gives the State from LARGER mesh with radius < m_interRadius 
//      // 2) closestStates[k]    gives the State from SMALLER mesh to be used for extrapolation
//      // 3) closestStateDist[k] gives the distance between 1) and 2)
//      // 4) closestStates[k]->getLocalID() gives the state ID to be used to fetch the gradient in uX, uY, uZ
//            
//      CFLog(VERBOSE, "size internalStates => " << internalStates[0]->size() << "\n");
//      CFLog(VERBOSE, "size closestStates  => " << closestStates[0]->size() << "\n");
//      CFLog(VERBOSE, "size internalStates => " << internalStates[0]->size() << "\n");
//      CFLog(VERBOSE, "size closestInternalStates  => " << closestInternalStates[0]->size() << "\n");
//      
//      for (CFuint iState = 0; iState < nbInternalStates; ++iState) {
//	const CFuint istateID = internalStates[iState]->getLocalID();
//	const CFuint ostateID = closestInternalStates[iState]->getLocalID();
//	flag[istateID] = true;
//	/// use ux_i[ostateID] to compute Br, Btheta, Bphi
//	(*states[istateID])[xVar] = ux[ostateID];
//	(*states[istateID])[yVar] = uy[ostateID];
//	if (dim == DIM_3D) {
//	  (*states[istateID])[zVar] = uz[ostateID];
//	}
//      }
//      
//      //<<<<<<<<<<< HERE FOLLOWS THE EXTRAPOLATION ONTO THE OUTER MESH >>>>>>>>>>>>>>>>>
//      
//      cout << "Number of external states: " << nbExternalStates << endl;
//      
//      // Loop over the external states that have not yet been assigned field values:
//      for (CFuint iState=0; iState<nbExternalStates; ++iState) {
//	const CFuint estateID = externalStates[iState]->getLocalID();
//	const CFuint ostateID = closestStates[iState]->getLocalID();
//	flag[estateID] = true;
//	
//	const CFreal Bx = ux[ostateID];
//	const CFreal By = uy[ostateID];
//	const CFreal Bz = uz[ostateID];
//	const CFreal normB = std::sqrt(Bx*Bx+By*By+Bz*Bz);
//	
//	// Query x,y,z of the current external state
//	const CFreal xe = (externalStates[iState]->getCoordinates())[0];
//	const CFreal ye = (externalStates[iState]->getCoordinates())[1];
//	const CFreal ze = (externalStates[iState]->getCoordinates())[2];
//	const CFreal re = externalStates[iState]->getCoordinates().norm2();
//	
//	// Query x,y,z from the closest state:
//	const CFreal x = (closestStates[iState]->getCoordinates())[0];
//	const CFreal y = (closestStates[iState]->getCoordinates())[1];
//	const CFreal z = (closestStates[iState]->getCoordinates())[2];
//	const CFreal r = closestStates[iState]->getCoordinates().norm2();
//	
//	// Compute theta and phi of the current external state from x,y,z
//	// theta = std::atan2(std::sqrt(x*x + y*y),z);
//	// phi = std::atan2(y,x);
//	
//	// Compute components Br, Btheta, Bphi at the source surface of ostate:
//	// Br = sin(theta)*cos(phi)*Bx + sin(theta)*sin(phi)*By + cos(theta)*Bz;
//	// Btheta = cos(theta)*cos(phi)*Bx + cos(theta)*sin(phi)*By - sin(theta)*Bz;
//	// Bphi = -sin(phi)*Bx + cos(phi)*By;
//	const CFreal Br = (x*Bx + y*By + z*Bz)/r;
//	
//	// Check whether Btheta and Bphi are small (ideally exactly zero)
//	// cout << "Btheta: " << Btheta << endl;
//	// cout << "Bphi: " << Bphi << endl;
//	
//	// Do the extrapolation:
//	const CFreal dist = closestStateDist[iState];
//	// cout << "dist: " << dist << endl;
//	const CFreal Bre = Br/(dist*dist);
//      
//	// From extrapolated Bre compute Bxe, Bye, Bze assuming a purely radial field
//	// Bxe = sin(theta)*cos(phi)*Bre;
//	// Bye = sin(theta)*sin(phi)*Bre;
//	// Bze = cos(theta)*Bre;
//	const CFreal Bxe = (xe/re)*Bre;
//	const CFreal Bye = (ye/re)*Bre;
//	const CFreal Bze = (ze/re)*Bre;
//	
//	// Write the extrapolated cartesian field values to the extended mesh:
//	// Just for checking where the patching error comes from: write 1.0 to states
//	(*states[estateID])[xVar] = Bxe;
//	(*states[estateID])[yVar] = Bye;
//	if (dim == DIM_3D) {
//	  (*states[estateID])[zVar] = Bze;
//	}   
//      }
//      
//      cf_assert(nbInternalStates + nbExternalStates == states.size());
//      
//      for (CFuint i = 0; i < nbStates; ++i) {
//	if (!flag[i]) {
//	  CFLog(INFO, "TransferField::execute() => stateID [" << i<< "] not found during extrapolation\n");
//	  cf_assert(flag[i]);
//	}
//      }
//      
//      CFLog(INFO,"TransferField::execute() => extrapolation took " << stp.read() << "s\n");
        }

        // AL: this has to be changed when doing unsteady cases!!!!!!
        // We need to implement a generic logic to apply the processing only when Poisson is solved
        m_applyProcessing = false;
        if (SubSystemStatusStack::getActive()->getNbIter()%this->getProcessRate()==0) 
        {
          m_applyProcessing = true;
        }
      } 
      
      //release the GeometricEntity
      m_cellBuilder->releaseGE();
    }
  }
  
  CFLog(INFO, "TransferField::execute() => END\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void TransferField::computeFlxPntStates(std::vector< RealVector > states, std::vector< RealVector >& statesFlxPnt)
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

void TransferField::setup()
{
  CFAUTOTRACE;
  
  FluxReconstructionSolverCom::setup();

  // get number of equations
  m_nbrEqs = PhysicalModelStack::getActive()->getNbEq();
  m_dim    = PhysicalModelStack::getActive()->getDim ();

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
  
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();

  // get the maximum number of flux points
  m_maxNbrFlxPnts = 0;
  
  // get the elementTypeData
  SafePtr< vector<ElementTypeData> > elemType = MeshDataStack::getActive()->getElementTypeData();
  
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    const CFuint nbrFlxPnts = frLocalData[iElemType]->getNbrOfFlxPnts();
    m_maxNbrFlxPnts = m_maxNbrFlxPnts > nbrFlxPnts ? m_maxNbrFlxPnts : nbrFlxPnts;
  }
  
  // get the coefs for extrapolation of the states to the flx pnts
  m_solPolyValsAtFlxPnts = frLocalData[0]->getCoefSolPolyInFlxPnts();
  
  for (CFuint iFlx = 0; iFlx < m_maxNbrFlxPnts; ++iFlx)
  {
    RealVector temp(m_nbrEqs);
    m_cellStatesFlxPnt.push_back(temp);
  }
  
  // sets the var IDs to a default in case it is not set in the CFcase
  if (m_variableIDs.size() == 0) 
  {
    m_variableIDs.resize(m_dim);
    for (CFuint i = 0; i < m_dim; ++i) 
    {
      m_variableIDs[i] = i;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void TransferField::unsetup()
{
  CFAUTOTRACE;

  //DataHandle< RealVector > nodeNghbCellMinAvgStates = socket_nodeNghbCellMinAvgStates.getDataHandle();
  //DataHandle< RealVector > nodeNghbCellMaxAvgStates = socket_nodeNghbCellMaxAvgStates.getDataHandle();
  //nodeNghbCellMinAvgStates.resize(0);
  //nodeNghbCellMaxAvgStates.resize(0);
  
  FluxReconstructionSolverCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
TransferField::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_otherStates);
  result.push_back(&socket_states);
  result.push_back(&socket_pastStates);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSource> >
TransferField::providesSockets()
{
  vector<SafePtr<BaseDataSocketSource> > result;

  //result.push_back(&socket_nodeNghbCellMinAvgStates);
  //result.push_back(&socket_nodeNghbCellMaxAvgStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
