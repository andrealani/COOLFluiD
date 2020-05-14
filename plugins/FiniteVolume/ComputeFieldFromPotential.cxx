#include "Common/PE.hh"
#include "Common/BadValueException.hh"
#include "Common/Stopwatch.hh"

#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/PhysicalModel.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/ComputeFieldFromPotential.hh"

#include <cmath>

/////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeFieldFromPotential,
		      CellCenterFVMData, FiniteVolumeModule>
ComputeFieldFromPotentialProvider("ComputeFieldFromPotential");

//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotential::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFuint> >
    ("VariableIDs", "IDs of the variable to be assigned to the newly computed field.");
  options.addConfigOption< string >
    ("OtherNamespace", "Name of the other namespace (providing the potential).");
  options.addConfigOption< CFreal >
    ("InterRadius",
     "Radius corresponding to the internal boundary between donor and current grids (<= 0 assumes one mesh).");
  options.addConfigOption< CFreal >
    ("DeltaSelection",
     "Distance within which points in the smaller mesh are selected.");
}
      
//////////////////////////////////////////////////////////////////////////////

ComputeFieldFromPotential::ComputeFieldFromPotential(const std::string& name) :
  CellCenterFVMCom(name),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_pastStates("pastStates"),
  socket_nodes("nodes"),
  socket_otherUX("uX"),
  socket_otherUY("uY"),
  socket_otherUZ("uZ"),
  socket_otherStates("states"),
  m_applyProcessing(true)
{
  addConfigOptionsTo(this);
  
  m_variableIDs = vector<CFuint>();
  setParameter("VariableIDs",&m_variableIDs);
  
  m_otherNamespace = "";
  setParameter("OtherNamespace", &m_otherNamespace);
  
  m_interRadius = -1.; // 2.5;
  setParameter("InterRadius", &m_interRadius);

  m_deltaSelection = 0.;
  setParameter("DeltaSelection", &m_deltaSelection);
}

//////////////////////////////////////////////////////////////////////////////

ComputeFieldFromPotential::~ComputeFieldFromPotential()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ComputeFieldFromPotential::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_nodes);
  result.push_back(&socket_otherUX);
  result.push_back(&socket_otherUY);
  result.push_back(&socket_otherUZ);
  result.push_back(&socket_otherStates);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotential::setup()
{
  CFAUTOTRACE;

  CellCenterFVMCom::setup();

  if (m_variableIDs.size() == 0) {
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    m_variableIDs.resize(dim);
    for (CFuint i = 0; i < dim; ++i) {
      m_variableIDs[i] = i;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotential::configure ( Config::ConfigArgs& args )
{
  CellCenterFVMCom::configure(args);

  cf_assert(m_otherNamespace != "");
  CFLog(VERBOSE, "ComputeFieldFromPotential::configure() => m_otherNamespace = " <<
	m_otherNamespace << "\n");
  socket_otherUX.setDataSocketNamespace(m_otherNamespace);
  socket_otherUY.setDataSocketNamespace(m_otherNamespace);
  socket_otherUZ.setDataSocketNamespace(m_otherNamespace);
  socket_otherStates.setDataSocketNamespace(m_otherNamespace);
}
      
//////////////////////////////////////////////////////////////////////////////
 
void ComputeFieldFromPotential::execute()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "ComputeFieldFromPotential::execute() => START\n");

  if (SubSystemStatusStack::getActive()->getNbIter() >= 1 && m_applyProcessing) {
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    
    Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
      (SubSystemStatusStack::getCurrentName()).getNamespace(m_otherNamespace);
    Common::SafePtr<SubSystemStatus> otherSubSystemStatus =
      SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);

    DataHandle<CFreal> ux = socket_otherUX.getDataHandle();
    DataHandle<CFreal> uy = socket_otherUY.getDataHandle();
    DataHandle<CFreal> uz = socket_otherUZ.getDataHandle();
    
    // States on the extended corona MHD mesh
    DataHandle<State*, GLOBAL> states = socket_states.getDataHandle(); // LARGER mesh
    
    const CFuint nbStates = states.size();
    cf_assert(dim >= DIM_2D);
    cf_assert(m_variableIDs.size() >= 2);
    const CFuint xVar = m_variableIDs[0];
    cf_assert(xVar < nbEqs);
    const CFuint yVar = m_variableIDs[1];
    cf_assert(yVar < nbEqs);
    const CFuint zVar = (dim == DIM_3D) ? m_variableIDs[2] : 0;
    cf_assert(zVar < nbEqs);
    
    Stopwatch<WallTime> stp;
    stp.start();
    
    if (m_interRadius <= 0.) {
      CFLog(INFO, "ComputeFieldFromPotential::execute() => transferring field\n");
      
      for (CFuint iState = 0; iState < nbStates; ++iState) {
	(*states[iState])[xVar] = ux[iState];
	(*states[iState])[yVar] = uy[iState];
	if (dim == DIM_3D) {
	  (*states[iState])[zVar] = uz[iState];
	}
      }
    }
    else {
      DataHandle<State*, GLOBAL> otherStates = socket_otherStates.getDataHandle(); // SMALLER mesh
      
      // loop over internal faces
      //   if face has all its nodes radius < (m_interRadius + eps) and > (m_interRadius - eps)
      //     store the face states (or stateIDs) whose radius < m_interRadius as interStates
      //     those states will be those within which looking for the matching "fronteer" otherStates 
      vector<State*> interStates;
      interStates.reserve(otherStates.size()/4); // rough estimation
      
      SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");
      SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> > faceBuilder = 
	this->getMethodData().getFaceTrsGeoBuilder();
      FaceTrsGeoBuilder::GeoData& geoData = faceBuilder->getDataGE();
      SafePtr<FaceTrsGeoBuilder> faceBuilderPtr = faceBuilder->getGeoBuilder();
      faceBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
      geoData.trs = faces;
      geoData.isBFace = false;

      cf_assert(m_deltaSelection > 0.);
      const CFreal rMax = m_interRadius + m_deltaSelection;
      const CFreal rMin = m_interRadius - m_deltaSelection;
      
      const CFuint nbFaces = faces->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	geoData.idx = iFace;
	const GeometricEntity *const face = faceBuilder->buildGE();
	const vector<Node*>& nodesInFace = face->getNodes();
	const CFuint nbNodesInFace = nodesInFace.size();
	cf_assert(nbNodesInFace > 1);
	CFuint countInterNodes = 0;
	CFreal radius = 0.;
	for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
	  radius = nodesInFace[iNode]->norm2();
	  if (radius < rMax && radius > rMin) {
	    countInterNodes++;
	  }
	}
	if (countInterNodes == nbNodesInFace) {
	  // store state on the inner side of the internal surface
	  const CFreal radius0 = face->getState(0)->getCoordinates().norm2();
	  const CFreal radius1 = face->getState(1)->getCoordinates().norm2();
	  (radius0 < radius1) ? interStates.push_back(face->getState(0)) : interStates.push_back(face->getState(1)); 
	  CFLog(DEBUG_MIN, "ComputeFieldFromPotential::execute() => #" << interStates.size() <<  " face with radius [" << radius << "] detected\n");
	}
	faceBuilder->releaseGE();
      }

      const CFreal maxFloat = std::numeric_limits<double>::max();
      
      // find the closest otherStates for each of the selected interStates
      const CFuint nbOtherStates = otherStates.size();
      const CFuint nbInterStates = interStates.size();
      vector<State*> interOtherStates(nbInterStates); // internal states (smaller mesh) attached to m_interRadius boundary
      vector<CFreal> distMin(nbInterStates, maxFloat);
      for (CFuint iState = 0; iState < nbInterStates; ++iState) {
	const Node& coord = interStates[iState]->getCoordinates();
	for (CFuint jState = 0; jState < nbOtherStates; ++jState) {
	  const CFreal dist2 = MathFunctions::getSquaredDistance(coord, otherStates[jState]->getCoordinates());
	  if (dist2 < distMin[iState]) {
	    distMin[iState] = dist2;
	    interOtherStates[iState] = otherStates[jState];
	  }
	}
	
	CFLog(DEBUG_MIN, "ComputeFieldFromPotential::execute() => distMin[" << iState << "] = " << std::sqrt(distMin[iState]) << "\n");
      }
      
      // count the external states and store the corresponding state IDs in a list
      vector<State*> externalStates;
      externalStates.reserve(nbStates); // size overestimated to avoid reallocation
      vector<State*> internalStates;
      internalStates.reserve(nbStates); // size overestimated to avoid reallocation
      
      const CFreal radiusPFSS = m_interRadius + m_deltaSelection;
      for (CFuint iState = 0; iState < nbStates; ++iState) {
	if (states[iState]->getCoordinates().norm2() > radiusPFSS) {
	  cf_assert(iState == states[iState]->getLocalID());
	  externalStates.push_back(states[iState]);
	}
	else {
	  internalStates.push_back(states[iState]);
	}
      }
      cf_assert(externalStates.size() + internalStates.size() == nbStates);
      
      CFLog(INFO, "ComputeFieldFromPotential::execute() => detected [" << externalStates.size() << "] external states\n");
      
      // loop over all otherStates
      //   if state radius > m_interRadius
      //     loop over all interStates to find the closest state and store the minimal distance
      //   else increment internal state counter (which will have to match otherStates.size() at the end)
      const CFuint nbExternalStates = externalStates.size();

      vector<State*> closestStates(nbExternalStates);
      RealVector closestStateDist(maxFloat, nbExternalStates);
      
      for (CFuint iState = 0; iState < nbExternalStates; ++iState) {
	const Node& coord = externalStates[iState]->getCoordinates();
	for (CFuint jState = 0; jState < nbInterStates; ++jState) {
	  State* const intState = interOtherStates[jState];
	  // use square distance for comparison since it is faster to compute 
	  const CFreal dist2 =
	    MathFunctions::getSquaredDistance(coord, intState->getCoordinates());
	  if (dist2 < closestStateDist[iState]) {
	    closestStateDist[iState] = dist2;
	    closestStates[iState] = intState;
	  }
	}
      }            
      closestStateDist = sqrt(closestStateDist);
      
      const CFuint nbInternalStates = internalStates.size();
      vector<State*> closestInternalStates(nbInternalStates);

      for (CFuint iState = 0; iState < nbInternalStates; ++iState) {
	const Node& coord = internalStates[iState]->getCoordinates();
	CFreal distMin = maxFloat;
	for (CFuint jState = 0; jState < nbOtherStates; ++jState) {
	  State* const intState = otherStates[jState];
	  // use square distance for comparison since it is faster to compute 
	  const CFreal dist2 =
	    MathFunctions::getSquaredDistance(coord, intState->getCoordinates());
	  if (dist2 < distMin) {
	    distMin = dist2;
	    closestInternalStates[iState] = intState;
	  }
	}
      } 
      
      // from here on:
      // 1) externalStates[k]   gives the State from LARGER mesh with radius > m_interRadius 
      // 2) closestStates[k]    gives the State from SMALLER mesh to be used for extrapolation
      // 3) closestStateDist[k] gives the distance between 1) and 2)
      // 4) closestInternalStates[k]->getLocalID() gives the state ID to be used to fetch the gradient in uX, uY, uZ
      CFLog(VERBOSE, "size internalStates => " << internalStates[0]->size() << "\n");
      CFLog(VERBOSE, "size closestInternalStates  => " << closestInternalStates[0]->size() << "\n");

      // array to keep track of the state IDs that correspond to internal and external states
      vector<bool> flag(nbStates, false);
      
      // for (CFuint iState = 0; iState < nbExternalStates; ++iState) {
      //	const CFuint estateID = externalStates[iState]->getLocalID();
      //	const CFuint ostateID = closestStates[iState]->getLocalID();
      /// use ux_i[ostateID] to compute Br, Btheta, Bphi
      //	flag[estateID] = true;
      //	(*states[estateID])[xVar] = ux[ostateID];
      //	(*states[estateID])[yVar] = uy[ostateID];
      //	if (dim == DIM_3D) {
      //	  (*states[estateID])[zVar] = uz[ostateID];
      //	}
      // }

      // from here on:
      // 1) internalStates[k]   gives the State from LARGER mesh with radius < m_interRadius 
      // 2) closestStates[k]    gives the State from SMALLER mesh to be used for extrapolation
      // 3) closestStateDist[k] gives the distance between 1) and 2)
      // 4) closestStates[k]->getLocalID() gives the state ID to be used to fetch the gradient in uX, uY, uZ
            
      CFLog(VERBOSE, "size internalStates => " << internalStates[0]->size() << "\n");
      CFLog(VERBOSE, "size closestStates  => " << closestStates[0]->size() << "\n");
      CFLog(VERBOSE, "size internalStates => " << internalStates[0]->size() << "\n");
      CFLog(VERBOSE, "size closestInternalStates  => " << closestInternalStates[0]->size() << "\n");
      
      for (CFuint iState = 0; iState < nbInternalStates; ++iState) {
	const CFuint istateID = internalStates[iState]->getLocalID();
	const CFuint ostateID = closestInternalStates[iState]->getLocalID();
	flag[istateID] = true;
	/// use ux_i[ostateID] to compute Br, Btheta, Bphi
	(*states[istateID])[xVar] = ux[ostateID];
	(*states[istateID])[yVar] = uy[ostateID];
	if (dim == DIM_3D) {
	  (*states[istateID])[zVar] = uz[ostateID];
	}
      }
      
      //<<<<<<<<<<< HERE FOLLOWS THE EXTRAPOLATION ONTO THE OUTER MESH >>>>>>>>>>>>>>>>>
      
      cout << "Number of external states: " << nbExternalStates << endl;
      
      // Loop over the external states that have not yet been assigned field values:
      for (CFuint iState=0; iState<nbExternalStates; ++iState) {
	const CFuint estateID = externalStates[iState]->getLocalID();
	const CFuint ostateID = closestStates[iState]->getLocalID();
	flag[estateID] = true;
	
	const CFreal Bx = ux[ostateID];
	const CFreal By = uy[ostateID];
	const CFreal Bz = uz[ostateID];
	const CFreal normB = std::sqrt(Bx*Bx+By*By+Bz*Bz);
	
	// Query x,y,z of the current external state
	const CFreal xe = (externalStates[iState]->getCoordinates())[0];
	const CFreal ye = (externalStates[iState]->getCoordinates())[1];
	const CFreal ze = (externalStates[iState]->getCoordinates())[2];
	const CFreal re = externalStates[iState]->getCoordinates().norm2();
	
	// Query x,y,z from the closest state:
	const CFreal x = (closestStates[iState]->getCoordinates())[0];
	const CFreal y = (closestStates[iState]->getCoordinates())[1];
	const CFreal z = (closestStates[iState]->getCoordinates())[2];
	const CFreal r = closestStates[iState]->getCoordinates().norm2();
	
	// Compute theta and phi of the current external state from x,y,z
	// theta = std::atan2(std::sqrt(x*x + y*y),z);
	// phi = std::atan2(y,x);
	
	// Compute components Br, Btheta, Bphi at the source surface of ostate:
	// Br = sin(theta)*cos(phi)*Bx + sin(theta)*sin(phi)*By + cos(theta)*Bz;
	// Btheta = cos(theta)*cos(phi)*Bx + cos(theta)*sin(phi)*By - sin(theta)*Bz;
	// Bphi = -sin(phi)*Bx + cos(phi)*By;
	const CFreal Br = (x*Bx + y*By + z*Bz)/r;
	
	// Check whether Btheta and Bphi are small (ideally exactly zero)
	// cout << "Btheta: " << Btheta << endl;
	// cout << "Bphi: " << Bphi << endl;
	
	// Do the extrapolation:
	const CFreal dist = closestStateDist[iState];
	// cout << "dist: " << dist << endl;
	const CFreal Bre = Br/(dist*dist);
      
	// From extrapolated Bre compute Bxe, Bye, Bze assuming a purely radial field
	// Bxe = sin(theta)*cos(phi)*Bre;
	// Bye = sin(theta)*sin(phi)*Bre;
	// Bze = cos(theta)*Bre;
	const CFreal Bxe = (xe/re)*Bre;
	const CFreal Bye = (ye/re)*Bre;
	const CFreal Bze = (ze/re)*Bre;
	
	// Write the extrapolated cartesian field values to the extended mesh:
	// Just for checking where the patching error comes from: write 1.0 to states
	(*states[estateID])[xVar] = Bxe;
	(*states[estateID])[yVar] = Bye;
	if (dim == DIM_3D) {
	  (*states[estateID])[zVar] = Bze;
	}   
      }
      
      cf_assert(nbInternalStates + nbExternalStates == states.size());
      
      for (CFuint i = 0; i < nbStates; ++i) {
	if (!flag[i]) {
	  CFLog(INFO, "ComputeFieldFromPotential::execute() => stateID [" << i<< "] not found during extrapolation\n");
	  cf_assert(flag[i]);
	}
      }
      
      CFLog(INFO,"ComputeFieldFromPotential::execute() => extrapolation took " << stp.read() << "s\n");
    }

    // AL: this has to be changed when doing unsteady cases!!!!!!
    // We need to implement a generic logic to apply the processing only when Poisson is solved
    m_applyProcessing = false;
    if (SubSystemStatusStack::getActive()->getNbIter()%this->getProcessRate()==0) {
      m_applyProcessing = true;
    }

    
    //CFLog(INFO, "State[20415] after B field =  " << (*(states[20415]))[4] << " " << (*(states[20415]))[5] << " " << (*(states[20415]))[6] << ", [x y z] = "<< states[20415]->getCoordinates() << "\n");
    
    // update the pastStates
    DataHandle<State*> pastStates  = socket_pastStates.getDataHandle();
    // Set Initial States to current states
    if(SubSystemStatusStack::getActive()->isSubIterationFirstStep()){
      for(CFuint i = 0; i < states.size(); ++i) {
	*(pastStates[i]) = *(states[i]);
      }
    }
  } 
  
  CFLog(VERBOSE, "ComputeFieldFromPotential::execute() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

 /*void ComputeFieldFromPotential::execute()
{
  CFAUTOTRACE;
  
  CFLog(INFO, "ComputeFieldFromPotential::execute() => START\n");

  if (SubSystemStatusStack::getActive()->getNbIter() >= 1) {
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    
    Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
      (SubSystemStatusStack::getCurrentName()).getNamespace(m_otherNamespace);
    Common::SafePtr<SubSystemStatus> otherSubSystemStatus =
      SubSystemStatusStack::getInstance().getEntryByNamespace(nsp);

    DataHandle<CFreal> ux = socket_otherUX.getDataHandle();
    DataHandle<CFreal> uy = socket_otherUY.getDataHandle();
    DataHandle<CFreal> uz = socket_otherUZ.getDataHandle();
    
    // States on the extended corona MHD mesh
    DataHandle<State*, GLOBAL> states = socket_states.getDataHandle(); // LARGER mesh
    
    const CFuint nbStates = states.size();
    cf_assert(dim >= DIM_2D);
    cf_assert(m_variableIDs.size() >= 2);
    const CFuint xVar = m_variableIDs[0];
    cf_assert(xVar < nbEqs);
    const CFuint yVar = m_variableIDs[1];
    cf_assert(yVar < nbEqs);
    const CFuint zVar = (dim == DIM_3D) ? m_variableIDs[2] : 0;
    cf_assert(zVar < nbEqs);
    
    Stopwatch<WallTime> stp;
    stp.start();
    
    if (m_interRadius <= 0.) {
      for (CFuint iState = 0; iState < nbStates; ++iState) {
	(*states[iState])[xVar] = ux[iState];
	(*states[iState])[yVar] = uy[iState];
	if (dim == DIM_3D) {
	  (*states[iState])[zVar] = uz[iState];
	}
      }
    }
    else {
      DataHandle<State*, GLOBAL> otherStates = socket_otherStates.getDataHandle(); // SMALLER mesh
      
      // loop over internal faces
      //   if face has all its nodes radius < (m_interRadius + eps) and > (m_interRadius - eps)
      //     store the face states (or stateIDs) whose radius < m_interRadius as interStates
      //     those states will be those within which looking for the matching "fronteer" otherStates 
      vector<State*> interStates;
      interStates.reserve(otherStates.size()/4); // rough estimation
      
      SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");
      SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> > faceBuilder = 
	this->getMethodData().getFaceTrsGeoBuilder();
      FaceTrsGeoBuilder::GeoData& geoData = faceBuilder->getDataGE();
      SafePtr<FaceTrsGeoBuilder> faceBuilderPtr = faceBuilder->getGeoBuilder();
      faceBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
      geoData.trs = faces;
      geoData.isBFace = false;

      cf_assert(m_deltaSelection > 0.);
      const CFreal rMax = m_interRadius + m_deltaSelection;
      const CFreal rMin = m_interRadius - m_deltaSelection;
      
      const CFuint nbFaces = faces->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
	geoData.idx = iFace;
	const GeometricEntity *const face = faceBuilder->buildGE();
	const vector<Node*>& nodesInFace = face->getNodes();
	const CFuint nbNodesInFace = nodesInFace.size();
	cf_assert(nbNodesInFace > 1);
	CFuint countInterNodes = 0;
	CFreal radius = 0.;
	for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
	  radius = nodesInFace[iNode]->norm2();
	  if (radius < rMax && radius > rMin) {
	    countInterNodes++;
	  }
	}
	if (countInterNodes == nbNodesInFace) {
	  // store state on the inner side of the internal surface
	  const CFreal radius0 = face->getState(0)->getCoordinates().norm2();
	  const CFreal radius1 = face->getState(1)->getCoordinates().norm2();
	  (radius0 < radius1) ? interStates.push_back(face->getState(0)) : interStates.push_back(face->getState(1)); 
	  CFLog(DEBUG_MIN, "ComputeFieldFromPotential::execute() => #" << interStates.size() <<  " face with radius [" << radius << "] detected\n");
	}
	faceBuilder->releaseGE();
      }

      const CFreal maxFloat = std::numeric_limits<double>::max();
      
      // find the closest GLOBAL otherStates for each of the selected LOCAL interStates
      const CFuint nbOtherStates = otherStates.size();
      const CFuint nbInterStates = interStates.size();
      // internal states (smaller mesh) attached to m_interRadius boundary
      vector<State*> interOtherStates(nbInterStates); 
      vector<CFreal> distMin(nbInterStates, maxFloat);
      
      
      MPIStruct ms;
      CFuint sizeOtherStates = nbOtherStates*otherStates[0]->size();
      CFuint sizeOtherCoords = nbOtherStates*dim;
      
      for (CFuint root = 0; root < m_nbProc; ++root) {
	if (m_nbProc > 1) {
	  int ln[2];
	  ln[0] = ln[1] = 1;
	  
	  MPIStructDef::buildMPIStruct(&sizeOtherStates, &sizeOtherCoords, ln, ms);
	  MPI_Bcast(ms.start, 1, ms.type, root, m_comm);
	} 
	
	vector<CFreal> otherStatesArray;
	vector<CFreal> otherCoordsArray;

	// current rank copies its own other states and coordinates
	if (m_myRank == root) {
	  otherStatesArray.resize(sizeOtherStates);
	  otherCoordsArray.resize(sizeOtherCoords);
	  for (CFuint s = 0; s < nbOtherStates; ++s) {
	    const CFuint starts = s*nbOtherEqs;
	    for (CFuint iEq = 0; iEq < nbOtherEqs; ++iEq) {
	      otherStatesArray[start+iEq] = (*otherStates[s])[iEq]; 
	    }
	    for (CFuint iEq = 0; iEq < dim; ++iEq) {
	      otherCoordsArray[start+iEq] = (otherStates[s]->getCoordinates())[iEq]; 
	    }
	  }
	}

	if (m_nbProc > 1) {
	  int le[2];
	  // AL: to be tested, here we are communicating just a part of the arrays
	  le[0] = sizeOtherStates;
	  le[1] = sizeOtherCoords;
	  
	  MPIStructDef::buildMPIStruct<CFuint, CFuint>
	    (&otherStatesArray[0], &otherCoordsArray[0], le, ms);
	  MPI_Bcast(ms.start, 1, ms.type, root, m_comm);
	}
	
	for (CFuint iState = 0; iState < nbInterStates; ++iState) {
	const Node& coord = interStates[iState]->getCoordinates();
	for (CFuint jState = 0; jState < nbOtherStates; ++jState) {
	  const CFreal dist2 =
	    MathFunctions::getSquaredDistance(coord, otherStates[jState]->getCoordinates());
	  if (dist2 < distMin[iState]) {
	    distMin[iState] = dist2;
	    interOtherStates[iState] = otherStates[jState];
	  }
	}
	
	CFLog(DEBUG_MIN, "ComputeFieldFromPotential::execute() => distMin[" << iState << "] = " << std::sqrt(distMin[iState]) << "\n");
      }
      








      
      
      // count the external states and store the corresponding state IDs in a list
      vector<State*> externalStates;
      externalStates.reserve(nbStates); // size overestimated to avoid reallocation
      vector<State*> internalStates;
      internalStates.reserve(nbStates); // size overestimated to avoid reallocation
      const CFreal radiusPFSS = m_interRadius + m_deltaSelection;
      for (CFuint iState = 0; iState < nbStates; ++iState) 
      if (states[iState]->getCoordinates().norm2() > radiusPFSS) {
      cf_assert(iState == states[iState]->getLocalID());
      externalStates.push_back(states[iState]);
	}
	else {
	  internalStates.push_back(states[iState]);
	}
      }
      cf_assert(externalStates.size() + internalStates.size() == nbStates);
      
      CFLog(INFO, "ComputeFieldFromPotential::execute() => detected [" << externalStates.size() << "] external states\n");
      
      // loop over all otherStates
      //   if state radius > m_interRadius
      //     loop over all interStates to find the closest state and store the minimal distance
      //   else increment internal state counter (which will have to match otherStates.size() at the end)
      const CFuint nbExternalStates = externalStates.size();

      vector<State*> closestStates(nbExternalStates);
      RealVector closestStateDist(maxFloat, nbExternalStates);
      
      for (CFuint iState = 0; iState < nbExternalStates; ++iState) {
	const Node& coord = externalStates[iState]->getCoordinates();
	for (CFuint jState = 0; jState < nbInterStates; ++jState) {
	  State* const intState = interOtherStates[jState];
	  // use square distance for comparison since it is faster to compute 
	  const CFreal dist2 =
	    MathFunctions::getSquaredDistance(coord, intState->getCoordinates());
	  if (dist2 < closestStateDist[iState]) {
	    closestStateDist[iState] = dist2;
	    closestStates[iState] = intState;
	  }
	}
      }            
      closestStateDist = sqrt(closestStateDist);
      
      const CFuint nbInternalStates = internalStates.size();
      vector<State*> closestInternalStates(nbInternalStates);

      for (CFuint iState = 0; iState < nbInternalStates; ++iState) {
	const Node& coord = internalStates[iState]->getCoordinates();
	CFreal distMin = maxFloat;
	for (CFuint jState = 0; jState < nbOtherStates; ++jState) {
	  State* const intState = otherStates[jState];
	  // use square distance for comparison since it is faster to compute 
	  const CFreal dist2 =
	    MathFunctions::getSquaredDistance(coord, intState->getCoordinates());
	  if (dist2 < distMin) {
	    distMin = dist2;
	    closestInternalStates[iState] = intState;
	  }
	}
      } 
      
      // from here on:
      // 1) externalStates[k]   gives the State from LARGER mesh with radius > m_interRadius 
      // 2) closestStates[k]    gives the State from SMALLER mesh to be used for extrapolation
      // 3) closestStateDist[k] gives the distance between 1) and 2)
      // 4) closestInternalStates[k]->getLocalID() gives the state ID to be used to fetch the gradient in uX, uY, uZ
      CFLog(VERBOSE, "size internalStates => " << internalStates[0]->size() << "\n");
      CFLog(VERBOSE, "size closestInternalStates  => " << closestInternalStates[0]->size() << "\n");

      // array to keep track of the state IDs that correspond to internal and external states
      vector<bool> flag(nbStates, false);
      
      // for (CFuint iState = 0; iState < nbExternalStates; ++iState) {
      //	const CFuint estateID = externalStates[iState]->getLocalID();
      //	const CFuint ostateID = closestStates[iState]->getLocalID();
      /// use ux_i[ostateID] to compute Br, Btheta, Bphi
      //	flag[estateID] = true;
      //	(*states[estateID])[xVar] = ux[ostateID];
      //	(*states[estateID])[yVar] = uy[ostateID];
      //	if (dim == DIM_3D) {
      //	  (*states[estateID])[zVar] = uz[ostateID];
      //	}
      // }

      // from here on:
      // 1) internalStates[k]   gives the State from LARGER mesh with radius < m_interRadius 
      // 2) closestStates[k]    gives the State from SMALLER mesh to be used for extrapolation
      // 3) closestStateDist[k] gives the distance between 1) and 2)
      // 4) closestStates[k]->getLocalID() gives the state ID to be used to fetch the gradient in uX, uY, uZ
            
      CFLog(VERBOSE, "size internalStates => " << internalStates[0]->size() << "\n");
      CFLog(VERBOSE, "size closestStates  => " << closestStates[0]->size() << "\n");
      CFLog(VERBOSE, "size internalStates => " << internalStates[0]->size() << "\n");
      CFLog(VERBOSE, "size closestInternalStates  => " << closestInternalStates[0]->size() << "\n");
      
      for (CFuint iState = 0; iState < nbInternalStates; ++iState) {
	const CFuint istateID = internalStates[iState]->getLocalID();
	const CFuint ostateID = closestInternalStates[iState]->getLocalID();
	flag[istateID] = true;
	/// use ux_i[ostateID] to compute Br, Btheta, Bphi
	(*states[istateID])[xVar] = ux[ostateID];
	(*states[istateID])[yVar] = uy[ostateID];
	if (dim == DIM_3D) {
	  (*states[istateID])[zVar] = uz[ostateID];
	}
      }
      
      //<<<<<<<<<<< HERE FOLLOWS THE EXTRAPOLATION ONTO THE OUTER MESH >>>>>>>>>>>>>>>>>
      
      cout << "Number of external states: " << nbExternalStates << endl;
      
      // Loop over the external states that have not yet been assigned field values:
      for (CFuint iState=0; iState<nbExternalStates; ++iState) {
	const CFuint estateID = externalStates[iState]->getLocalID();
	const CFuint ostateID = closestStates[iState]->getLocalID();
	flag[estateID] = true;
	
	const CFreal Bx = ux[ostateID];
	const CFreal By = uy[ostateID];
	const CFreal Bz = uz[ostateID];
	const CFreal normB = std::sqrt(Bx*Bx+By*By+Bz*Bz);
	
	// Query x,y,z of the current external state
	const CFreal xe = (externalStates[iState]->getCoordinates())[0];
	const CFreal ye = (externalStates[iState]->getCoordinates())[1];
	const CFreal ze = (externalStates[iState]->getCoordinates())[2];
	const CFreal re = externalStates[iState]->getCoordinates().norm2();
	
	// Query x,y,z from the closest state:
	const CFreal x = (closestStates[iState]->getCoordinates())[0];
	const CFreal y = (closestStates[iState]->getCoordinates())[1];
	const CFreal z = (closestStates[iState]->getCoordinates())[2];
	const CFreal r = closestStates[iState]->getCoordinates().norm2();
	
	// Compute theta and phi of the current external state from x,y,z
	// theta = std::atan2(std::sqrt(x*x + y*y),z);
	// phi = std::atan2(y,x);
	
	// Compute components Br, Btheta, Bphi at the source surface of ostate:
	// Br = sin(theta)*cos(phi)*Bx + sin(theta)*sin(phi)*By + cos(theta)*Bz;
	// Btheta = cos(theta)*cos(phi)*Bx + cos(theta)*sin(phi)*By - sin(theta)*Bz;
	// Bphi = -sin(phi)*Bx + cos(phi)*By;
	const CFreal Br = (x*Bx + y*By + z*Bz)/r;
	
	// Check whether Btheta and Bphi are small (ideally exactly zero)
	// cout << "Btheta: " << Btheta << endl;
	// cout << "Bphi: " << Bphi << endl;
	
	// Do the extrapolation:
	const CFreal dist = closestStateDist[iState];
	// cout << "dist: " << dist << endl;
	const CFreal Bre = Br/(dist*dist);
      
	// From extrapolated Bre compute Bxe, Bye, Bze assuming a purely radial field
	// Bxe = sin(theta)*cos(phi)*Bre;
	// Bye = sin(theta)*sin(phi)*Bre;
	// Bze = cos(theta)*Bre;
	const CFreal Bxe = (xe/re)*Bre;
	const CFreal Bye = (ye/re)*Bre;
	const CFreal Bze = (ze/re)*Bre;
	
	// Write the extrapolated cartesian field values to the extended mesh:
	// Just for checking where the patching error comes from: write 1.0 to states
	(*states[estateID])[xVar] = Bxe;
	(*states[estateID])[yVar] = Bye;
	if (dim == DIM_3D) {
	  (*states[estateID])[zVar] = Bze;
	}   
      }
      
      cf_assert(nbInternalStates + nbExternalStates == states.size());
      
      for (CFuint i = 0; i < nbStates; ++i) {
	if (!flag[i]) {
	  CFLog(INFO, "ComputeFieldFromPotential::execute() => stateID [" << i<< "] not found during extrapolation\n");
	  cf_assert(flag[i]);
	}
      }
      
      CFLog(INFO,"ComputeFieldFromPotential::execute() => extrapolation took " << stp.read() << "s\n");
      /// missing states inside the m_interRadius
    }
  } 
  
  CFLog(INFO, "ComputeFieldFromPotential::execute() => END\n");
}
 */
//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotential::unsetup()
{
  CFAUTOTRACE;
  
  CellCenterFVMCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

