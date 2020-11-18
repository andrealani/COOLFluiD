#include "Common/PE.hh"
#include "Common/BadValueException.hh"
#include "Common/Stopwatch.hh"

#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/PhysicalModel.hh"

#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "FiniteVolumeMHD/ComputeFieldFromPotentialMHD.hh"

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

MethodCommandProvider<ComputeFieldFromPotentialMHD,
		      CellCenterFVMData, FiniteVolumeMHDModule>
ComputeFieldFromPotentialMHDProvider("ComputeFieldFromPotentialMHD");

//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotentialMHD::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFuint> >
    ("VariableIDs", "IDs of the variable to be assigned to the newly computed field.");
  options.addConfigOption< string >
    ("OtherNamespace", "Name of the other namespace (providing the potential).");
}
      
//////////////////////////////////////////////////////////////////////////////

ComputeFieldFromPotentialMHD::ComputeFieldFromPotentialMHD(const std::string& name) :
  CellCenterFVMCom(name),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_pastStates("pastStates"),
  socket_nodes("nodes"),
  socket_otherUX("uX"),
  socket_otherUY("uY"),
  socket_otherUZ("uZ"),
  socket_faceCenters("faceCenters"),
  socket_otherStates("states"),
  socket_BPFSS("BPFSS"),
  socket_BPFSSCells("BPFSSCells"),
  m_applyProcessing(true)
{
  addConfigOptionsTo(this);
  
  m_variableIDs = vector<CFuint>();
  setParameter("VariableIDs",&m_variableIDs);
  
  m_otherNamespace = "";
  setParameter("OtherNamespace", &m_otherNamespace);
}

//////////////////////////////////////////////////////////////////////////////

ComputeFieldFromPotentialMHD::~ComputeFieldFromPotentialMHD()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ComputeFieldFromPotentialMHD::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_BPFSS);
  result.push_back(&socket_BPFSSCells);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ComputeFieldFromPotentialMHD::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_nodes);
  result.push_back(&socket_otherUX);
  result.push_back(&socket_otherUY);
  result.push_back(&socket_otherUZ);
  result.push_back(&socket_faceCenters);
  result.push_back(&socket_otherStates);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotentialMHD::setup()
{
  CFAUTOTRACE;

  CellCenterFVMCom::setup();

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  if (m_variableIDs.size() == 0) {
    m_variableIDs.resize(dim);
    for (CFuint i = 0; i < dim; ++i) {
      m_variableIDs[i] = i;
    }
  }
  assert(m_variableIDs.size() == dim);

  const CFuint nbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  DataHandle<std::vector<CFreal> > BPFSS = socket_BPFSS.getDataHandle();
  BPFSS.resize(nbFaces);
  for (CFuint i = 0; i < nbFaces; ++i) {
    BPFSS[i].resize(dim);
  }
  
  DataHandle<State*, GLOBAL> states = socket_states.getDataHandle();
  const CFuint nbStates = states.size();
  DataHandle<CFreal> cellBPFSS = socket_BPFSSCells.getDataHandle();
  cellBPFSS.resize(nbStates*dim);
}
      
//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotentialMHD::configure ( Config::ConfigArgs& args )
{
  CellCenterFVMCom::configure(args);

  cf_assert(m_otherNamespace != "");
  CFLog(VERBOSE, "ComputeFieldFromPotentialMHD::configure() => m_otherNamespace = " <<
	m_otherNamespace << "\n");
  socket_otherUX.setDataSocketNamespace(m_otherNamespace);
  socket_otherUY.setDataSocketNamespace(m_otherNamespace);
  socket_otherUZ.setDataSocketNamespace(m_otherNamespace);
  socket_otherStates.setDataSocketNamespace(m_otherNamespace);
}
      
//////////////////////////////////////////////////////////////////////////////
 
void ComputeFieldFromPotentialMHD::execute()
{
  CFAUTOTRACE;

  CFLog(VERBOSE, "ComputeFieldFromPotentialMHD::execute() => START\n");

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
    
    DataHandle<State*, GLOBAL> states = socket_states.getDataHandle(); // MHD states
    DataHandle<CFreal> cellBPFSS = socket_BPFSSCells.getDataHandle();
    
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
       
    CFLog(INFO, "ComputeFieldFromPotentialMHD::execute() => transferring field\n");
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      const CFuint startID = iState*dim;
      cellBPFSS[startID]   = (*states[iState])[xVar] = ux[iState];
      cellBPFSS[startID+1] = (*states[iState])[yVar] = uy[iState];
      if (dim == DIM_3D) {
	cellBPFSS[startID+2] = (*states[iState])[zVar] = uz[iState];
      }
    }
    
    DataHandle<CFreal> faceCenters = socket_faceCenters.getDataHandle();
    DataHandle<std::vector<CFreal> > BPFSSFace = socket_BPFSS.getDataHandle();
    
    // prepare the building of the faces
    Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> > geoBuilder =
      getMethodData().getFaceTrsGeoBuilder();
    geoBuilder->getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
    FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
    
    // get all the TRSs
    vector<SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
    const CFuint nbTRSs = trs.size();
    
    // temporary wrappers for pointers
    RealVector xcFacePtr(dim, static_cast<CFreal*>(CFNULL));
    
    for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS) {
      SafePtr<TopologicalRegionSet> currTrs = trs[iTRS];
      
      bool isBFace = false;
      CFLog(VERBOSE, "TRS name = " << currTrs->getName() << "\n");
      if (currTrs->getName() != "PartitionFaces" && currTrs->getName() != "InnerCells") {
	if (currTrs->hasTag("writable")) {
	  isBFace = geoData.isBFace = true;
	}
	else {
	  isBFace = geoData.isBFace = false;
	}
	
	// set the current TRS in the geoData
	geoData.trs = currTrs;
	
	const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
	for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
	  // build the GeometricEntity
	  geoData.idx = iFace;
	  GeometricEntity *const face = geoBuilder->buildGE();
	  const vector<Node*>& nodes = *face->getNodes();
	  const CFuint nbNodesInFace = nodes.size();
	  const CFuint faceID  = face->getID();
	  const CFuint startID = faceID*dim;
	  xcFacePtr.wrap(dim, &faceCenters[startID]);
	  
	  State* const stateL = face->getState(0);
	  State* const stateR = face->getState(1);
	  
	  if (!isBFace) {
	    // on an internal face, we use a weighted average to have some sort of linear extrapolation
	    CFreal weightL =
	      1./MathFunctions::getDistance(face->getState(0)->getCoordinates(), xcFacePtr);
	    CFreal weightR =
	      1./MathFunctions::getDistance(face->getState(1)->getCoordinates(), xcFacePtr);
	    const CFreal weightSum = weightL + weightR; 
	    weightL /= weightSum;
	    weightR /= weightSum;
	    
	    for (CFuint i = 0; i < dim; ++i) {
	      const CFuint varID = m_variableIDs[i];
	      BPFSSFace[faceID][i] = (*stateL)[varID]*weightL + (*stateR)[varID]*weightR;
	    }	  
	  }
	  else {
	    // on a boundary face we constantly extrapolate the inner cell value
	    for (CFuint i = 0; i < dim; ++i) {
	      const CFuint varID = m_variableIDs[i];
	      BPFSSFace[faceID][i] = (*stateL)[varID];
	    }
	  }
	  
	  geoBuilder->releaseGE(); 
	}
      }
    }  
  }
  
  CFLog(VERBOSE, "ComputeFieldFromPotentialMHD::execute() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void ComputeFieldFromPotentialMHD::unsetup()
{
  CFAUTOTRACE;
  
  CellCenterFVMCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

