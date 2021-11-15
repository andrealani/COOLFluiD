#include "Common/PE.hh"
#include "Common/EventHandler.hh"
#include "Common/PEFunctions.hh"

#include "Environment/FileHandlerInput.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "MathTools/MathFunctions.hh"
#include "MathTools/MathConsts.hh"

#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "FiniteVolume/CellData.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolume/ComputeBFieldCom.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////
      
MethodCommandProvider<ComputeBFieldCom, 
		      DataProcessingData, 
		      FiniteVolumeModule>
computeBFieldComProvider("ComputeBFieldCom");

//////////////////////////////////////////////////////////////////////////////

void ComputeBFieldCom::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFuint> >
    ("VariableIDs", "IDs of the variable to be assigned to the newly computed field.");
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption< std::string >("DefFileName","Name of file where Functions are defined.");
  options.addConfigOption< string >
    ("OtherNamespace", "Name of the other namespace (providing the potential).");
}

//////////////////////////////////////////////////////////////////////////////

ComputeBFieldCom::ComputeBFieldCom(const std::string& name) :
  Framework::DataProcessingCom(name),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_otherStates("states"),
  socket_gstates("gstates"),
  socket_faceCenters("faceCenters"),
  socket_Bfield("Bfield"),
  socket_BfieldFaces("BfieldFaces"),
  m_input()
{
  this->addConfigOptionsTo(this);

  m_variableIDs = vector<CFuint>();
  setParameter("VariableIDs",&m_variableIDs);
  
  m_functions = std::vector<std::string>();
  setParameter("Def",&m_functions);
  
  m_vars = std::vector<std::string>();
  setParameter("Vars",&m_vars);

  m_functionsFileName = "";
  setParameter("DefFileName",&m_functionsFileName);
  
  m_otherNamespace = "";
  setParameter("OtherNamespace", &m_otherNamespace);
}
      
//////////////////////////////////////////////////////////////////////////////

ComputeBFieldCom::~ComputeBFieldCom()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > 
ComputeBFieldCom::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_otherStates);
  result.push_back(&socket_gstates);
  result.push_back(&socket_faceCenters);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > 
ComputeBFieldCom::providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result;
  result.push_back(&socket_Bfield);
  result.push_back(&socket_BfieldFaces);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeBFieldCom::setup()
{
  CFAUTOTRACE;
  
  DataProcessingCom::setup();
    
  m_input.resize(3);

  if (m_variableIDs.size() == 0) {
    m_variableIDs.resize(3);
    for (CFuint i = 0; i < 3; ++i) {m_variableIDs[i] = i;}
  }
  
  // B in cells
  const CFuint nbStates = socket_states.getDataHandle().size();
  Framework::DataHandle< CFreal> Bfield = socket_Bfield.getDataHandle();
  Bfield.resize(nbStates*3);
  
  // B in faces
  const CFuint nbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  Framework::DataHandle< CFreal> BfieldFaces = socket_BfieldFaces.getDataHandle();
  BfieldFaces.resize(nbFaces*3);
}
      
//////////////////////////////////////////////////////////////////////////////

void ComputeBFieldCom::unsetup()
{
  CFAUTOTRACE;
  
  Framework::DataProcessingCom::unsetup();
}
      
//////////////////////////////////////////////////////////////////////////////

void ComputeBFieldCom::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

  if (m_functionsFileName != "" && m_functions.size() == 0) {
    cf_assert(m_otherNamespace == "");
    const std::string name = getMethodData().getNamespace();
    Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
      (SubSystemStatusStack::getCurrentName()).getNamespace(name);
    runSerial<void, ComputeBFieldCom, &ComputeBFieldCom::readFunctionsFile>(this, name);
  }
  
  if (m_functions.size() > 0) {
    m_vFunction.setFunctions(m_functions);
    m_vFunction.setVariables(m_vars);
    try {
      m_vFunction.parse();
    }
    catch (Common::ParserException& e) {
      CFout << e.what() << "\n";
      throw; // retrow the exception to signal the error to the user
    }
  }
  
  if (m_otherNamespace != "") {
    cf_assert(m_functions.size() == 0);
    CFLog(VERBOSE, "ComputeBFieldCom::configure() => m_otherNamespace = " <<
	  m_otherNamespace << "\n");
    socket_otherStates.setDataSocketNamespace(m_otherNamespace);
  }
}
      
//////////////////////////////////////////////////////////////////////////////

void ComputeBFieldCom::execute()
{
  CFLog(VERBOSE, "ComputeBFieldCom::execute() => START\n");
  
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  Framework::DataHandle< CFreal> Bfield = socket_Bfield.getDataHandle();
  const CFuint nbStates = states.size(); 
  
  if (m_functions.size() > 0) {
    computeFieldFromFunction();
  }

  if (m_otherNamespace != "") {
    computeFieldFromOtherStates();
  }
  
  // fill in the face centered B field values
  DataHandle<CFreal> faceCenters = socket_faceCenters.getDataHandle();
  DataHandle<CFreal> BfieldFaces = socket_BfieldFaces.getDataHandle();
  
  // prepare the building of the faces
  GeometricEntityPool<FaceTrsGeoBuilder> geoBuilder;
  geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  geoBuilder.setup();
  
  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder.getDataGE();
  
  // get all the TRSs
  vector<SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  const CFuint nbTRSs = trs.size();
  const CFuint dim    = PhysicalModelStack::getActive()->getDim();
  const CFuint dimBBB = DIM_3D; // always 3 components for B
  
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
	GeometricEntity *const face = geoBuilder.buildGE();
	const vector<Node*>& nodes = *face->getNodes();
	const CFuint nbNodesInFace = nodes.size();
	const CFuint faceID  = face->getID();
	const CFuint startFaceID = faceID*dim;
	xcFacePtr.wrap(dim, &faceCenters[startFaceID]);

	const CFuint startID = faceID*dimBBB;
	const CFuint lstateID = face->getState(0)->getLocalID()*dimBBB;
	CFLog(DEBUG_MIN, "BfieldFaces[" << faceID << "] = ");
	
	if (!isBFace) {
	  // on an internal face, we use a weighted average to have some sort of linear extrapolation
	  CFreal weightL = 1./MathFunctions::getDistance(face->getState(0)->getCoordinates(), xcFacePtr);
	  CFreal weightR = 1./MathFunctions::getDistance(face->getState(1)->getCoordinates(), xcFacePtr);
	  const CFreal weightSum = weightL + weightR; 
	  weightL /= weightSum;
	  weightR /= weightSum;
	  const CFuint rstateID = face->getState(1)->getLocalID()*dimBBB;
	  
	  for (CFuint i = 0; i < dimBBB; ++i) {
	    const CFuint startIDi = startID + i;
	    cf_assert(startIDi < BfieldFaces.size());
	    cf_assert(lstateID+i < Bfield.size());
	    cf_assert(rstateID+i < Bfield.size());
	    BfieldFaces[startIDi] = Bfield[lstateID+i]*weightL + Bfield[rstateID+i]*weightR;
	    CFLog(DEBUG_MIN, BfieldFaces[startIDi] << " ");
	  }
	}
	else {
	  // on a boundary face we constantly extrapolate the inner cell value
	  for (CFuint i = 0; i < dimBBB; ++i) {
	    const CFuint startIDi = startID + i;
	    cf_assert(startIDi < BfieldFaces.size());
	    cf_assert(lstateID+i < Bfield.size());
	    BfieldFaces[startIDi] = Bfield[lstateID+i];
	    CFLog(DEBUG_MIN, BfieldFaces[startIDi] << " ");
	  }
	}
	CFLog(DEBUG_MIN, "\n");
	
	geoBuilder.releaseGE(); 
      }
    }
  }

  CFLog(VERBOSE, "ComputeBFieldCom::execute() => END\n"); 
}

///////////////////////////////////////////////////////////////

void ComputeBFieldCom::readFunctionsFile()
{
  boost::filesystem::path fpath(m_functionsFileName);

  CFLog(INFO, "ComputeBFieldCom::readFunctionsFile() " << m_functionsFileName << "\n");
  
  SelfRegistPtr<Environment::FileHandlerInput> fhandle = 
    Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
    
  string line = "";
  CFuint nbLines = 0;
  ifstream& inputFile1 = fhandle->open(fpath);
  
  if (inputFile1.is_open()) {
    while (!inputFile1.eof()) {
      getline (inputFile1,line);
      ++nbLines;
    }
  }
  fhandle->close();
    
  ifstream& inputFile2 = fhandle->open(fpath);
  string functionName = "";
  for (CFuint iLine = 0; iLine < nbLines-1; ++iLine) {
    getline(inputFile2, line);
    // if you find "\", remove it and update the string
    if (line.find('\\') != std::string::npos) {
      CFLog(VERBOSE, "BEFORE InitState::readFunctionsFile() => " << line << "\n");
      line.erase(remove(line.begin(), line.end(), '\\'), line.end());
      CFLog(VERBOSE, "AFTER  InitState::readFunctionsFile() => " << line << "\n");
      StringOps::trim(line);
      functionName += line;
    }
    else {
      StringOps::trim(line);
      functionName += line;
      m_functions.push_back(functionName);
      CFLog(VERBOSE, "ComputeBFieldCom::readFunctionsFile() => [" << functionName << "]\n");
      functionName = "";
    }
  }
  
  fhandle->close();
 }

//////////////////////////////////////////////////////////////////////////////

void ComputeBFieldCom::computeFieldFromFunction()
{
  CFLog(VERBOSE, "ComputeBFieldCom::computeFieldFromFunction() => START\n");
  
  Framework::DataHandle< CFreal> Bfield = socket_Bfield.getDataHandle();
  const CFuint nbStates = socket_states.getDataHandle().size(); 
  
  // fill in the cell centered B field values
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    const CFuint startID = iState*3;
    
    Bfield[startID]   = m_input[0];
    Bfield[startID+1] = m_input[1];
    Bfield[startID+2] = m_input[2];
    CFLog(VERBOSE, "ComputeBFieldCom::computeFieldFromFunction() => [ Bx  By Bz ] = [ "
	  << m_input[0] << " " << m_input[1] << " " << m_input[2] << "]\n");
  }
}

//////////////////////////////////////////////////////////////////////////////

void ComputeBFieldCom::computeFieldFromOtherStates()
{
  CFLog(INFO, "ComputeBFieldCom::computeFieldFromFunction() => START\n");
  
  DataHandle<State*, GLOBAL> otherStates = socket_otherStates.getDataHandle();
  Framework::DataHandle< CFreal> Bfield = socket_Bfield.getDataHandle();
  const CFuint nbStates = socket_states.getDataHandle().size(); 
  
  cf_assert(nbStates == otherStates.size());
  cf_assert(m_variableIDs.size() == 3);
  
  // fill in the cell centered B field values
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    const CFuint startID = iState*3;
    CFLog(VERBOSE, "ComputeBFieldCom::computeFieldFromOtherStates() => [ Bx By Bz ] = [ ");
    for (CFuint i = 0; i < 3; ++i) {
      Bfield[startID+i] = (*otherStates[iState])[m_variableIDs[i]];
      CFLog(VERBOSE, Bfield[startID+i] << " ");
    }
    CFLog(VERBOSE, "]\n");
  }
}

//////////////////////////////////////////////////////////////////////////////
    
    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
