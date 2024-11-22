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
#include "FiniteVolumePoissonNEQ/ComputeBField.hh"
#include "FiniteVolumePoissonNEQ/FiniteVolumePoissonNEQ.hh"

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
      
MethodCommandProvider<ComputeBField, 
		      DataProcessingData, 
		      FiniteVolumePoissonNEQModule>
computeBFieldProvider("ComputeBField");

//////////////////////////////////////////////////////////////////////////////

void ComputeBField::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption< std::string >("DefFileName","Name of file where Functions are defined.");
  options.addConfigOption< std::string >("DataFileName","Name of file where Bx By Bz are given."); // Vatsalya: NEw
  options.addConfigOption< std::string >("OtherNamespace","Name of OtherNamespace where Bx By Bz are present."); // Vatsalya: NEw
}

//////////////////////////////////////////////////////////////////////////////

ComputeBField::ComputeBField(const std::string& name) :
  Framework::DataProcessingCom(name),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_gstates("gstates"),
//  socket_otherBX("BX"),
//  socket_otherBY("BY"),
//  socket_otherBZ("BZ"),
  socket_faceCenters("faceCenters"),
  socket_Bfield("Bfield"),
  socket_BfieldFaces("BfieldFaces"), 
  socket_otherStates("states"),
  m_input()
{
  this->addConfigOptionsTo(this);

  m_functions = std::vector<std::string>();
  setParameter("Def",&m_functions);
  
  m_vars = std::vector<std::string>();
  setParameter("Vars",&m_vars);

  m_functionsFileName = "";
  setParameter("DefFileName",&m_functionsFileName);

  m_dataFileName = "";
  setParameter("DataFileName",&m_dataFileName);

  m_otherNamespace = "";
  setParameter("OtherNamespace", &m_otherNamespace);
  
}
      
//////////////////////////////////////////////////////////////////////////////

ComputeBField::~ComputeBField()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > 
ComputeBField::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_faceCenters);
  result.push_back(&socket_otherStates);
 // result.push_back(&socket_otherBX);
 // result.push_back(&socket_otherBY);
 // result.push_back(&socket_otherBZ);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > 
ComputeBField::providesSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > result;
  result.push_back(&socket_Bfield);
  result.push_back(&socket_BfieldFaces);
 
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeBField::setup()
{
  CFAUTOTRACE;
  
  DataProcessingCom::setup();
  
  
  
  // B in cells
  const CFuint nbStates = socket_states.getDataHandle().size();

  m_input.resize(3*nbStates);  // Vatsalya : This is new, *3 to resize it for B

  Framework::DataHandle< CFreal> Bfield = socket_Bfield.getDataHandle();
  Bfield.resize(nbStates*3);
  
  // B in faces
  const CFuint nbFaces = MeshDataStack::getActive()->Statistics().getNbFaces();
  Framework::DataHandle< CFreal> BfieldFaces = socket_BfieldFaces.getDataHandle();
  BfieldFaces.resize(nbFaces*3);
  
}
      
//////////////////////////////////////////////////////////////////////////////

void ComputeBField::unsetup()
{
  CFAUTOTRACE;
  
  Framework::DataProcessingCom::unsetup();
}
      
//////////////////////////////////////////////////////////////////////////////

void ComputeBField::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);
  
  if (m_functionsFileName != "" && m_functions.size() == 0) {
    const std::string name = getMethodData().getNamespace();
    Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
      (SubSystemStatusStack::getCurrentName()).getNamespace(name);
    runSerial<void, ComputeBField, &ComputeBField::readFunctionsFile>(this, name);
  }  
  m_vFunction.setFunctions(m_functions);
  m_vFunction.setVariables(m_vars);
  try {
    m_vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
  if(m_otherNamespace != ""){
  CFLog(INFO, "ComputeFieldFromPotential::configure() => m_otherNamespace = " <<	m_otherNamespace << "\n");
 // socket_otherBX.setDataSocketNamespace(m_otherNamespace);
 // socket_otherBY.setDataSocketNamespace(m_otherNamespace);
 // socket_otherBZ.setDataSocketNamespace(m_otherNamespace);
	socket_otherStates.setDataSocketNamespace(m_otherNamespace);
  
}
  
}
      
//////////////////////////////////////////////////////////////////////////////

void ComputeBField::execute()
{
  Framework::DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
  Framework::DataHandle< CFreal> Bfield = socket_Bfield.getDataHandle();


  const CFuint nbStates = states.size(); 
  const CFuint totState = nbStates*3 ;
  if(m_dataFileName != "" && m_otherNamespace == ""){    // Vatsalya : New addition. datafile is present, namespace is not present
     this->readFile();   // put this in if else
  // fill in the cell centered B field values
    cout << m_input.size() << "\t = m_input.size() \n";
    cout<< "m_input[3] = \t " << m_input[3] <<"\n";
    for (CFuint startID = 0; startID < totState; startID++) {   
      Bfield[startID]   = m_input[startID];
    }
    cout<< "Bfield[3] = \t " << Bfield[3] <<"\n";
  }
  if(m_dataFileName == "" && m_otherNamespace == ""){  // no datafile and no other namespace, use expression
    // fill in the cell centered B field values
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      const CFuint startID = iState*3;
      m_vFunction.evaluate(states[iState]->getCoordinates(), m_input);
      
      Bfield[startID]   = m_input[0];
      Bfield[startID+1] = m_input[1];
      Bfield[startID+2] = m_input[2];
     // CFLog(INFO, "[Bx, By, By] = " << m_input << "\n");
    }
	cout<< "Bfield[3] = \t " << Bfield[3] <<"\n";
  }
  if(m_dataFileName == "" && m_otherNamespace != ""){  // Vatsalya : no data file but copy from othernamespace
  //  this->CopyFromNameSpace(); // copy data from namespace
    m_input.resize(3*(nbStates));
   	// call data from other namespace

	DataHandle<State*, GLOBAL> otherStates = socket_otherStates.getDataHandle();
	const CFuint nbOtherStates = otherStates.size();
	//cout<<"nbOtherStates = " << nbOtherStates <<"\n";
	//cout<<"nbStates = " << nbStates <<"\n";
	//cout<<"SubSystemStatusStack::getActive()->getName() = " << SubSystemStatusStack::getActive()->getName() <<"\n";

	//const RealVector& linearData = getModel()->getPhysicalData();
	//cout<<"*states[iState])[xVar] = " << (*otherStates[0])[0] <<"\n"; // istate is xVar is variable number
	
 
   // cout << m_input.size() << "\t = m_input.size() \n";
   // cout<< "m_input[3] = \t " << m_input[3] <<"\n";
    for (CFuint iState = 0; iState < nbOtherStates; ++iState) {
      const CFuint startID = iState*3;
        Bfield[startID]   = ((*otherStates[iState])[0])/1000; //m_input[startID];
	Bfield[startID+1]   = ((*otherStates[iState])[1])/1000;
	Bfield[startID+2]   = ((*otherStates[iState])[2])/1000;
    }
  //  cout<< "Bfield[3] = \t " << Bfield[3] <<"\n";
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
}

///////////////////////////////////////////////////////////////

void ComputeBField::readFunctionsFile()
{
  boost::filesystem::path fpath(m_functionsFileName);

  CFLog(INFO, "ComputeBField::readFunctionsFile() " << m_functionsFileName << "\n");
  
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
      CFLog(VERBOSE, "ComputeBField::readFunctionsFile() => [" << functionName << "]\n");
      functionName = "";
    }
  }
  
  fhandle->close();
}

 void ComputeBField::readFile()
{
  boost::filesystem::path fpath(m_dataFileName);

  CFLog(INFO, "ComputeBField::readFile() " << m_dataFileName << "\n");
  
  SelfRegistPtr<Environment::FileHandlerInput> fhandle = 
    Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
    
  string line = "";
  CFuint nbLines = 0;
  ifstream& inputFile1 = fhandle->open(fpath);
  //Vatsalya: The following don't work here, I need to get the sockets in function definition, else use proxy m_input
 // Framework::DataHandle< CFreal> Bfield = socket_Bfield.getDataHandle();
 // Framework::DataHandle<Framework::State*, Framework::GLOBAL> states = socket_states.getDataHandle();
 // const CFuint nbStates = states.size(); 

  if (inputFile1.is_open()) {
    while (!inputFile1.eof()) {
      getline (inputFile1,line);
      ++nbLines;
    }
  }
  cout<<"nblines = "<< nbLines<<"\n";
  //cout<<"nbStates = "<< nbStates<<"\n";
  fhandle->close();
  double B_x, B_y, B_z;  
  ifstream& inputFile2 = fhandle->open(fpath);
  // Vatsalya: This should be tru --> nbLines = iState
  //assert(nbLines == nbStates);
  CFLog(INFO, "ComputeBField::readFile() " << m_dataFileName << "\n");
  //cout << m_functions.size() << "\t = m_functions.size() \n";
  m_input.resize(3*(nbLines-1));
  cout << m_input.size() << "\t = m_input.size() \n";
  //const CFuint startID = 0;
  int startID = 0;
  while( inputFile2 >> B_x >> B_y >> B_z){
  // cout << B_x <<" = B_X \t"<< B_y<<" = B_y \t" << B_z<<"= B_z \n" ;
    m_input[startID] = B_x/1000; // vatsalya : Magpylib gives in milli-Tesla
    startID++;
    m_input[startID] = B_y/1000;
    startID++;
    m_input[startID] = B_z/1000;
    startID++;
  }
  cout<< "m_input[3] = \t " << m_input[3] <<"\n";
  cout<< "startID = \t " << startID <<"\n";
  cout<< "m_input[3] = \t " << m_input[3] <<"\n";
  fhandle->close();
}
void ComputeBField::CopyFromNameSpace()
{
 CFLog(INFO, "ComputeBField::CopyFromNameSpace() " << m_dataFileName << "\n");


}
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////