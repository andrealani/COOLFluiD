#include "Common/BadValueException.hh"
#include "Common/SafePtr.hh"
#include "Common/PE.hh"

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"
#include "Environment/FileHandlerInput.hh"

#include "Framework/PathAppender.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/LoadSourceData.hh"
#include "LinEuler/LinearizedEuler.hh"

#include "Framework/PhysicalModel.hh"

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/FileHandlerInput.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Physics {

      namespace LinearizedEuler {

RealVector LoadSourceData::m_indexTable = RealVector();

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<LoadSourceData,
                      DataProcessingData,
                      LinearizedEulerModule>
aLoadSourceDataProvider("LoadSourceData");
//////////////////////////////////////////////////////////////////////////////

void LoadSourceData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("InputFile","Name of Input File.");
  options.addConfigOption< std::string >("InputType","Type of Input File.");
  options.addConfigOption< CFuint >("AveragingPeriod","Number of iterations used for averaging.");
  options.addConfigOption< CFreal >("StartTime","Time at first time-step.");
  options.addConfigOption< CFreal >("avgStartTime","Time to start averaging.");
  options.addConfigOption< CFreal >("TimeStep","Time-step.");
  options.addConfigOption< CFuint >("Interpolation","Interpolation of meshes.");
}

//////////////////////////////////////////////////////////////////////////////

LoadSourceData::LoadSourceData( const std::string& name) :
  DataProcessingCom(name),
  socket_states("states"),
  socket_sourceData("source_Data"),
  socket_tecplot_states("tecplot_states"),
  socket_tecplot_coords("tecplot_coords")
{
  addConfigOptionsTo(this);

  m_InputFile = "sourceData";
  setParameter("InputFile",&m_InputFile);

  m_InputType = "CoolFluid";
  setParameter("InputType",&m_InputType);

  m_AveragingPeriod = 1;
  setParameter("AveragingPeriod",&m_AveragingPeriod);

  m_StartTime = 0;
  setParameter("StartTime",&m_StartTime);

  m_avgStartTime = 0;
  setParameter("avgStartTime",&m_avgStartTime);

  m_TimeStep = 1;
  setParameter("TimeStep",&m_TimeStep);

  m_Interpolation = 0;
  setParameter("Interpolation",&m_Interpolation);
}

//////////////////////////////////////////////////////////////////////////////

LoadSourceData::~LoadSourceData()
{
}

//////////////////////////////////////////////////////////////////////////////

void LoadSourceData::setup()
{

  DataProcessingCom::setup(); // first call setup of parent class

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealVector> sourceData = socket_sourceData.getDataHandle();
  DataHandle<RealVector> tecplot_states = socket_tecplot_states.getDataHandle();
  sourceData.resize(states.size());

  LoadSourceData::m_indexTable.resize(states.size());

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // Initialize the sourceData with zeros
  for (CFuint iState = 0; iState < sourceData.size(); ++iState)
  {
	sourceData[iState].resize(nbEqs+6); // (u',v',rhou'u'_avg,rhou'v'_avg,rhov'v'_avg,rho,u_avg,v_avg)
  }


  for (CFuint iState = 0; iState < sourceData.size(); ++iState)
  {
	sourceData[iState]=0.;
	LoadSourceData::m_indexTable[iState]=0;
  }


  //Interpolate with coordinates in first source file and mesh coordinate
  DataHandle<RealVector> tecplot_coords = socket_tecplot_coords.getDataHandle();
  std::ostringstream filename;
  if (m_InputType=="OFSourceData")
  	filename << m_InputFile << m_avgStartTime << ".plt";
  else 
	filename << m_InputFile << 0 << ".plt";
  boost::filesystem::path fname = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(filename.str());
  Common::SelfRegistPtr<FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<FileHandlerInput>::getInstance().create();
  
  ifstream& file = fhandle->open(fname);
  string dummy;
  float fdummy;
  char ch;
  CFuint NNode = 0;
  CFuint NElem = 0;

  for(int i=1; i<3; i++) getline (file,dummy);
  
  if (m_InputType=="OFSourceData")
  	for(int i=1; i<=6;i++) {file >> ch; }//from OpenFoam
  else 
	for(int i=1; i<24; i++) { file >>ch;};//from coolfluid
  
  file >> NNode;
  file >> ch >> ch >> ch; 
  file >> NElem;
  getline (file,dummy);

  const CFuint nbdim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbstates=states.size();
  const CFuint nbvariables=nbdim;
  tecplot_coords.resize(NNode);

  for (CFuint i = 0; i < NNode; i++){
    tecplot_coords[i].resize(nbdim);
  }

  cout << "Reading coordinates for interpolation from" << " Nodes: " << NNode << " Elements: " << NElem <<"... \n";
  for(CFuint i=0; i<NNode; i++){	
	for (CFuint j=0; j<nbdim; j++) {
		file >> fdummy;
		tecplot_coords[i][j]=fdummy;
	};
	file>>fdummy;
   	for (CFuint j=0; j<nbvariables; j++){
 		file >> fdummy;
	};
	if (m_InputType=="CoolFluid") file>>fdummy;
  }
  cout << "Done.\n";

  if (m_Interpolation==1){ 
  cout << "Interpolating meshes (grab a coffee, this can take some time)...";
  fhandle->close();
  // end interpolation

  //build index table for insertion of values in case of different node numbering / mesh structure
  for (CFuint iState = 0; iState < nbstates; ++iState) {
    Node& coord = states[iState]->getCoordinates();
    CFreal mindist = 10000.0;
    for (CFuint tecplotState = 0; tecplotState < NNode; ++tecplotState) {
      RealVector d=tecplot_coords[tecplotState]-coord;
      CFreal dist=d.norm2();
      if(dist<mindist) {
        mindist = dist;
        LoadSourceData::m_indexTable[iState]=tecplotState;
      }
    }
  }
  cout << "Done.\n" << "Index table written from " << filename.str() << " Nodes: " << NNode << " States: " << nbstates << "\n";

  	std::ostringstream filename;
  	filename << "index.dat";
  	cout << "Save index to " << filename <<".\n";
  	boost::filesystem::path fname = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(filename.str());
  	Common::SelfRegistPtr<FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<FileHandlerOutput>::getInstance().create();
  	ofstream& file = fhandle->open(fname);
  	for (CFuint iState = 0; iState < nbstates; ++iState) {
    		file<<m_indexTable[iState] << "\n";
    	};
  	fhandle->close();
  }
  else
  {
	cout << "Index table written from " << filename.str() << " Nodes: " << NNode << " States: " << nbstates << "\n";
	std::ostringstream filename;
  	filename << "index.dat";
  	cout << "Load index from " << filename.str() <<".\n";
  	boost::filesystem::path fname = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(filename.str());
  	Common::SelfRegistPtr<FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<FileHandlerInput>::getInstance().create();
  	ifstream& file = fhandle->open(fname);
  	for (CFuint iState = 0; iState < nbstates; ++iState) {
    		file>>m_indexTable[iState];
    	};
  	fhandle->close();
  };  //build index table

  // calculate averages of velocities
  for(CFuint k=1; k<=m_AveragingPeriod; k++){
 	std::ostringstream filename;
        if (m_InputType=="OFSourceData")
  		filename <<  m_InputFile << m_avgStartTime+(k-1)*m_TimeStep << ".plt";
  	else
		filename << m_InputFile << k << ".plt";
  	cout << "Averaging velocities with " << filename.str() << "...";
  	boost::filesystem::path fname = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(filename.str());
	Common::SelfRegistPtr<FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<FileHandlerInput>::getInstance().create();
  	ifstream& file = fhandle->open(fname);
 	NNode = 0;
 	NElem = 0;
  	for(int i=1; i<3; i++) getline (file,dummy);
  	if (m_InputType=="OFSourceData")
  		for(int i=1; i<=6;i++) {file >> ch; }//from OpenFoam
  	else
		for(int i=1; i<24; i++) { file >>ch;};//from coolfluid
  	file >> NNode;
  	file >> ch >> ch >> ch; 
	file >> NElem;
  	getline (file,dummy);
	tecplot_states.resize(NNode);
	for (CFuint i = 0; i < NNode; i++){
    		tecplot_states[i].resize(nbvariables+1); //u, v, rho
  	}
	for(CFuint i=0; i<NNode; i++){
		if (m_InputType != "OFSourceData" || k==1) 		
			for (CFuint j=0; j<nbdim; j++) {
 				file >> fdummy;
 			};
 		file >> fdummy;
    		for (CFuint j=0; j<nbvariables; j++){
  			file >> fdummy;
         		tecplot_states[i][j]=fdummy;
 		};
		if (m_InputType=="CoolFluid") file>>fdummy;
   	}
	for (CFuint iState = 0; iState < nbstates; ++iState) {
        	RealVector& sourceData_state = sourceData[iState];
		CFuint tecplotState=LoadSourceData::m_indexTable[iState];
		sourceData_state[6] = 0;// removed to test: 1./(k) * ((k-1)*sourceData_state[6]+tecplot_states[tecplotState][0]); // u_avg
		sourceData_state[7] = 0;// removed to test: 1./(k) * ((k-1)*sourceData_state[7]+tecplot_states[tecplotState][1]); // v_avg
		
   	}
	cout << "Done.\n";
   	fhandle->close();
   }//calculating averages of velocities


   //calculating averages of source terms
   for(CFuint k=1; k<=m_AveragingPeriod; k++){
 	std::ostringstream filename;
        if (m_InputType=="OFSourceData")
  		filename <<  m_InputFile << m_avgStartTime+(k-1)*m_TimeStep << ".plt";
  	else
		filename << m_InputFile << k << ".plt";
  	cout << "Averaging source terms with " << filename.str() << "...";
  	boost::filesystem::path fname = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(filename.str());
	Common::SelfRegistPtr<FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<FileHandlerInput>::getInstance().create();
  	ifstream& file = fhandle->open(fname);
 	NNode = 0;
 	NElem = 0;
  	for(int i=1; i<3; i++) getline (file,dummy);
  	if (m_InputType=="OFSourceData")
  		for(int i=1; i<=6;i++) {file >> ch; }//from OpenFoam
  	else
		for(int i=1; i<24; i++) { file >>ch;};//from coolfluid
  	file >> NNode;
  	file >> ch >> ch >> ch; 
	file >> NElem;
  	getline (file,dummy);
	tecplot_states.resize(NNode);
	for (CFuint i = 0; i < NNode; i++){
    		tecplot_states[i].resize(nbvariables+1); //u, v, rho
  	}
	for(CFuint i=0; i<NNode; i++){
		if (m_InputType != "OFSourceData" || k==1) 		
		for (CFuint j=0; j<nbdim; j++) {
 			file >> fdummy;
 		};
 		file >> fdummy;
		tecplot_states[i][nbvariables] = fdummy; // rho
    		for (CFuint j=0; j<nbvariables; j++){
  			file >> fdummy;
         		tecplot_states[i][j]=fdummy;
 		};
		if (m_InputType=="CoolFluid") file>>fdummy;
   	}
  	for (CFuint iState = 0; iState < nbstates; ++iState) {
        	RealVector& sourceData_state = sourceData[iState];
		CFuint tecplotState=LoadSourceData::m_indexTable[iState];
	
		sourceData_state[2] = 1./(k) * ((k-1)*sourceData_state[2]+tecplot_states[tecplotState][nbvariables]*(tecplot_states[tecplotState][0]-sourceData_state[6])*(tecplot_states[tecplotState][0]-sourceData_state[6])); // rhou'u'_avg
		sourceData_state[3] = 1./(k) * ((k-1)*sourceData_state[3]+tecplot_states[tecplotState][nbvariables]*(tecplot_states[tecplotState][0]-sourceData_state[6])*(tecplot_states[tecplotState][1]-sourceData_state[7])); // rhou'v'_avg
		sourceData_state[4] = 1./(k) * ((k-1)*sourceData_state[4]+tecplot_states[tecplotState][nbvariables]*(tecplot_states[tecplotState][1]-sourceData_state[7])*(tecplot_states[tecplotState][1]-sourceData_state[7]));  //rhov'v'_avg
	}
	cout << "Done.\n";
   	fhandle->close();
   }//calculating averages of source terms
}// setup()

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
LoadSourceData::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
LoadSourceData::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_sourceData);
  result.push_back(&socket_tecplot_states);
  result.push_back(&socket_tecplot_coords);
  return result;
}


//////////////////////////////////////////////////////////////////////////////

void LoadSourceData::execute(){

  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealVector> sourceData = socket_sourceData.getDataHandle();
  DataHandle<RealVector> tecplot_states = socket_tecplot_states.getDataHandle();
  CFreal time =SubSystemStatusStack::getActive()->getNbIter();
  
  /* read data from Tecplot file */
  std::ostringstream filename;
  if (m_InputType=="OFSourceData")
  	filename <<  m_InputFile << m_StartTime+(time)*m_TimeStep << ".plt";
  else
	filename << m_InputFile << (time) << ".plt";
  boost::filesystem::path fname = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(filename.str());
  Common::SelfRegistPtr<FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<FileHandlerInput>::getInstance().create();
  ifstream& file = fhandle->open(fname);
  
  string dummy; 
  float fdummy;
  char ch;
  CFuint NNode = 0;
  CFuint NElem = 0;

  /* header */
  for(int i=1; i<3; i++) getline (file,dummy);
  if (m_InputType=="OFSourceData")
  	for(int i=1; i<=6;i++) {file >> ch; }//from OpenFoam
  else
	for(int i=1; i<24; i++) { file >>ch;};//from coolfluid
  file >> NNode;
  file >> ch >> ch >> ch; 
  file >> NElem;
  getline (file,dummy);
  cout << "Data loaded from " << filename.str() << " Nodes: " << NNode << " Elements: " << NElem << "\n";
  const CFuint nbdim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbstates=states.size();
  const CFuint nbvariables=nbdim;

  /* variables */
  for(CFuint i=0; i<NNode; i++){
		if (m_InputType != "OFSourceData" || m_StartTime+(time)*m_TimeStep==m_avgStartTime) 		
		for (CFuint j=0; j<nbdim; j++) {
			file >> fdummy;
		};
		file >> fdummy;
		tecplot_states[i][nbvariables] = fdummy; // rho
		for (CFuint j=0; j<nbvariables; j++){
			file >> fdummy;
			tecplot_states[i][j]=fdummy;
		};
		if (m_InputType=="CoolFluid") file>>fdummy;
   }


  fhandle->close();

  for (CFuint iState = 0; iState < nbstates; ++iState) {
        RealVector& sourceData_state = sourceData[iState];
	CFuint tecplotState=LoadSourceData::m_indexTable[iState];
        sourceData_state[0] = tecplot_states[tecplotState][0] -sourceData_state[6]; //u'
	sourceData_state[1] = tecplot_states[tecplotState][1] -sourceData_state[7]; //v'
	sourceData_state[5] = tecplot_states[tecplotState][nbvariables]; // rho
   }
}//execute

//////////////////////////////////////////////////////////////////////////////

void LoadSourceData::unsetup()
{

  DataProcessingCom::unsetup(); // at last call setup of parent class
  DataHandle<RealVector> sourceData = socket_sourceData.getDataHandle();
  DataHandle<RealVector> tecplot_states = socket_tecplot_states.getDataHandle();
  DataHandle<RealVector> tecplot_coords = socket_tecplot_coords.getDataHandle();

  sourceData.resize(0);
  tecplot_states.resize(0);
  tecplot_coords.resize(0);

}

//////////////////////////////////////////////////////////////////////////////

void LoadSourceData::configure ( Config::ConfigArgs& args )
{
  DataProcessingCom::configure( args );

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////



