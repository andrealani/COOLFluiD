#include "LinEuler/LinearizedEuler.hh"
#include "LinEuler/ReadMeanFlowIncomp.hh"
#include "LinEuler/LinEulerTerm.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PathAppender.hh"
#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"


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

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ReadMeanFlowIncomp,
                      DataProcessingData,
                      LinearizedEulerModule>
aReadMeanFlowIncompProvider("ReadMeanFlowIncomp");

//////////////////////////////////////////////////////////////////////////////

void ReadMeanFlowIncomp::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("MeanFlowFile","File stores the mean flow in tecplot format.");
  options.addConfigOption< CFreal >("MeanDensity","Meanflow density (constant for incompressible flows.)");
  options.addConfigOption< CFreal >("BasePressure","Meanflow base pressure (to correct the gauge pressure.)");
  options.addConfigOption< bool >("Interpolate","Interpolate the meanflow to the LEE grid?)");
  options.addConfigOption< bool >("Write","Write the meanflow to file?)");
  options.addConfigOption< bool >("WriteCT","Write the connectivity table between grids?)");
  
}

//////////////////////////////////////////////////////////////////////////////
ReadMeanFlowIncomp::ReadMeanFlowIncomp(const std::string& name) :
  DataProcessingCom(name),
  socket_states("states"),
  socket_meanflow("meanflow"),
  socket_fluent_states("fluent_states"),
  socket_fluent_coords("fluent_coords"),
  nbNodesPerElemFluent()
{
  addConfigOptionsTo(this);
  
  setParameter("MeanFlowFile",&m_file_name);
  
  MeanDensity = 0.0;
  setParameter("MeanDensity",&MeanDensity);
  
  baseP = 0.0;
  setParameter("BasePressure",&baseP);  
  
  _Interpolate = false;
  setParameter("Interpolate",&_Interpolate);
  
  _Write = false;
  setParameter("Write",&_Write);
  
    _WriteCT = false;
  setParameter("WriteCT",&_WriteCT);
}

//////////////////////////////////////////////////////////////////////////////

ReadMeanFlowIncomp::~ReadMeanFlowIncomp()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ReadMeanFlowIncomp::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ReadMeanFlowIncomp::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_meanflow);
  result.push_back(&socket_fluent_states);
  result.push_back(&socket_fluent_coords);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ReadMeanFlowIncomp::configure( Config::ConfigArgs& args )

{
  DataProcessingCom::configure(args);

  if(m_file_name.empty())
     throw BadValueException(FromHere(),"ReadMeanFlowIncomp::setFile(): no tecplot file for mean flow provided.");

}

//////////////////////////////////////////////////////////////////////////////

void ReadMeanFlowIncomp::setup()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  DataHandle<RealVector> meanflow = socket_meanflow.getDataHandle();
  meanflow.resize(states.size());

  const CFuint nb_vars = PhysicalModelStack::getActive()->getNbEq();

  for (CFuint i = 0; i < states.size(); ++i)
    meanflow[i].resize(nb_vars);
  
  executeOnTrs();
  
}

//////////////////////////////////////////////////////////////////////////////

void ReadMeanFlowIncomp::executeOnTrs()
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMin( "Meanflow::executeOnTrs() called for TRS: " << trs->getName() << "\n");

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealVector> meanflow = socket_meanflow.getDataHandle();
  DataHandle<RealVector> fluent_states = socket_fluent_states.getDataHandle();
  DataHandle<RealVector> fluent_coords = socket_fluent_coords.getDataHandle();
  
  const CFuint nbdim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbstates=PhysicalModelStack::getActive()->getNbEq();

/* read data from Fluent-Tecplot file */

  boost::filesystem::path fname = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(m_file_name);
  Common::SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();

  ifstream& file = fhandle->open(fname);

  std::string dummy;
  char ch;
  CFuint NNode = 0;

/* header */
  getline (file,dummy);

  const CFuint ns = 2*nbdim + 2;
  for(CFuint i=1; i< ns; i++) {
    getline (file,dummy);
  }
  
  getline (file,dummy);
  getline (file,dummy);

  file >> ch >> ch; 
  file >> NNode;
  
  fluent_coords.resize(NNode);
  fluent_states.resize(NNode);

  for (CFuint i = 0; i < NNode; ++i) {
    fluent_states[i].resize(nbstates);
    fluent_coords[i].resize(nbdim);
  }

  for(CFint i=1; i<4; i++)
    getline (file,dummy);

/* variables */
  for(CFuint dim=0; dim<nbdim; dim++) {
    for(CFuint i=0; i<NNode; i++)
      file >> fluent_coords[i][dim];
  }
  
  //density is not loaded
  for(CFuint st=1; st<nbstates; st++) {
    for(CFuint i=0; i<NNode; i++)
      file >> fluent_states[i][st];
    }

/* connectivity */
//   nbNodesPerElemFluent.resize(NElem);
//   for(CFuint i=0; i<NElem; i++)
//     nbNodesPerElemFluent[i]=3;
// 
//   Common::ConnectivityTable<CFuint> FluentMesh;
//   FluentMesh.resize(nbNodesPerElemFluent);
// 
// 
// 
//   for(CFuint i=0; i<NElem; i++) {
//     for(CFuint j=0; j<3; j++)
//       file >> (FluentMesh) (i, j);
//   }

  fhandle->close();

if (_Interpolate) {
  CTable.resize(states.size());
  
/* interpolate data to the RDS grid */
  for (CFuint iState = 0; iState < states.size(); ++iState) {
    Node& coord = states[iState]->getCoordinates();
    CFreal mindist = 10000.0;
    for (CFuint FluentState = 0; FluentState < NNode; ++FluentState) {
      RealVector d=fluent_coords[FluentState]-coord;
      CFreal dist=d.norm2();
      if(dist<mindist) {
	CTable[iState] = FluentState;
        mindist = dist;
        RealVector& meanflow_state = meanflow[iState];
        meanflow_state[nbdim+1] = fluent_states[FluentState][nbdim+1]+101325.;
	meanflow_state[0] = MeanDensity;
	for (CFuint iDim = 1; iDim < nbdim+1; iDim++)
	  meanflow_state[iDim] = fluent_states[FluentState][iDim];
      }
    }
  }
  
// write connectivity table for the source interpolation
 
 if(_WriteCT) {
  writeConnectivityFile();
 }
 
}
else {
    for (CFuint iState = 0; iState < states.size(); ++iState) {
      RealVector& meanflow_state = meanflow[iState];
      meanflow_state[nbdim+1] = fluent_states[iState][nbdim+1]+ baseP;
      meanflow_state[0] = MeanDensity;
	for (CFuint iDim = 1; iDim < nbdim+1; iDim++)
	  meanflow_state[iDim] = fluent_states[iState][iDim];
      }
}

if(_Write) {
  writeMeanflowFile();
}


  SafePtr<LinEulerTerm> lterm = PhysicalModelStack::getActive()->getImplementor()->	getConvectiveTerm().d_castTo<LinEulerTerm>();

  lterm->setMeanFlowArray(meanflow);
  
  
}

//////////////////////////////////////////////////////////////////////////////
void ReadMeanFlowIncomp::writeMeanflowFile()
{

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbdim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbstates=PhysicalModelStack::getActive()->getNbEq();
    
  boost::filesystem::path fname = Environment::DirPaths::getInstance().getResultsDir() / ("meanflowCF.dat");
  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  
  DataHandle<RealVector> meanflow = socket_meanflow.getDataHandle();

  ofstream& file = fhandle->open(fname);
  
  CFuint nbnodes = states.size();

  std::string onestring;
  onestring = "DUMMY";
  
  char ch;
  ch = 'X';
  
  /* header */
  file << onestring << "\n";

  const CFuint ns = 2*nbdim + 2;
  for(CFuint i=1; i< ns; i++) {
    file << onestring << "\n";
  }
  
  file << onestring << "\n";
  file << onestring << "\n";

  file << ch << ch; 
  file << nbnodes;
 
  for(CFint i=1; i<4; i++)
    file << onestring << "\n";
  
  /* coordinates */
  for(CFuint dim=0; dim<nbdim; dim++) {
    for(CFuint i=0; i<nbnodes; i++) {
      Node& coord = states[i]->getCoordinates();
      file << coord[dim] << "\t";
    }
   file << "\n";
    
  }
  file << "\n";
   
//density is not written
  for(CFuint st=1; st<nbstates; st++) {
    for(CFuint i=0; i<nbnodes; i++)
      file << meanflow[i][st] << "\t";
    file << "\n";
    }

  fhandle->close();


}


//////////////////////////////////////////////////////////////////////////////
void ReadMeanFlowIncomp::writeConnectivityFile()
{

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
    
  boost::filesystem::path fname = Environment::DirPaths::getInstance().getResultsDir() / ("CTable.dat");
  Common::SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  
  ofstream& file = fhandle->open(fname);
  
  for(CFuint st=0; st<states.size(); st++) {
      file << CTable[st] << "\n";
    }

  fhandle->close();  
}

//////////////////////////////////////////////////////////////////////////////


void ReadMeanFlowIncomp::unsetup()
{
  DataHandle<RealVector> meanflow = socket_meanflow.getDataHandle();
  DataHandle<RealVector> fluent_states = socket_fluent_states.getDataHandle();
  DataHandle<RealVector> fluent_coords = socket_fluent_coords.getDataHandle();

  meanflow.resize(0);
  fluent_states.resize(0);
  fluent_coords.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinearizedEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

