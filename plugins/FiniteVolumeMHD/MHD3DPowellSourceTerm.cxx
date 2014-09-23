#include "Framework/MethodStrategyProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/DataHandle.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/PathAppender.hh"
#include "Environment/DirPaths.hh"
#include "MHD/MHD3DVarSet.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "FiniteVolumeMHD/MHD3DPowellSourceTerm.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<MHD3DPowellSourceTerm,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeMHDModule>
mHD3DPowellSTFVMCCProvider("MHD3DPowellST");

//////////////////////////////////////////////////////////////////////////////

MHD3DPowellSourceTerm::MHD3DPowellSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  _varSet(CFNULL),
  socket_divBCellCenter("divBCellCenter"),
  socket_avgBxFace("avgBxFace"),
  socket_avgByFace("avgByFace"),
  socket_avgBzFace("avgBzFace"),
  _BDipole(),
  _divB(CFNULL),
  _divBMax(0.0),
  _divBMin(0.0),
  _saveRate(),
  _nameOutputFile(),
  _physicalData(),
  _dataLeftState(),
  _dataRightState()
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DPowellSourceTerm::~MHD3DPowellSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

boost::filesystem::path MHD3DPowellSourceTerm::constructFilename()
{
  boost::filesystem::path fpath(_nameOutputFile);
  fpath = PathAppender::getInstance().appendAllInfo( fpath );

  CFout << "Writing divB Errors to : " << fpath.string() << "\n";

  return fpath;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DPowellSourceTerm::prepareOutputFile(std::ofstream& outputFile)
{
  outputFile << "TITLE  =  divB Errors" << "\n";
  outputFile << "VARIABLES = elementID divB" << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DPowellSourceTerm::writeOutputFile()
{
  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // Execute and save file if needed...
  if(!(iter % _saveRate)) {

    SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    ofstream& outputFile = fhandle->open(constructFilename());

    prepareOutputFile(outputFile);
    
    DataHandle<State*, GLOBAL> state = _globalSockets.getSocketSink<State*>("states")->getDataHandle();
    const CFuint nbStates = state.size();
    
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      // Output to File
      outputFile << iState
		  << " "
		  << _divB[iState]
		  << "\n";
    }
    outputFile << "\ndivBMax = "
		  << _divBMax
		  << "\n"
		  << "divBMin = "
		  << _divBMin
		  << "\n";
    fhandle->close();
  }
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DPowellSourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DVarSet>();
  cf_assert(_varSet.isNotNull());
  
  _varSet->getModel()->resizePhysicalData(_physicalData);
  _varSet->getModel()->resizePhysicalData(_dataLeftState);
  _varSet->getModel()->resizePhysicalData(_dataRightState);
  _nameOutputFile =_varSet->getNameOutputFile();
  _saveRate = _varSet->getOutputFileSaveRate();
  
  DataHandle<State*, GLOBAL> state = _globalSockets.getSocketSink<State*>("states")->getDataHandle();
  const CFuint nbStates = state.size();

  DataHandle<CFreal> divBCellCenter  = socket_divBCellCenter.getDataHandle();
  divBCellCenter.resize(nbStates);

  _BDipole.resize(PhysicalModelStack::getActive()->getDim());

  _divB.resize(nbStates);
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DPowellSourceTerm::computeSource(Framework::GeometricEntity *const element,
					  RealVector& source,
					  RealMatrix& jacobian)
{
  cf_assert(_varSet.isNotNull());
  
  DataHandle<State*, GLOBAL> state = _globalSockets.getSocketSink<State*>("states")->getDataHandle();
  DataHandle<CFreal> normals = this->socket_normals.getDataHandle();
  DataHandle<CFint> isOutward =  this->socket_isOutward.getDataHandle();
  DataHandle<CFreal> divBCellCenter = socket_divBCellCenter.getDataHandle();
  DataHandle<CFreal> avgBxFace = socket_avgBxFace.getDataHandle();
  DataHandle<CFreal> avgByFace = socket_avgByFace.getDataHandle();
  DataHandle<CFreal> avgBzFace = socket_avgBzFace.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  const CFuint nbStates = state.size();

  // this source term is for MHD flows
  const vector<State*>* const states = element->getStates();
  const CFuint elementID = element->getID();

  // all elements in FVM should have only one state
  cf_assert(states->size() == 1);

  State *const currState = (*states)[0];

  _varSet->computePhysicalData(*currState, _physicalData);

  const GeomEntList *const faces = element->getNeighborGeos();
  const CFuint nbFacesInElem = faces->size();
  CFreal sumBnFace = 0.0;

  for (CFuint iFace = 0; iFace < nbFacesInElem; ++iFace) {

    const CFuint faceID = (*faces)[iFace]->getID();
    //State *const leftState = (*faces)[iFace]->getState(0);
    //State *const rightState = (*faces)[iFace]->getState(1);

    const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();

    // dimensional normal vector whose magnitude is equal to the face area

    CFreal nx = normals[startID];
    CFreal ny = normals[startID + 1];
    CFreal nz = normals[startID + 2];

    // to ensure that the normal is outward from the face

    if (isOutward[faceID] != static_cast<CFint>(elementID)) {
      nx *= -1.;
      ny *= -1.;
      nz *= -1.;
    }

    //_varSet->computePhysicalData(*leftState, _dataLeftState);
    //_varSet->computePhysicalData(*rightState, _dataRightState);

    //const CFreal bnFace = 0.5*
      //((_dataLeftState[MHDTerm::BX]*nx + _dataLeftState[MHDTerm::BY]*ny +
	//_dataLeftState[MHDTerm::BZ]*nz) +
       //(_dataRightState[MHDTerm::BX]*nx + _dataRightState[MHDTerm::BY]*ny +
	//_dataRightState[MHDTerm::BZ]*nz));
	
    const CFreal bnFace = avgBxFace[faceID]*nx + avgByFace[faceID]*ny + avgBzFace[faceID]*nz;

    sumBnFace += bnFace;

  }

  _divB[elementID] = sumBnFace / volumes[elementID];

  if (!getMethodData().isPerturb()) {
     divBCellCenter[elementID] = sumBnFace / volumes[elementID];
  }

  if ((!getMethodData().isPerturb()) && (elementID == (nbStates-1))) {
    _divBMax = 0.0;
    _divBMin = 0.0;
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      _divBMax = max(_divBMax,_divB[iState]);
      _divBMin = min(_divBMin,_divB[iState]);
    }
    writeOutputFile();
  }

  // In the source term, we take Btotal in the component corresponding to the momentum
  // equations according to Powell's JCP paper Vol.154 pp.284--309, 1999

  _BDipole = _varSet->getMagneticDipole(currState->getCoordinates()[XX],
					currState->getCoordinates()[YY],
					currState->getCoordinates()[ZZ]);
  const CFreal BxTotal = _physicalData[MHDTerm::BX] + _BDipole[0];
  const CFreal ByTotal = _physicalData[MHDTerm::BY] + _BDipole[1];
  const CFreal BzTotal = _physicalData[MHDTerm::BZ] + _BDipole[2];
  
  source[0] = 0.0;
  source[1] = -_divB[elementID]*BxTotal*volumes[elementID];
  source[2] = -_divB[elementID]*ByTotal*volumes[elementID];
  source[3] = -_divB[elementID]*BzTotal*volumes[elementID];
  source[4] = -_divB[elementID]*_physicalData[MHDTerm::VX]*volumes[elementID];
  source[5] = -_divB[elementID]*_physicalData[MHDTerm::VY]*volumes[elementID];
  source[6] = -_divB[elementID]*_physicalData[MHDTerm::VZ]*volumes[elementID];
  source[7] = -_divB[elementID]*(_physicalData[MHDTerm::VX]*_physicalData[MHDTerm::BX] +
		     _physicalData[MHDTerm::VY]*_physicalData[MHDTerm::BY] +
		     _physicalData[MHDTerm::VZ]*_physicalData[MHDTerm::BZ])*volumes[elementID];
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
MHD3DPowellSourceTerm::providesSockets()
{ 
  std::vector<Common::SafePtr<BaseDataSocketSource> > result =
          ComputeSourceTermFVMCC::providesSockets();
  result.push_back(&socket_divBCellCenter);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<SafePtr<BaseDataSocketSink> > MHD3DPowellSourceTerm::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = ComputeSourceTermFVMCC::needsSockets();
  result.push_back(&socket_avgBxFace);
  result.push_back(&socket_avgByFace);
  result.push_back(&socket_avgBzFace);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
