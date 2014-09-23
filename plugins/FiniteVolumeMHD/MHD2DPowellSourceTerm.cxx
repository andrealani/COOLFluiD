#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/DataHandle.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/PathAppender.hh"
#include "MHD/MHD2DVarSet.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "FiniteVolumeMHD/MHD2DPowellSourceTerm.hh"
#include "Framework/MethodStrategyProvider.hh"

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

MethodStrategyProvider<MHD2DPowellSourceTerm,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeMHDModule>
mHD2DPowellSTFVMCCProvider("MHD2DPowellST");

//////////////////////////////////////////////////////////////////////////////

MHD2DPowellSourceTerm::MHD2DPowellSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  socket_divBCellCenter("divBCellCenter"),
  socket_avgBxFace("avgBxFace"),
  socket_avgByFace("avgByFace"),
  _divB(),
  _varSet(CFNULL),
  _BDipole(),
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

MHD2DPowellSourceTerm::~MHD2DPowellSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

boost::filesystem::path MHD2DPowellSourceTerm::constructFilename()
{
  boost::filesystem::path fpath(_nameOutputFile);
  fpath = PathAppender::getInstance().appendAllInfo( fpath );

  CFout << "Writing divB Errors to : " << fpath.string() << "\n";

  return fpath;
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DPowellSourceTerm::prepareOutputFile(std::ofstream& outputFile)
{
  outputFile << "TITLE  =  divB Errors" << "\n";
  outputFile << "VARIABLES = elementID divB" << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DPowellSourceTerm::writeOutputFile()
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

void MHD2DPowellSourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<MHD2DVarSet>();
  
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

void MHD2DPowellSourceTerm::computeSource(Framework::GeometricEntity *const element,
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
    const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
    
    // dimensional normal vector whose magnitude is equal to the face length
    
    CFreal nx = normals[startID];
    CFreal ny = normals[startID + 1];
    
    // to ensure that the normal is outward from the face
    
    if (isOutward[faceID] != static_cast<CFint>(elementID)) {
      nx *= -1.;
      ny *= -1.;
    }

    const CFreal bnFace = avgBxFace[faceID]*nx + avgByFace[faceID]*ny;	
    sumBnFace += bnFace;
  }
  
  _divB[elementID] = sumBnFace / volumes[elementID];
  
  if (!getMethodData().isPerturb()) {
    divBCellCenter[elementID] = sumBnFace / volumes[elementID];
  }
  
  if (!getMethodData().isPerturb()) {
    static CFuint counter = 0;
    if (counter == nbStates-1) {
      _divBMax = 0.0;
      _divBMin = 0.0;
      for (CFuint iState = 0; iState < nbStates; ++iState) {
	_divBMax = max(_divBMax,_divB[iState]);
	_divBMin = min(_divBMin,_divB[iState]);
      }
      writeOutputFile();
      counter = 0;
    }
    else {
      counter++;
    }
  }
  
  // In the source term, we take Btotal in the component corresponding to the momentum
  // equations according to Powell's JCP paper Vol.154 pp.284--309, 1999
  
  _BDipole = _varSet->getMagneticDipole(currState->getCoordinates()[XX], currState->getCoordinates()[YY]);
  const CFreal BxTotal = _physicalData[MHDTerm::BX] + _BDipole[0];
  const CFreal ByTotal = _physicalData[MHDTerm::BY] + _BDipole[1];
  
  source[0] = 0.0;
  source[1] = -_divB[elementID]*BxTotal*volumes[elementID];
  source[2] = -_divB[elementID]*ByTotal*volumes[elementID];
  source[3] = -_divB[elementID]*_physicalData[MHDTerm::BZ]*volumes[elementID];
  source[4] = -_divB[elementID]*_physicalData[MHDTerm::VX]*volumes[elementID];
  source[5] = -_divB[elementID]*_physicalData[MHDTerm::VY]*volumes[elementID];
  source[6] = -_divB[elementID]*_physicalData[MHDTerm::VZ]*volumes[elementID];
  source[7] = -_divB[elementID]*(_physicalData[MHDTerm::VX]*_physicalData[MHDTerm::BX] +
				 _physicalData[MHDTerm::VY]*_physicalData[MHDTerm::BY] +
				 _physicalData[MHDTerm::VZ]*_physicalData[MHDTerm::BZ])*volumes[elementID];

}
  
//////////////////////////////////////////////////////////////////////////////

//void MHD2DPowellSourceTerm::computeSource(GeometricEntity *const element,
					 // RealVector& source)
//{
  /// My Implementation
  
  //cf_assert(_varSet.isNotNull());
  
  //DataHandle<State*> state = _sockets.getSocketSink<State*>("states")->getDataHandle();
  // DataHandle<CFreal> normals = this->socket_normals.getDataHandle();
  //DataHandle<CFint> isOutward = this->socket_isOutward.getDataHandle();
  //DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  //DataHandle<RealVector> nstates = _sockets.getSocketSink<RealVector>("nstates")->getDataHandle();
//   const CFuint nbStates = state.size();
  
  // this source term is for MHD flows
  //const vector<State*>* const states = element->getStates();
  //const CFuint elementID = element->getID();
  
  // all elements in FVM should have only one state
  //cf_assert(states->size() == 1);
  
  //State *const currState = (*states)[0];
  
  //_varSet->computePhysicalData(*currState, _physicalData);

  //const vector<GeometricEntity*>& faces = *element->getNeighborGeos();
  //const vector<Node*>& nodes = *element->getNodes();
  //const CFuint nbNodesInElem = nodes.size();

  //cf_assert(faces.size() == nbNodesInElem);
  //const CFuint n0 = nodes[0]->getLocalID();
  
  //const CFreal Bx_n0 = nstates[n0][4];
  //const CFreal By_n0 = nstates[n0][5];
  
  // compute the gradients by applying Green Gauss in the
  // cell volume
  //CFreal divB = 0.0;
  //for (CFuint i = 0; i < nbNodesInElem; ++i) {
    // get the face normal
    //const CFuint faceID = faces[i]->getID();
    //const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
    //CFreal nx = normals[startID];
    //CFreal ny = normals[startID + 1];
    //if (static_cast<CFuint>(isOutward[faceID]) != elementID) {
      //nx *= -1.;
      //ny *= -1.;
    //}
    
    //const CFuint n1 = nodes[i]->getLocalID();
    //const CFreal Bx_n1 = nstates[n1][4];
    //const CFreal By_n1 = nstates[n1][5];
    //if (i < (nbNodesInElem - 1)) {
      //const CFuint n2 = nodes[i+1]->getLocalID();
      //const CFreal Bx_n2 = nstates[n2][4];
      //const CFreal By_n2 = nstates[n2][5];
      //divB += nx*(Bx_n1 + Bx_n2) + ny*(By_n1 + By_n2);
      //divB += nx*(nstates[n1][4] + nstates[n2][4]) + ny*(nstates[n1][5] + nstates[n2][5]);
    //}
    //else {
      //divB += nx*(Bx_n1 + Bx_n0) + ny*(By_n1 + By_n0);
      //divB += nx*(nstates[n1][4] + nstates[n0][4]) + ny*(nstates[n1][5] + nstates[n0][5]);
    //}
  //}
  //divB *= 0.5/volumes[elementID];

  //if (elementID == (nbStates-1)) {
//     _divBMax = 0.0;
//     _divBMin = 0.0;
//     for (CFuint iState = 0; iState < nbStates; ++iState) {
//       _divBMax = max(_divBMax,_divB[iState]);
//       _divBMin = min(_divBMin,_divB[iState]);
//     }
//     writeOutputFile();
//   }
    
  // source[0] = 0.0;
//   source[1] = -_divB[elementID]*_physicalData[MHDTerm::BX]*volumes[elementID];
//   source[2] = -_divB[elementID]*_physicalData[MHDTerm::BY]*volumes[elementID];
//   source[3] = -_divB[elementID]*_physicalData[MHDTerm::BZ]*volumes[elementID];
//   source[4] = -_divB[elementID]*_physicalData[MHDTerm::VX]*volumes[elementID];
//   source[5] = -_divB[elementID]*_physicalData[MHDTerm::VY]*volumes[elementID];
//   source[6] = -_divB[elementID]*_physicalData[MHDTerm::VZ]*volumes[elementID];
//   source[7] = -_divB[elementID]*(_physicalData[MHDTerm::VX]*_physicalData[MHDTerm::BX] +
// 				 _physicalData[MHDTerm::VY]*_physicalData[MHDTerm::BY] +
// 				 _physicalData[MHDTerm::VZ]*_physicalData[MHDTerm::BZ])*volumes[elementID];

  // In the source term, we take Btotal in the component corresponding to the momentum
  // equations according to Powell's JCP paper Vol.154 pp.284--309, 1999

  //_BDipole = _varSet->getMagneticDipole(*currState);
  //const CFreal BxTotal = _physicalData[MHDTerm::BX] + _BDipole[0];
  //const CFreal ByTotal = _physicalData[MHDTerm::BY] + _BDipole[1];
 
  //source[0] = 0.0;
  //source[1] = -divB*BxTotal*volumes[elementID];
  //source[2] = -divB*ByTotal*volumes[elementID];
  //source[3] = -divB*_physicalData[MHDTerm::BZ]*volumes[elementID];
  //source[4] = -divB*_physicalData[MHDTerm::VX]*volumes[elementID];
  //source[5] = -divB*_physicalData[MHDTerm::VY]*volumes[elementID];
  //source[6] = -divB*_physicalData[MHDTerm::VZ]*volumes[elementID];
  //source[7] = -divB*(_physicalData[MHDTerm::VX]*_physicalData[MHDTerm::BX] +
//				 _physicalData[MHDTerm::VY]*_physicalData[MHDTerm::BY] +
//				 _physicalData[MHDTerm::VZ]*_physicalData[MHDTerm::BZ])*volumes[elementID];

//}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> > MHD2DPowellSourceTerm::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result =
        ComputeSourceTermFVMCC::providesSockets();
  result.push_back(&socket_divBCellCenter);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<SafePtr<BaseDataSocketSink> > MHD2DPowellSourceTerm::needsSockets()
{
 vector<SafePtr<BaseDataSocketSink> > result = ComputeSourceTermFVMCC::needsSockets();
 result.push_back(&socket_avgBxFace);
 result.push_back(&socket_avgByFace);
 return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
