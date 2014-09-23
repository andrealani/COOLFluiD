#include "Common/PE.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/DirPaths.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/DataHandle.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/PathAppender.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "FiniteVolumeMHD/MHD3DProjectionPowellSourceTerm.hh"

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

Environment::ObjectProvider<MHD3DProjectionPowellSourceTerm,
               ComputeSourceTerm,
               FiniteVolumeMHDModule,
               1>
mHD3DProjectionPowellSTFVMCCProvider("MHD3DProjectionPowellST");

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionPowellSourceTerm::MHD3DProjectionPowellSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  _varSet(CFNULL),
  _divB(CFNULL),
  _divBMax(0.0),
  _divBMin(0.0),
  _leftEvals(),
  _rightEvals(),
  _unitNormal(),
  _physicalData(),
  _dataLeftState(),
  _dataRightState(),
  _saveRate(),
  _nameOutputFile()
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionPowellSourceTerm::~MHD3DProjectionPowellSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

boost::filesystem::path MHD3DProjectionPowellSourceTerm::constructFilename()
{
  boost::filesystem::path fpath(_nameOutputFile);
  fpath = PathAppender::getInstance().appendAllInfo( fpath );

  CFout << "Writing divB Errors to : " << fpath.string() << "\n";

  return fpath;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPowellSourceTerm::prepareOutputFile(std::ofstream& outputFile)
{
  outputFile << "TITLE  =  divB Errors" << "\n";
  outputFile << "VARIABLES = elementID divB" << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPowellSourceTerm::writeOutputFile()
{
  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // Execute and save file if needed...
  if(!(iter % _saveRate)) {

    SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    ofstream& outputFile = fhandle->open(constructFilename());

    prepareOutputFile(outputFile);
    
    DataHandle<State*, GLOBAL> state = 
      _globalSockets.getSocketSink<State*>("states")->getDataHandle();
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

void MHD3DProjectionPowellSourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();

  _leftEvals.resize(PhysicalModelStack::getActive()->getNbEq());
  _rightEvals.resize(PhysicalModelStack::getActive()->getNbEq());
  _unitNormal.resize(PhysicalModelStack::getActive()->getDim());

  cf_assert(_varSet.isNotNull());

  _varSet->getModel()->resizePhysicalData(_physicalData);
  _varSet->getModel()->resizePhysicalData(_dataLeftState);
  _varSet->getModel()->resizePhysicalData(_dataRightState);
  _nameOutputFile =_varSet->getNameOutputFile();
  _saveRate = _varSet->getOutputFileSaveRate();

  DataHandle<State*> state = _globalSockets.getSocketSink<State*>("states")->getDataHandle();
  const CFuint nbStates = state.size();

  _divB.resize(nbStates);
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPowellSourceTerm::setVarSet(Common::SafePtr<ConvectiveVarSet> varSet)
{
  _varSet = varSet.d_castTo<MHD3DProjectionVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPowellSourceTerm::computeSource(GeometricEntity *const element,
           RealVector& source)
{
  cf_assert(_varSet.isNotNull());
  
  DataHandle<State*, GLOBAL> state = _globalSockets.getSocketSink<State*>("states")->getDataHandle();
  DataHandle<CFreal> normals = this->socket_normals.getDataHandle();
  DataHandle<CFint> isOutward =  this->socket_isOutward.getDataHandle();
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
  CFreal bnFace;
  CFreal sumFluxPhiFace = 0.0;

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal refSpeed = _varSet->getModel()->getRefSpeed();
  const CFreal refSpeedSq = refSpeed*refSpeed;

  for (CFuint iFace = 0; iFace < nbFacesInElem; ++iFace) {

    const CFuint faceID = (*faces)[iFace]->getID();
    State *const leftState = (*faces)[iFace]->getState(0);
    State *const rightState = (*faces)[iFace]->getState(1);

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

    const CFreal faceLength = sqrt(nx*nx + ny*ny + nz*nz);
    const CFreal invFaceLength = 1.0/faceLength;

    const CFreal nxUnit = nx * invFaceLength;
    const CFreal nyUnit = ny * invFaceLength;
    const CFreal nzUnit = nz * invFaceLength;

    _unitNormal[0] = nxUnit;
    _unitNormal[1] = nyUnit;
    _unitNormal[2] = nzUnit;

    _varSet->computeEigenValues(*leftState, _unitNormal, _leftEvals);
    _varSet->computeEigenValues(*rightState, _unitNormal, _rightEvals);

    _varSet->computePhysicalData(*leftState, _dataLeftState);
    _varSet->computePhysicalData(*rightState, _dataRightState);

    bnFace = refSpeedSq*(_dataLeftState[MHDTerm::BX]*nxUnit +
                         _dataLeftState[MHDTerm::BY]*nyUnit +
			 _dataLeftState[MHDTerm::BZ]*nzUnit);
    bnFace += refSpeedSq*(_dataRightState[MHDTerm::BX]*nxUnit +
                          _dataRightState[MHDTerm::BY]*nyUnit +
			  _dataRightState[MHDTerm::BZ]*nzUnit);

    const CFreal diffPhi = _dataRightState[MHDProjectionTerm::PHI] -
      _dataLeftState[MHDProjectionTerm::PHI];

    CFreal a = 0.0;
    for (CFuint i = 0; i < nbEqs; ++i) {
      CFreal a_Av = 0.5*(_leftEvals[i] + _rightEvals[i]);
      a = max(a, std::abs(a_Av));
    }

    const CFreal fluxPhiFace = 0.5*(bnFace - a*diffPhi);

    sumFluxPhiFace += fluxPhiFace*faceLength;
  }

  _divB[elementID] = sumFluxPhiFace / (refSpeedSq*volumes[elementID]);

  if (elementID == (nbStates-1)) {
    _divBMax = 0.0;
    _divBMin = 0.0;
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      _divBMax = max(_divBMax,_divB[iState]);
      _divBMin = min(_divBMin,_divB[iState]);
    }
    writeOutputFile();
  }

  for (CFuint i = 0; i < (nbEqs-1); ++i) {
    source[i] = 0.0;
  }

  const std::string correctionType = _varSet->getModel()->getCorrectionType();

  if (correctionType == "Mixed") {
    // mixed (hyperbolic and parabolic) correction
    const CFreal dissipCoeff = _varSet->getModel()->getDissipationCoefficient();
    const CFreal dissipCoeffSq = dissipCoeff*dissipCoeff;

    source[8] = -(refSpeedSq/dissipCoeffSq)*
      _physicalData[MHDProjectionTerm::PHI]*
      volumes[elementID];
  }
  else {
    // hyperbolic correction
    source[8] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
