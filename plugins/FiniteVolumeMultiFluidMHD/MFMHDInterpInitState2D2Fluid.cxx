#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"

#include "FiniteVolumeMultiFluidMHD/MFMHDInterpInitState2D2Fluid.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/DataHandle.hh"
#include "Common/BadValueException.hh"
#include "Common/FilesystemException.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MultiFluidMHD;
using namespace COOLFluiD::Physics::Maxwell;


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MFMHDInterpInitState2D2Fluid, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule>
mFMHDInterpInitState2D2FluidProvider("MFMHDInterpInitState2D2Fluid");

//////////////////////////////////////////////////////////////////////////////

void MFMHDInterpInitState2D2Fluid::defineConfigOptions(Config::OptionList& options)
{
//   options.addConfigOption< CFreal >("mass","Mass of the external object for which the escape velocity will be computed.");
//   options.addConfigOption< CFreal >("TRef","Reference temperature to compute the isothermal speed of sound for Parker's solar wind.");
//   options.addConfigOption< CFreal >("epsilon","Desired accuracy for Parker solution.");
//   options.addConfigOption< CFreal >("rMin","Radius of the inner boundary.");
//   options.addConfigOption< CFreal >("rMax","Radius of the outer boundary.");
//   options.addConfigOption< CFreal >("rSource","Radius of the source surface for the PFSS model.");
//   options.addConfigOption< CFreal >("n","Polytropic index.");
}

//////////////////////////////////////////////////////////////////////////////

MFMHDInterpInitState2D2Fluid::MFMHDInterpInitState2D2Fluid(const std::string& name) :
  InitState(name),
  _varSet(CFNULL),
  socket_IonsDens("IonsDens"),
//  socket_Almimg("Almimg"),
//  socket_Blmreal("Blmreal"),
//  socket_Blmimg("Blmimg"),
  _nameDataFile()
{
  addConfigOptionsTo(this);

//  // the mass of the Sun is assigned by default
//  _mass = 1.98892e30;
//  setParameter("mass",&_mass);

//  _TRef = 1.5e6;
//  setParameter("TRef",&_TRef);

//  _epsilon = 1.0e-6;
//  setParameter("epsilon",&_epsilon);

//  _rMin = 1.03;
//  setParameter("rMin",&_rMin);

//  _rMax = 20.0;
//  setParameter("rMax",&_rMax);

//  _rSource = 2.5;
//  setParameter("rSource",&_rSource);

//  _n = 1.05;
//  setParameter("n",&_n);

}

//////////////////////////////////////////////////////////////////////////////

MFMHDInterpInitState2D2Fluid::~MFMHDInterpInitState2D2Fluid()
{
}

//////////////////////////////////////////////////////////////////////////////

void MFMHDInterpInitState2D2Fluid::readInputFile()
{
  DataHandle<std::vector<CFreal> > IonsDens  = socket_IonsDens.getDataHandle();
//  DataHandle<std::vector<CFreal> > Almimg  = socket_Almimg.getDataHandle();
//  DataHandle<std::vector<CFreal> > Blmreal  = socket_Blmreal.getDataHandle();
//  DataHandle<std::vector<CFreal> > Blmimg  = socket_Blmimg.getDataHandle();

  SelfRegistPtr<Environment::FileHandlerInput> fhandle =
     Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& inputFile = fhandle->open(constructFilename(_nameDataFile));

//  std::string line, tmpstr;
//  vector<CFreal> tmprow;



//  const CFuint nbLModes = _varSet->getNbLModes();

//  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
//      getline(inputFile,line);
//      istringstream liness( line );
  
//      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
//       	getline(liness, tmpstr, ' ');
//       	tmprow.push_back(atof(tmpstr.c_str()));
//      }
//      Almreal[iLine] = tmprow;
//      tmprow.clear();
//  }

//  // skip the undesired modes and the blank line separating the coefficients
//  while (getline(inputFile,line)) {
//        if (line == "") {
//           break;
//        }
//  }

//  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
//      getline(inputFile,line);
//      istringstream liness( line );

//      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
//        getline(liness, tmpstr, ' ');
//        tmprow.push_back(atof(tmpstr.c_str()));
//      }
//      Almimg[iLine] = tmprow;
//      tmprow.clear();
//  }

//  // skip the undesired modes and the blank line separating the coefficients
//  while (getline(inputFile,line)) {
//        if (line == "") {
//           break;
//        }
//  }

//  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
//      getline(inputFile,line);
//      istringstream liness( line );

//      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
//        getline(liness, tmpstr, ' ');
//        tmprow.push_back(atof(tmpstr.c_str()));
//      }
//      Blmreal[iLine] = tmprow;
//      tmprow.clear();
//  }

//  // skip the undesired modes and the blank line separating the coefficients
//  while (getline(inputFile,line)) {
//        if (line == "") {
//           break;
//        }
//  }

//  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
//      getline(inputFile,line);
//      istringstream liness( line );

//      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
//        getline(liness, tmpstr, ' ');
//        tmprow.push_back(atof(tmpstr.c_str()));
//      }
//      Blmimg[iLine] = tmprow;
//      tmprow.clear();
//  }

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void MFMHDInterpInitState2D2Fluid::executeOnTrs()
{
//  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
//  CFLogDebugMax( "InitState::executeOnTrs() called for TRS: "
//  << trs->getName() << "\n");

//  if (trs->getName() != "InnerFaces") {
//    throw BadValueException (FromHere(),"InitState not applied to InnerFaces!!!");
//  }

//  // this cannot be used for FV boundary faces because
//  // ghost state and inner state could have the same local ID !!!
//  SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();
//  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

//  std::vector<CFuint>::iterator itd;
//  CFreal rhoParker, pParker;
//  RealVector stateCoordsCartesian(PhysicalModelStack::getActive()->getDim()),
//        stateCoordsSpherical(PhysicalModelStack::getActive()->getDim()),
//        BPFSSCartesian(PhysicalModelStack::getActive()->getDim()),
//        velParkerSpherical(PhysicalModelStack::getActive()->getDim()),
//        velParkerCartesian(PhysicalModelStack::getActive()->getDim());
//  RealMatrix carSphTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim()),
//        sphCarTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim());
//  const CFreal gamma = _varSet->getModel()->getGamma();
//  /*const CFreal lRef = _varSet->getModel()->getLRef();
//  const CFreal BRef = _varSet->getModel()->getBRef();*/
//  if(_inputAdimensionalValues)
//  {
//    for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
//      State* const currState = states[(*itd)];
//      stateCoordsCartesian = currState->getCoordinates();
//      //stateCoordsCartesian *= lRef;
//      _varSet->setTransformationMatrices(stateCoordsCartesian,stateCoordsSpherical,carSphTransMat,sphCarTransMat);
//      computePFSSMagneticField(stateCoordsSpherical,BPFSSCartesian,sphCarTransMat);
//      computeParkerSolution(stateCoordsSpherical[0]/*/lRef*/,velParkerSpherical,rhoParker,pParker);
//      velParkerCartesian = sphCarTransMat*velParkerSpherical;
//      //_vFunction.evaluate(currState->getCoordinates(), *_input);
//      // input and update variables are assumed to be conservative
//      (*_input)[0] = rhoParker;
//      (*_input)[1] = rhoParker*velParkerCartesian[0];
//      (*_input)[2] = rhoParker*velParkerCartesian[1];
//      (*_input)[3] = rhoParker*velParkerCartesian[2];
//      (*_input)[4] = BPFSSCartesian[0]/*/BRef*/;
//      (*_input)[5] = BPFSSCartesian[1]/*/BRef*/;
//      (*_input)[6] = BPFSSCartesian[2]/*/BRef*/;
//      (*_input)[7] = pParker/(gamma-1.0)+(0.5*((*_input)[1]*(*_input)[1]+(*_input)[2]*(*_input)[2]+
//		(*_input)[3]*(*_input)[3])/(*_input)[0])+0.5*((*_input)[4]*(*_input)[4]+(*_input)[5]*(*_input)[5]+
//		(*_input)[6]*(*_input)[6]);
//      *currState = *_inputToUpdateVar->transform(_input);
//    }
//  }
//  else
//  {
//    State dimState;
//    for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
//      State* const currState = states[(*itd)];
//      stateCoordsCartesian = currState->getCoordinates();
//      //stateCoordsCartesian *= lRef;
//      _varSet->setTransformationMatrices(stateCoordsCartesian,stateCoordsSpherical,carSphTransMat,sphCarTransMat);
//      computePFSSMagneticField(stateCoordsSpherical,BPFSSCartesian,sphCarTransMat);
//      computeParkerSolution(stateCoordsSpherical[0]/*/lRef*/,velParkerSpherical,rhoParker,pParker);
//      velParkerCartesian = sphCarTransMat*velParkerSpherical;
//      //_vFunction.evaluate(currState->getCoordinates(), *_input);
//      // input and update variables are assumed to be conservative
//      (*_input)[0] = rhoParker;
//      (*_input)[1] = rhoParker*velParkerCartesian[0];
//      (*_input)[2] = rhoParker*velParkerCartesian[1];
//      (*_input)[3] = rhoParker*velParkerCartesian[2];
//      (*_input)[4] = BPFSSCartesian[0]/*/BRef*/;
//      (*_input)[5] = BPFSSCartesian[1]/*/BRef*/;
//      (*_input)[6] = BPFSSCartesian[2]/*/BRef*/;
//      (*_input)[7] = pParker/(gamma-1.0)+(0.5*((*_input)[1]*(*_input)[1]+(*_input)[2]*(*_input)[2]+
//                (*_input)[3]*(*_input)[3])/(*_input)[0])+0.5*((*_input)[4]*(*_input)[4]+(*_input)[5]*(*_input)[5]+
//                (*_input)[6]*(*_input)[6]);
//      dimState = *_inputToUpdateVar->transform(_input);
//      _varSet->setAdimensionalValues(dimState, *currState);
//    }
//  }
}


//////////////////////////////////////////////////////////////////////////////

void MFMHDInterpInitState2D2Fluid::setup()
{

  InitState::setup();

  Common::SafePtr<Common::ConnectivityTable<CFuint> > cells =
    MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");
  const CFuint nbCells = cells->nbRows();
  cout << "MFMHDInterpInitState2D2Fluid::setup():: nbCells = " << nbCells << endl;


  _varSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
  cf_assert(_varSet.isNotNull());

  DataHandle<std::vector<CFreal> > IonsDens  = socket_IonsDens.getDataHandle();
  IonsDens.resize(nbCells);

//  DataHandle<std::vector<CFreal> > Almimg  = socket_Almimg.getDataHandle();
//  Almimg.resize(_varSet->getNbLModes()+1);

//  DataHandle<std::vector<CFreal> > Blmreal  = socket_Blmreal.getDataHandle();
//  DataHandle<std::vector<CFreal> > Blmimg  = socket_Blmimg.getDataHandle();
//  Blmreal.resize(_varSet->getNbLModes()+1);
//  Blmimg.resize(_varSet->getNbLModes()+1);

//  _namePFSSDataFile = _varSet->getNamePFSSCoeffFile();

//  // if this is a parallel simulation, only ONE process at a time reads the file
//  if (PE::GetPE().IsParallel()) {
    
//    PE::GetPE().setBarrier();

//    for (CFuint i = 0; i < PE::GetPE().GetProcessorCount(); ++i) {
//      if (i == PE::GetPE().GetRank ()) {
//        readInputFile();
//      }
      
//      PE::GetPE().setBarrier();
//    }
//  }
//  else {
//    readInputFile();
//  }
}

//////////////////////////////////////////////////////////////////////////////

boost::filesystem::path MFMHDInterpInitState2D2Fluid::constructFilename(std::string fileName)
{
//  boost::filesystem::path fpath(fileName);

//  CFout << "Reading PFSS spherical harmonics coefficients from: " << fpath.string() << "\n";

//  return fpath;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
MFMHDInterpInitState2D2Fluid::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result =
          InitState::providesSockets();
  result.push_back(&socket_IonsDens);
//  result.push_back(&socket_Almimg);
//  result.push_back(&socket_Blmreal);
//  result.push_back(&socket_Blmimg);
//  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
