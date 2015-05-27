#include "Framework/MeshData.hh"
#include "Common/ShouldNotBeHereException.hh"
#include "Environment/Factory.hh"

#include "FluctSplit/WeakBC.hh"
#include "FluctSplit/InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

void WeakBC::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("alpha","distribution coefficient");
  options.addConfigOption< CFuint >("FaceGeoOrder","HO boundary face - order of geometry");
  options.addConfigOption< CFuint >("SolOrder","HO boundary face - order of solution");
  options.addConfigOption<std::string>("StateInterpolator", "State interpolator for BCs");
}
    
//////////////////////////////////////////////////////////////////////////////

WeakBC::WeakBC(const std::string& name) :
  FluctuationSplitCom(name),
  m_sInterpolator(),
  socket_rhs("rhs"),
  socket_normals("normals"),
  socket_updateCoeff("updateCoeff"),
  socket_faceNeighCell("faceNeighCell"),
  socket_states("states"),
  m_solutionToDistMatTrans(CFNULL),
  m_distToSolutionMatTrans(CFNULL),
  m_linearToDistMatTrans(CFNULL),
  m_updateToLinearVecTrans(CFNULL),
  m_linearizer(CFNULL),
  m_distribVar(CFNULL),
  m_faceNormal(),
  m_adimNormal(),
  m_rightEv(),
  m_leftEv(),
  m_eValues(),
  m_eValuesP(),
  m_kPlus(),
  m_tFlux(),
  m_twoStates(2),
  m_flagState(),
  m_faceGeoOrder(1)
{
  addConfigOptionsTo(this);
  
  m_sInterpolatorStr = "Null";
  this->setParameter("StateInterpolator",&m_sInterpolatorStr); 
  
  m_alpha = 1.0;
  setParameter("alpha",&m_alpha);

  setParameter("FaceGeoOrder",&m_faceGeoOrder);

  m_solorder = 1;
  setParameter("SolOrder",&m_solorder);
}

//////////////////////////////////////////////////////////////////////////////

WeakBC::~WeakBC()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WeakBC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_normals);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC::unsetup()
{
  FluctuationSplitCom::unsetup();
  
  cout << "WeakBC::unsetup()" << endl; exit(1);
  m_sInterpolator->unsetup();
  
  // allocating data for the temporary local linearized states
  for (CFuint i = 0; i < MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell(); ++i) {
    deletePtr(m_zStates[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC::setup()
{
  FluctuationSplitCom::setup();
  
  // set the interpolator
  m_sInterpolator->setup();
  
  m_faceNormal.resize(PhysicalModelStack::getActive()->getDim());
  m_adimNormal.resize(PhysicalModelStack::getActive()->getDim());
  m_rightEv.resize(PhysicalModelStack::getActive()->getNbEq(), PhysicalModelStack::getActive()->getNbEq());
  m_leftEv.resize(PhysicalModelStack::getActive()->getNbEq(), PhysicalModelStack::getActive()->getNbEq());
  m_eValues.resize(PhysicalModelStack::getActive()->getNbEq());
  m_eValuesP.resize(PhysicalModelStack::getActive()->getNbEq());
  m_kPlus.resize(PhysicalModelStack::getActive()->getNbEq(), PhysicalModelStack::getActive()->getNbEq());
  m_tFlux.resize(PhysicalModelStack::getActive()->getNbEq());

  ///HO Hack
//   m_hoFaceNormal0.resize(PhysicalModelStack::getActive()->getDim());
//   m_hoFaceNormal1.resize(PhysicalModelStack::getActive()->getDim());
//   m_hoFaceNormal2.resize(PhysicalModelStack::getActive()->getDim());

  // get the distribution var set
  m_distribVar = getMethodData().getDistribVar();

  // get the linearizer
  m_linearizer = getMethodData().getLinearizer();

  m_solutionToDistMatTrans = getMethodData().getSolutionToDistribMatTrans();
  m_distToSolutionMatTrans = getMethodData().getDistribToSolutionMatTrans();
  m_linearToDistMatTrans = getMethodData().getLinearToDistribMatTrans();
  m_updateToLinearVecTrans = getMethodData().getUpdateToLinearVecTrans();

  // get the data handle for the states
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  m_flagState.resize(states.size());
  m_flagState = false;

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  m_zStates.resize(maxNbStatesInCell);
  m_tStates.resize(maxNbStatesInCell);

  // Resizing tStates
  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    m_tStates[i].resize(PhysicalModelStack::getActive()->getNbEq());
  }

  // allocating data for the temporary local linearized states
  for (CFuint i = 0; i < MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell(); ++i) {
    m_zStates[i]  = new RealVector(PhysicalModelStack::getActive()->getNbEq());
    *m_zStates[i] = 0.0;
  }

 }

//////////////////////////////////////////////////////////////////////////////

void WeakBC::setFaceNormal(const CFuint faceID)
{
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  const CFuint cellTrsID = faceNeighCell[faceID].first;
  const CFuint iFaceLocal = faceNeighCell[faceID].second;
  Common::SafePtr<TopologicalRegionSet> cellTrs = faceNeighCell[faceID].third;
  const CFuint cellLocalID = cellTrs->getLocalGeoID(cellTrsID);

  const CFuint dim = PhysicalModelStack::getActive()->getDim();

//   if(m_faceGeoOrder == 2)
//   {
// 
//     CFuint hoFace0, hoFace1;
// 
//     switch(iFaceLocal) {
//     case 0: {
//               hoFace0 = 2;
//               hoFace1 = 5;
//             }
//             break;
//     case 1: {
//               hoFace0 = 3;
//               hoFace1 = 6;
//             }
//             break;
//     case 2: {
//               hoFace0 = 7;
//               hoFace1 = 1;
//             }
//             break;
//     default:
//       throw Common::ShouldNotBeHereException(FromHere(),"iFaceLocal != 0 && iFaceLocal != 1 && iFaceLocal != 2");
//     }



///Another try (begin)

//       for (CFuint iDim = 0; iDim < dim; ++iDim) {
//       // external normal is already inward for ghost state
//       m_hoFaceNormal0[iDim] = normals[cellLocalID]->getFaceNormComp(hoFace0,iDim);
//       m_hoFaceNormal1[iDim] = normals[cellLocalID]->getFaceNormComp(hoFace1,iDim);
// 
//       }
// 
//       J0 = m_hoFaceNormal0.norm2();
//       J1 = m_hoFaceNormal1.norm2();
// 
//       totFaceLen = J0 + J1;
// 
//       m_hoFaceNormal0 = 1.0/J0 * m_hoFaceNormal0;
//       m_hoFaceNormal1 = 1.0/J1 * m_hoFaceNormal1;
// 
//       for (CFuint iDim = 0; iDim < dim; ++iDim) {
//       m_faceNormal[iDim] = 0.5*(m_hoFaceNormal0[iDim] + m_hoFaceNormal1[iDim]);
//       }
// 
// 
//       m_faceNormal *= -totFaceLen;

///Another try (end)



//       for (CFuint iDim = 0; iDim < dim; ++iDim) {
//       // external normal is already inward for ghost state
//       m_faceNormal[iDim] = -normals[cellLocalID]->getFaceNormComp(hoFace0,iDim) - normals[cellLocalID]->getFaceNormComp(hoFace1,iDim);
//       }

  ///Just a try: normal at arbitrary point
//   Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
//   geoBuilder = getMethodData().getStdTrsGeoBuilder();
//
//   StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
//   geoData.trs = getCurrentTRS();
//
//   geoData.idx = faceID;
//   GeometricEntity& face = *geoBuilder->buildGE();
//
//   m_nodes = *face.getNodes();
//
//   m_CP2N.ComputeBNormal(m_nodes,0.5,m_faceNormal);



//   } //if geometry order == 2

// else if(m_faceGeoOrder == 3)
//   {
// 
//     CFuint hoFace0, hoFace1, hoFace2;
// 
//     switch(iFaceLocal) {
//     case 0: {
//               hoFace0 = 2;
//               hoFace1 = 17;
//               hoFace2 = 5;
//             }
//             break;
//     case 1: {
//               hoFace0 = 3;
//               hoFace1 = 9;
//               hoFace2 = 6;
//             }
//             break;
//     case 2: {
//               hoFace0 = 7;
//               hoFace1 = 13;
//               hoFace2 = 1;
//             }
//             break;
//     default:
//       throw Common::ShouldNotBeHereException(FromHere(),"iFaceLocal != 0 && iFaceLocal != 1 && iFaceLocal != 2");
//     }
// 
// 
//       for (CFuint iDim = 0; iDim < dim; ++iDim) {
//       // external normal is already inward for ghost state
//       m_hoFaceNormal0[iDim] = normals[cellLocalID]->getFaceNormComp(hoFace0,iDim);
//       m_hoFaceNormal1[iDim] = normals[cellLocalID]->getFaceNormComp(hoFace1,iDim);
//       m_hoFaceNormal2[iDim] = normals[cellLocalID]->getFaceNormComp(hoFace2,iDim);
// 
//       }
// 
//       J0 = m_hoFaceNormal0.norm2();
//       J1 = m_hoFaceNormal1.norm2();
//       J2 = m_hoFaceNormal2.norm2();
// 
//       totFaceLen = J0 + J1 + J2;
// 
//       m_hoFaceNormal0 = 1.0/J0 * m_hoFaceNormal0;
//       m_hoFaceNormal1 = 1.0/J1 * m_hoFaceNormal1;
//       m_hoFaceNormal2 = 1.0/J2 * m_hoFaceNormal2;
// 
//       for (CFuint iDim = 0; iDim < dim; ++iDim) {
//       m_faceNormal[iDim] = 1.0/3.0*(m_hoFaceNormal0[iDim] + m_hoFaceNormal1[iDim] + m_hoFaceNormal2[iDim]);
//       }
// 
// 
//       m_faceNormal *= -totFaceLen;
// 
// 
//  }  ///geometry order == 3





//   else
  for (CFuint iDim = 0; iDim < dim; ++iDim) {
    // external normal is already inward for ghost state
    m_faceNormal[iDim] = -normals[cellLocalID]
      ->getFaceNormComp(iFaceLocal,iDim);
  }

//     CFout << "Inlet normal = [" << m_faceNormal[0] << "," << m_faceNormal[1] << "]\n";
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC::computeFlux(const vector<State*>& states, RealVector& flux)
{
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  vector<State*> *const linearStates =
    m_updateToLinearVecTrans->transform(const_cast<vector<State*>*>(&states));

  /// @todo this is needed because linearStates->size() > 2
  m_twoStates[0] = (*linearStates)[0];
  m_twoStates[1] = (*linearStates)[1];

  // linearize the states in the cell
  // includes the transformation from solution to linearization
  // variables to evaluate the jacobians in the average state
  m_linearizer->linearize(m_twoStates);

  // transformation from solution to consistent variables
  vector<State*> *const tStates =
    m_linearToDistMatTrans->transformFromRef(&m_twoStates);

  const CFreal kCoeff = 1./ PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  // adimensionalize the normal
  const CFreal faceArea = m_faceNormal.norm2();
  m_adimNormal = m_faceNormal/faceArea;

  m_distribVar->computeEigenValuesVectors(m_rightEv,
             m_leftEv,
             m_eValues,
             m_adimNormal);

  for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
    m_eValuesP[iEq] = max(0.,m_eValues[iEq]);
  }
  m_kPlus = m_rightEv*(m_eValuesP*m_leftEv);
  m_kPlus *= kCoeff*faceArea;

  // the first state is a ghost state and gives no contribution
  // to the update coefficient
  updateCoeff[states[1]->getLocalID()] += kCoeff*faceArea*m_eValuesP.max();

  // flux in the transformed variables
  m_tFlux = m_kPlus*(*(*tStates)[0] - *(*tStates)[1]);

  // transform the flux back to solution variables
  // here you have copying try to avoid it ...
  flux = *m_distToSolutionMatTrans->transformFromRef(&m_tFlux);
}

//////////////////////////////////////////////////////////////////////////////

void WeakBC::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);
  
  SafePtr<StateInterpolator::PROVIDER> stateInterpProv = CFNULL;
  try {
    stateInterpProv = Factory<StateInterpolator>::getInstance().getProvider(m_sInterpolatorStr);
  }
  catch (NoSuchValueException& e) {
    CFLog(WARN, "FVMCC_BaseBC::configure() => StateInterpolator \"" << 
	  m_sInterpolatorStr << "\" missing => reset to \"Null\"\n");
    
    m_sInterpolatorStr = "Null";
    stateInterpProv = Factory<StateInterpolator>::getInstance().getProvider(m_sInterpolatorStr);
  }
  
  m_sInterpolator.reset(stateInterpProv->create(m_sInterpolatorStr));
  
  configureNested(m_sInterpolator.getPtr(), args);
}
    
//////////////////////////////////////////////////////////////////////////////

} // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
