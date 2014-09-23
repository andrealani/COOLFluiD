#include <numeric>
#include <algorithm>

#include "StrongSubInletLinEuler3DCons.hh"
#include "FluctSplit/CreateBoundaryNodalNormals.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/FluctSplitLinEuler.hh"
#include "Framework/MeshData.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "MathTools/MatrixInverter.hh"
#include "Environment/ObjectProvider.hh"
#include "FluctSplit/SpaceTime_Splitter.hh"
#include "Framework/CFL.hh"
#include "MathTools/RealMatrix.hh"
#include "Framework/SubSystemStatus.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::LinearizedEuler;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongSubInletLinEuler3DCons, FluctuationSplitData, FluctSplitLinEulerModule> strongSubInletLinEuler3DConsProvider("StrongSubInletLinEuler3DCons");

//////////////////////////////////////////////////////////////////////////////

StrongSubInletLinEuler3DCons::StrongSubInletLinEuler3DCons(const std::string& name) :
  FluctuationSplitCom(name),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_nodes("nodes"),
  socket_isUpdated("isUpdated"),
  socket_updateCoeff("updateCoeff"),
  socket_volumes("volumes"),
  socket_pastStates("pastStates"),
  socket_isBState("isBState"),
  m_Statenormals(0),
  m_kPast(0),
  _varSet(),
  bndNod2Elm(0),
  bndNodNorm(0),
  m_adimNormal(),
  m_adimCharNormal(),
  m_nbEqs(),
  m_kPlus(0),
  m_k(0)
{
}

//////////////////////////////////////////////////////////////////////////////

StrongSubInletLinEuler3DCons::~StrongSubInletLinEuler3DCons()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSubInletLinEuler3DCons::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_nodes);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_volumes);

  result.push_back(&socket_pastStates);
  result.push_back(&socket_isBState);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletLinEuler3DCons::setup()
{
  FluctuationSplitCom::setup();

  m_nbEqs = PhysicalModelStack::getActive()->getNbEq();

// create boundary nodal normals, which pointing outwards
  _bcNormals.resize(getTrsList().size());
  _varSet->setup();

  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
  obj.setDataSockets(socket_normals,socket_faceNeighCell);
  obj.create(getTrsList(), _bcNormals);

/// create the boundary normals for the nodes
//   DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
//   DataHandle<bool> isBState = socket_isBState.getDataHandle();
//   DataHandle< InwardNormalsData*> m_normals = socket_normals.getDataHandle();
// 
// //   prepares to loop over cells by getting the GeometricEntityPool
//   Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
//     geoBuilder = getMethodData().getStdTrsGeoBuilder();
//   StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
//   geoData.trs = MeshDataStack::getActive()->getTrs("InnerCells");
// 
//   Common::SafePtr< vector<CFuint> > statesIdx = getCurrentTRS()->getStatesInTrs();
//   
//   SafePtr<TopologicalRegionSet> bndTrs = getCurrentTRS();
//   const CFuint bndNbNodes= bndTrs->getNbNodesInTrs();
//   bndNodNorm.resize(bndNbNodes);
// 
//   for(CFuint i=0; i<bndNbNodes; i++)
//     bndNodNorm[i].resize(3);
// // go through the nodes
// 
//   for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
//     const CFuint stateID = (*statesIdx)[iState];
// 
// // get the boundary elements containing the state
//     vector <CFuint> cells = bndNod2Elm[iState];
//     CFuint nbElemsOfState = cells.size();
// 
//     RealVector NodalNormal(3);
//     NodalNormal = 0.0;
// 
// // loop over the elements
//     for (CFuint elem=0; elem<nbElemsOfState; ++elem) {
// 
//         // build the GeometricEntity
//         CFuint cellID = cells[elem];
//         geoData.idx = cellID;
//         GeometricEntity& cell = *geoBuilder->buildGE();
//         vector<State*> *const statesInCell = cell.getStates();
//         const CFuint nbStatesInCell = statesInCell->size();
// 
//         CFuint NrBoundNodes = 0;
// 
//         for (CFuint elemState=0; elemState<nbStatesInCell; ++elemState) {
//            State *const currState = (*statesInCell)[elemState];
//            CFuint currStateID = currState->getLocalID();
// 
//            if(isBState[currStateID])
//              NrBoundNodes+= 1;
//         }
// 
//         if (NrBoundNodes==2) {
//           for (CFuint elemState=0; elemState<nbStatesInCell; ++elemState) {
//             State *const currState = (*statesInCell)[elemState];
//             CFuint currStateID = currState->getLocalID();
// 
//             if(!isBState[currStateID]) {
//                NodalNormal[0] += m_normals[cellID]->getNodalNormComp(elemState,0);
//                NodalNormal[1] += m_normals[cellID]->getNodalNormComp(elemState,1);
// 
//             }
//           }
//         }
// 
//       geoBuilder->releaseGE();
// 
//       }
/*
   CFreal Nlength = 0.0;
   for (CFuint dim=0; dim<2; dim++)
    Nlength += NodalNormal[dim]*NodalNormal[dim];
   Nlength = sqrt(Nlength);


   NodalNormal/=-Nlength;

   bndNodNorm[iState] = NodalNormal;*/

// cout << bndNodNorm[iState] << "\t" << states[stateID]->getCoordinates() <<"\n";

//   }


/// for the computation of the distribution matrix

  getMethodData().getDistributionData().computeBetas = true;

  CFuint m_maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  m_kPlus.resize(m_maxNbStatesInCell);
  m_k.resize(m_maxNbStatesInCell);
  m_kPast.resize(m_maxNbStatesInCell);
  m_Statenormals.resize(m_maxNbStatesInCell);

  for (CFuint i = 0; i < m_maxNbStatesInCell; ++i) {
    m_Statenormals[i] = new RealVector(3);
  }

  m_adimNormal.resize(3);
  m_adimCharNormal.resize(3);

}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletLinEuler3DCons::executeOnTrs()
{
  vector<RealVector>* bcNormalsInTrs = &(_bcNormals[getCurrentTrsID()]);

  DataHandle<CFreal> m_volumes = socket_volumes.getDataHandle();
  DataHandle< InwardNormalsData*> m_normals = socket_normals.getDataHandle();

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<bool> isBState = socket_isBState.getDataHandle();

  DataHandle<State*> pastStatesStorage = socket_pastStates.getDataHandle();

  // prepares to loop over cells by getting the GeometricEntityPool
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = MeshDataStack::getActive()->getTrs("InnerCells");

  Common::SafePtr< vector<CFuint> > statesIdx = getCurrentTRS()->getStatesInTrs();

  const CFreal dt = SubSystemStatusStack::getActive()->getDT();

// go through all the states involved in the boundary trs
   for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];

// see if the state is updated already
    if (!isUpdated[stateID]) {
      const CFreal updateCoeffValue = updateCoeff[stateID];

      RealVector Res(m_nbEqs);
      RealVector charRes(m_nbEqs);
      for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq) {
        Res[iEq] = 0.;
        charRes[iEq] = 0.;
      }

      const RealVector& linearData = _varSet->getModel()->getPhysicalData();
      const CFreal c     = linearData[LinEulerTerm::c];
      const CFreal oneoverc = 1./c;

      if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps()) {

// the right vectors are not computed yet, MUST BE FINISHED      
      
      RealVector bcNormalState = (*bcNormalsInTrs)[stateID];
      bcNormalState *= -1.;

      const CFreal ncharx = (bcNormalState)[0];
      const CFreal nchary = (bcNormalState)[1];
      const CFreal ncharz = (bcNormalState)[2];
      
//       RealVector bcNormalState =  bndNodNorm[iState];
      
//       CFreal ncharx = (bcNormalState)[0];
//       CFreal nchary = (bcNormalState)[1];
//       CFreal ncharz = (bcNormalState)[2];


      charRes[0] = 0.0;
      charRes[1] = 0.0;
      charRes[2] = 0.0;
      charRes[3] = -(rhs(stateID, 1, m_nbEqs)*ncharx + rhs(stateID, 2, m_nbEqs)*nchary + rhs(stateID, 3, m_nbEqs)*ncharz + oneoverc*rhs(stateID, 4, m_nbEqs));
      charRes[4] = 0.0;
      
// transform back to conservative
      Res[0] = 0.5/c*(charRes[3]+charRes[4]);
      // just becouse the vorticity waves are zero
      Res[1] = 0.5*ncharx*(charRes[3]-charRes[4]);
      Res[2] = 0.5*nchary*(charRes[3]-charRes[4]);
      Res[3] = 0.5*ncharz*(charRes[3]-charRes[4]);
      Res[4] = (0.5*c*(charRes[3]+charRes[4]));      

     for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq) {
       rhs(stateID, iEq, m_nbEqs) = -Res[iEq];
     }

    }

    isUpdated[stateID] = true; // flagging is important!!!!!
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletLinEuler3DCons::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "LinEuler3DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::LinearizedEuler::LinEuler3DCons>());

  cf_assert(_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletLinEuler3DCons::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
