#include <numeric>
#include <algorithm>

#include "StrongSubInletLinEuler2DCons.hh"
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

MethodCommandProvider<StrongSubInletLinEuler2DCons, FluctuationSplitData, FluctSplitLinEulerModule> strongSubInletLinEuler2DConsProvider("StrongSubInletLinEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

StrongSubInletLinEuler2DCons::StrongSubInletLinEuler2DCons(const std::string& name) :
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

StrongSubInletLinEuler2DCons::~StrongSubInletLinEuler2DCons()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSubInletLinEuler2DCons::needsSockets()
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

void StrongSubInletLinEuler2DCons::setup()
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
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isBState = socket_isBState.getDataHandle();
  DataHandle< InwardNormalsData*> m_normals = socket_normals.getDataHandle();

//   prepares to loop over cells by getting the GeometricEntityPool
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = MeshDataStack::getActive()->getTrs("InnerCells");

  Common::SafePtr< vector<CFuint> > statesIdx = getCurrentTRS()->getStatesInTrs();
  
  SafePtr<TopologicalRegionSet> bndTrs = getCurrentTRS();
  const CFuint bndNbNodes= bndTrs->getNbNodesInTrs();
  bndNodNorm.resize(bndNbNodes);

  for(CFuint i=0; i<bndNbNodes; i++)
    bndNodNorm[i].resize(2);
// go through the nodes

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];

// get the boundary elements containing the state
    vector <CFuint> cells = bndNod2Elm[iState];
    CFuint nbElemsOfState = cells.size();

    RealVector NodalNormal(2);
    NodalNormal = 0.0;

// loop over the elements
    for (CFuint elem=0; elem<nbElemsOfState; ++elem) {

        // build the GeometricEntity
        CFuint cellID = cells[elem];
        geoData.idx = cellID;
        GeometricEntity& cell = *geoBuilder->buildGE();
        vector<State*> *const statesInCell = cell.getStates();
        const CFuint nbStatesInCell = statesInCell->size();

        CFuint NrBoundNodes = 0;

        for (CFuint elemState=0; elemState<nbStatesInCell; ++elemState) {
           State *const currState = (*statesInCell)[elemState];
           CFuint currStateID = currState->getLocalID();

           if(isBState[currStateID])
             NrBoundNodes+= 1;
        }

        if (NrBoundNodes==2) {
          for (CFuint elemState=0; elemState<nbStatesInCell; ++elemState) {
            State *const currState = (*statesInCell)[elemState];
            CFuint currStateID = currState->getLocalID();

            if(!isBState[currStateID]) {
               NodalNormal[0] += m_normals[cellID]->getNodalNormComp(elemState,0);
               NodalNormal[1] += m_normals[cellID]->getNodalNormComp(elemState,1);
            }
          }
        }
        else if (NrBoundNodes==3) {
          for (CFuint elemState=0; elemState<nbStatesInCell; ++elemState) {
            State *const currState = (*statesInCell)[elemState];
            CFuint currStateID = currState->getLocalID();
	    
	    if (currStateID == stateID) {
	      
// 	      cout << "Here we are...stateID: "  << stateID << "\n";
               NodalNormal[0] = m_normals[cellID]->getNodalNormComp(elemState,0);
               NodalNormal[1] = m_normals[cellID]->getNodalNormComp(elemState,1);
            }
          }  
	  
	}

      geoBuilder->releaseGE();

      }

   CFreal Nlength = 0.0;
   for (CFuint dim=0; dim<2; dim++)
    Nlength += NodalNormal[dim]*NodalNormal[dim];
   Nlength = sqrt(Nlength);


   NodalNormal/=-Nlength;

   bndNodNorm[iState] = NodalNormal;

  }


/// for the computation of the distribution matrix

  getMethodData().getDistributionData().computeBetas = true;

  CFuint m_maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();

  m_kPlus.resize(m_maxNbStatesInCell);
  m_k.resize(m_maxNbStatesInCell);
  m_kPast.resize(m_maxNbStatesInCell);
  m_Statenormals.resize(m_maxNbStatesInCell);

  for (CFuint i = 0; i < m_maxNbStatesInCell; ++i) {
    m_Statenormals[i] = new RealVector(2);
  }

  m_adimNormal.resize(2);
  m_adimCharNormal.resize(2);

}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletLinEuler2DCons::executeOnTrs()
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

      RealVector bcNormalState =  bndNodNorm[iState];

      CFreal ncharx = (bcNormalState)[0];
      CFreal nchary = (bcNormalState)[1];

// this is entropy, good as it is.
      charRes[0] = 0.0;

// this is vorticity in the direction where acoustic wave is propagating, so they are decpoupled

      charRes[1] = 0.0;
      charRes[2] = -(rhs(stateID, 1, m_nbEqs)*ncharx + rhs(stateID, 2, m_nbEqs)*nchary + oneoverc*rhs(stateID, 3, m_nbEqs));;
      charRes[3] = 0.0;
      
// transform back to conservative
      Res[0] = (charRes[0]+0.5/c*(charRes[2]));
      Res[1] = (nchary*charRes[1]+0.5*ncharx*(charRes[2]));
      Res[2] = (-ncharx*charRes[1]+0.5*nchary*(charRes[2]));
      Res[3] = (0.5*c*(charRes[2]));

     for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq) {
       rhs(stateID, iEq, m_nbEqs) = -Res[iEq];
     }


    }

    isUpdated[stateID] = true; // flagging is important!!!!!
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletLinEuler2DCons::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "LinEuler2DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::LinearizedEuler::LinEuler2DCons>());

  cf_assert(_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletLinEuler2DCons::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
