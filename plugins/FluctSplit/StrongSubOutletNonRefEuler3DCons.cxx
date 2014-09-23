#include <numeric>
#include <algorithm>

#include "StrongSubOutletNonRefEuler3DCons.hh"
#include "FluctSplit/CreateBoundaryNodalNormals.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/FluctSplitSpaceTimeNavierStokes.hh"
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
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongSubOutletNonRefEuler3DCons, FluctuationSplitData, FluctSplitSpaceTimeNavierStokesModule> strongSubOutletNonRefEuler3DConsProvider("StrongSubOutletNonRefEuler3DCons");

//////////////////////////////////////////////////////////////////////////////

StrongSubOutletNonRefEuler3DCons::StrongSubOutletNonRefEuler3DCons(const std::string& name) :
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
  m_adimNormal(),
  m_adimCharNormal(),
  m_nbEqs(),
  m_kPlus(0),
  m_k(0)
{

}

//////////////////////////////////////////////////////////////////////////////

StrongSubOutletNonRefEuler3DCons::~StrongSubOutletNonRefEuler3DCons()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSubOutletNonRefEuler3DCons::needsSockets()
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

void StrongSubOutletNonRefEuler3DCons::setup()
{
  FluctuationSplitCom::setup();

  m_nbEqs = PhysicalModelStack::getActive()->getNbEq();
// create boundary nodal normals, which pointing outwards
  _bcNormals.resize(getTrsList().size());
  _varSet->setup();
  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
  obj.setDataSockets(socket_normals,socket_faceNeighCell);
  obj.create(getTrsList(), _bcNormals);

/// create the cell connectivity

  // find the inner trs
  vector< SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  SafePtr<TopologicalRegionSet> innerTrs=0;
  for (CFuint i = 0; i < trs.size(); ++i) {
    if (trs[i]->hasTag("inner")) {
      innerTrs = trs[i];
      break;
    }
  }
  if (innerTrs==0) Common::NoSuchValueException(FromHere(),"Trs with tag 'inner' not found.");
  // set up numbers
  DataHandle< Node*,GLOBAL > innerNodes = socket_nodes.getDataHandle();
  const CFuint innerNbNodes = innerNodes.size();
  CFVec<CFint> bndNodeGlobal2Local(innerNbNodes);
  const CFuint innerNbGeos=innerTrs->getLocalNbGeoEnts();
  const Common::SafePtr< Common::ConnectivityTable<CFuint> > innerGeo2Node=innerTrs->getGeo2NodesConn();

  // set the boundary trs
  SafePtr<TopologicalRegionSet> bndTrs = getCurrentTRS();
   // getting nodes and preliminary settings
  Common::SafePtr< std::vector< CFuint > > bndNodes = bndTrs->getNodesInTrs();
  const CFuint bndNbNodes= bndTrs->getNbNodesInTrs();
  bndNod2Elm.resize(bndNbNodes);
  bndNodeGlobal2Local=-1;
  // cycle all the nodes in the TRS to set isBndNode flag
  CFuint bndCheckNbNode=0;
  for (std::vector< CFuint >::iterator itd = bndNodes->begin(); itd != bndNodes->end(); ++itd)
     bndNodeGlobal2Local[*itd]=bndCheckNbNode++;

  // build boundary node -> innercells elem connectivity
  for (CFuint i=0; i<innerNbGeos; i++) {
    const CFuint innerGeoLocalID=innerTrs->getLocalGeoID(i);
    const CFuint innerNbGeoNodes=innerGeo2Node->nbCols(i);
    for (CFuint j=0; j<innerNbGeoNodes; j++)
      if (bndNodeGlobal2Local[(*innerGeo2Node)(i,j)]!=-1)
        bndNod2Elm[bndNodeGlobal2Local[(*innerGeo2Node)(i,j)]].push_back(innerGeoLocalID);
  }

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

  m_adimCharNormal.resize(3);
 
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletNonRefEuler3DCons::executeOnTrs()
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

      const CFreal gamma = _varSet->getModel()->getGamma();
// see if the state is updated already
    if (!isUpdated[stateID]) {
      const CFreal updateCoeffValue = updateCoeff[stateID];

      CFreal charResThree = 0.0;
      CFreal charResFour = 0.0;
      RealVector Res(m_nbEqs);
      RealVector charRes(m_nbEqs);
      for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq) {
        Res[iEq] = 0.;
        charRes[iEq] = 0.;
      }
 
      const RealVector& linearData = _varSet->getModel()->getPhysicalData();
      const CFreal c     = linearData[EulerTerm::A];
      const CFreal oneoverc = 1./c;

      if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps()) {

      // get the boundary elements containing the state
      vector <CFuint> cells = bndNod2Elm[iState];
      CFuint nbElemsOfState = cells.size();
      
      RealVector bcNormalState = (*bcNormalsInTrs)[stateID];
      bcNormalState *= -1.;

      const CFreal ncharx = (bcNormalState)[0];
      const CFreal nchary = (bcNormalState)[1];
      const CFreal ncharz = (bcNormalState)[2];
      
       
//  loop over the elements in which the state is involved
    for (CFuint elem=0; elem<nbElemsOfState; ++elem) {

        // build the GeometricEntity
        CFuint cellID = cells[elem];
        geoData.idx = cellID;
        GeometricEntity& cell = *geoBuilder->buildGE();
        vector<State*> *const statesInCell = cell.getStates();
        const CFuint nbStatesInCell = statesInCell->size();

// compute the characteristic variables and the derivatives needed
        RealVector acoustics(nbStatesInCell);
        RealVector omega(nbStatesInCell);

	CFuint boundaryState = 10000;
        for (CFuint elemState=0; elemState<nbStatesInCell; ++elemState) {
          State *const currState = (*statesInCell)[elemState];
          if(stateID == currState->getLocalID()) {
            boundaryState = elemState;
          }
        }

        for (CFuint elemState=0; elemState<nbStatesInCell; ++elemState) {
          State *const currState = (*statesInCell)[elemState];
	  CFreal p;
          p = (gamma-1.0)*((*currState)[4]-
				  0.5*((*currState)[1]*(*currState)[1]+(*currState)[2]*(*currState)[2]+(*currState)[3]*(*currState)[3])
				  /(*currState)[0]);
          acoustics[elemState]= 2.0*oneoverc*(p)/(*currState)[0];
          omega[elemState]= 0.0;

        }
        
        
        

/********************* Matrix distribution in characteristic ***********************/

        CFreal resCharElemThree=0.;
        CFreal resCharElemFour=0.;

        RealVector faceLength(nbStatesInCell);

        for (CFuint i = 0; i < nbStatesInCell; ++i) {
          (*m_Statenormals[i])[0] = m_normals[cellID]->getNodalNormComp(i,0);
          (*m_Statenormals[i])[1] = m_normals[cellID]->getNodalNormComp(i,1);
	  (*m_Statenormals[i])[2] = m_normals[cellID]->getNodalNormComp(i,2);
          faceLength[i] = m_normals[cellID]->getAreaNode(i);
        }
        CFreal Area = m_volumes[cellID];

        computeCharK(*statesInCell, m_kPlus, m_k, m_kPast, m_Statenormals, faceLength, Area);

//         compute the distribution coefficients betas in this cell
        m_betas = distributeLDA(m_kPlus, m_betas, boundaryState);

        RealVector acoustic_past(nbStatesInCell);
        RealVector omega_past(nbStatesInCell);

        for (CFuint i = 0; i < nbStatesInCell; ++i) {
          State *const currState = (*statesInCell)[i];
          CFuint stateIDlocal = currState->getLocalID();
          State const currState_past = (*pastStatesStorage[stateIDlocal]);

	  CFreal p = (gamma-1.0)*((currState_past)[4]-
				  0.5*((currState_past)[1]*(currState_past)[1]+(currState_past)[2]*(currState_past)[2]+(currState_past)[3]*(currState_past)[3])
				  /(currState_past)[0]);
          acoustic_past[i]= 2.0 * oneoverc*(p)/(currState_past)[0];
          omega_past[i]= 0.0;
        }
        
        for (CFuint i = 0; i < nbStatesInCell; ++i) {
          resCharElemThree += (m_k[i])*(acoustics[i])+(m_kPast[i])*(acoustic_past[i]);
          resCharElemFour += (m_k[i])*(omega[i])+(m_kPast[i])*(omega_past[i]);
        }

// compute the distributed residual
        charResThree += m_betas*resCharElemThree;
        charResFour += m_betas*resCharElemFour;

/****************** End of matrix distribution in characteristic *******************/


      geoBuilder->releaseGE();

     }  // end of looping over the elements
     
     

  const CFreal avRho = linearData[EulerTerm::RHO];
  const CFreal avU   = linearData[EulerTerm::VX];
  const CFreal avV   = linearData[EulerTerm::VY];
  const CFreal avW   = linearData[EulerTerm::VZ];  
  const CFreal avH   = linearData[EulerTerm::H];
  const CFreal UU = (avU*avU + avV*avV + avW*avW);  
  const CFreal oOrho = 1.0/avRho;
     
     
      const CFreal drho = rhs(stateID, 0, m_nbEqs);
      const CFreal drhou = rhs(stateID, 1, m_nbEqs);
      const CFreal drhov = rhs(stateID, 2, m_nbEqs);
      const CFreal drhow = rhs(stateID, 3, m_nbEqs);
      const CFreal drhoE = rhs(stateID, 4, m_nbEqs);
      
      const CFreal LL = 2.0*(drhoE-drhou*avU+drhov*avV+drhow*avW)+drho*UU;
      const CFreal KK = -(gamma-1.0)*oneoverc*oneoverc*LL;

  
      charRes[0] = -(-ncharx*KK+oOrho*(ncharz*(drhov-drho*avV)-nchary*(drhow-drho*avW)+ncharx*drho*avRho));
      charRes[1] = -(-nchary*KK+oOrho*(ncharx*(drhow-drho*avW)-ncharz*(drhov-drho*avU)+nchary*drho*avRho));
      charRes[2] = -(-ncharz*KK+oOrho*(nchary*(drhou-drho*avU)-ncharx*(drhov-drho*avV)+ncharz*drho*avRho));
      charRes[3] = charResThree;
      charRes[4] = 0.0;
  
     // transform back to conservative
      Res[0] = ncharx*charRes[0]+nchary*charRes[1]+ncharz*charRes[2]+0.5*oneoverc*avRho*(charRes[3]);
      Res[1] = 0.5*oneoverc*avRho*avU*(charRes[3]) + avU*(ncharx*charRes[0]+nchary*charRes[1]+ncharz*charRes[2]) +avRho*(nchary*charRes[2]-nchary*charRes[1]);
      Res[2] = 0.5*oneoverc*avRho*avV*(charRes[3]) + avV*(ncharx*charRes[0]+nchary*charRes[1]+ncharz*charRes[2]) +avRho*(ncharz*charRes[0]-ncharx*charRes[2]);
      Res[3] = 0.5*oneoverc*avRho*avW*(charRes[3]) + avW*(ncharx*charRes[0]+nchary*charRes[1]+ncharz*charRes[2]) +avRho*(ncharx*charRes[1]-nchary*charRes[0]);
      Res[4] = 0.5*oneoverc*avRho*avH*(charRes[3]) +0.5*ncharx*(UU*charRes[0]-2.0*avRho*(avV*charRes[2]-avW*charRes[1]))+0.5*nchary*(UU*charRes[1]-2.0*(avW*charRes[0]-avU*charRes[2])) + 0.5*ncharz*(UU*charRes[2]-2.0*(avU*charRes[1]-avV*charRes[0]));
      

     for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq) {
       rhs(stateID, iEq, m_nbEqs) = -Res[iEq];
     }
    }

    isUpdated[stateID] = true; // flagging is important!!!!!
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletNonRefEuler3DCons::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "Euler3DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler3DCons>());

  cf_assert(_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletNonRefEuler3DCons::computeCharK(std::vector<Framework::State*>& states,   RealVector& _kPlus,   RealVector& _k, RealVector& _kPast, std::vector<RealVector*>& normal, RealVector& faceLength, CFreal& Area)
  {

  CFuint _nbStatesInCell = states.size();
  RealVector _nodeArea(states.size());
  RealVector _kspace = _k;

  const CFreal kCoeff = 1./3.;
  const CFreal tCoeff = SubSystemStatusStack::getActive()->getDT()/2.;

  const RealVector& linearData = _varSet->getModel()->getPhysicalData();
    
  const CFreal U0     = linearData[EulerTerm::VX];
  const CFreal V0     = linearData[EulerTerm::VY];
  const CFreal W0     = linearData[EulerTerm::VZ];  

    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState){
      _nodeArea[iState] = faceLength[iState];
      for (CFuint iDim = 0; iDim < 3; ++iDim) {
        m_adimNormal[iDim] = (*normal[iState])[iDim];
      }
      m_adimNormal *= 1./_nodeArea[iState];

    _kspace[iState] = U0*m_adimNormal[0]+V0*m_adimNormal[1]+W0*m_adimNormal[2];
    _kspace[iState] *= _nodeArea[iState] * kCoeff;

    _kPast[iState]  = _kspace[iState]*tCoeff - Area/4.;
    _k[iState] =_kspace[iState]*tCoeff + Area/4.;

    _kPlus[iState] = max(0.,_k[iState]);

    }

}


//////////////////////////////////////////////////////////////////////////////

CFreal StrongSubOutletNonRefEuler3DCons::distributeLDA(RealVector& m_kPlus, CFreal m_betas, CFuint boundarystate){

  m_sumKplus = m_kPlus[0];
  for (CFuint iState = 1; iState < 4; ++iState) {
    m_sumKplus  += m_kPlus[iState];
  }

  m_betas = (m_kPlus[boundarystate])/m_sumKplus;

  return m_betas;
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletNonRefEuler3DCons::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
