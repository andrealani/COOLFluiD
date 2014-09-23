#include <numeric>
#include <algorithm>

#include "StrongSubOutletNonRefEuler2DCons.hh"
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

MethodCommandProvider<StrongSubOutletNonRefEuler2DCons, FluctuationSplitData, FluctSplitSpaceTimeNavierStokesModule> strongSubOutletNonRefEuler2DConsProvider("StrongSubOutletNonRefEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

StrongSubOutletNonRefEuler2DCons::StrongSubOutletNonRefEuler2DCons(const std::string& name) :
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

StrongSubOutletNonRefEuler2DCons::~StrongSubOutletNonRefEuler2DCons()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSubOutletNonRefEuler2DCons::needsSockets()
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

void StrongSubOutletNonRefEuler2DCons::setup()
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

  m_adimCharNormal.resize(2);
 
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletNonRefEuler2DCons::executeOnTrs()
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

      CFreal charResTwo = 0.0;
      CFreal charResThree = 0.0;
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

      RealVector bcNormalState =  bndNodNorm[iState];

      CFreal ncharx = (bcNormalState)[0];
      CFreal nchary = (bcNormalState)[1];

 
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

	CFuint boundaryState = 100;
        for (CFuint elemState=0; elemState<nbStatesInCell; ++elemState) {
          State *const currState = (*statesInCell)[elemState];
          if(stateID == currState->getLocalID()) {
            boundaryState = elemState;
          }
        }
	
        for (CFuint elemState=0; elemState<nbStatesInCell; ++elemState) {
          State *const currState = (*statesInCell)[elemState];
	  CFreal p;
          p = (gamma-1.0)*((*currState)[3]-
				  0.5*((*currState)[1]*(*currState)[1]+(*currState)[2]*(*currState)[2])
				  /(*currState)[0]);
          acoustics[elemState]= 2.0*oneoverc*(p)/(*currState)[0];
          omega[elemState]= 0.0;

        }
        
        
        

/********************* Matrix distribution in characteristic ***********************/

        CFreal resCharElemTwo=0.;
        CFreal resCharElemThree=0.;

        RealVector faceLength(nbStatesInCell);

        for (CFuint i = 0; i < nbStatesInCell; ++i) {
          (*m_Statenormals[i])[0] = m_normals[cellID]->getNodalNormComp(i,0);
          (*m_Statenormals[i])[1] = m_normals[cellID]->getNodalNormComp(i,1);
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

	  CFreal p = (gamma-1.0)*((currState_past)[3]-
				  0.5*((currState_past)[1]*(currState_past)[1]+(currState_past)[2]*(currState_past)[2])
				  /(currState_past)[0]);
          acoustic_past[i]= 2.0 * oneoverc*(p)/(currState_past)[0];
          omega_past[i]= 0.0;
        }
        
        for (CFuint i = 0; i < nbStatesInCell; ++i) {
          resCharElemTwo += (m_k[i])*(acoustics[i])+(m_kPast[i])*(acoustic_past[i]);
          resCharElemThree += (m_k[i])*(omega[i])+(m_kPast[i])*(omega_past[i]);
        }

// compute the distributed residual
        charResTwo += m_betas*resCharElemTwo;
        charResThree += m_betas*resCharElemThree;

/****************** End of matrix distribution in characteristic *******************/


      geoBuilder->releaseGE();

     }  // end of looping over the elements


// this is entropy, good as it is.
      charRes[0] = -(rhs(stateID, 0, m_nbEqs) - rhs(stateID, 3, m_nbEqs)*oneoverc*oneoverc);

// this is vorticity in the direction where acoustic wave is propagating, so they are decpoupled

      charRes[1] = -(rhs(stateID, 1, m_nbEqs)*nchary - rhs(stateID, 2, m_nbEqs)*ncharx);
      charRes[2] = charResTwo;
//       charRes[3] = charResThree;
      charRes[3] = 0.0;

// transform back to conservative
  const CFreal avRho = linearData[EulerTerm::RHO];
  const CFreal avU   = linearData[EulerTerm::VX];
  const CFreal avV   = linearData[EulerTerm::VY];
  const CFreal avW   = linearData[EulerTerm::VZ];  
  const CFreal avE   = linearData[EulerTerm::E];
  const CFreal UU = (avU*avU + avV*avV);
    
  
      Res[0] = (charRes[0]+0.5*oneoverc*avRho*(charRes[2]));
      Res[1] = 0.5*oneoverc*avRho*avU*(charRes[2])+nchary*avRho*charRes[1] + avU*charRes[0];
      Res[2] = 0.5*oneoverc*avRho*avV*(charRes[2])-ncharx*avRho*charRes[1] + avV*charRes[0];
      Res[3] = 0.25*oneoverc* ( -(gamma-1.)*avRho*(charRes[2])*UU +4.0*c*avRho*(nchary*avU - ncharx*avV)*charRes[1] +2.0*c*charRes[0]*UU );


// get average state variables
//   const CFreal avRho = linearData[EulerTerm::RHO];
//   const CFreal avU   = linearData[EulerTerm::VX];
//   const CFreal avV   = linearData[EulerTerm::VY];
//   const CFreal avW   = linearData[EulerTerm::VZ];  
//   const CFreal avE   = linearData[EulerTerm::E];
//   const CFreal UU = (avU*avU + avV*avV);
// 
//       charRes[0] = -((1-0.5*(avU*avU+avV*avV)*(gamma-1)*oneoverc*oneoverc)*rhs(stateID, 0, m_nbEqs)
//                    +avU*(gamma-1)*oneoverc*oneoverc*rhs(stateID, 1, m_nbEqs) +avV*(gamma-1)*oneoverc*oneoverc*rhs(stateID, 2, m_nbEqs)
// 		   -(gamma-1)*oneoverc*oneoverc*rhs(stateID, 3, m_nbEqs));
//       
//      
//       charRes[1] = -((-avU*nchary+avV*ncharx)/avRho*rhs(stateID, 0 ,m_nbEqs) + nchary/avRho*rhs(stateID, 1, m_nbEqs)
//                    -ncharx/avRho*rhs(stateID, 2, m_nbEqs));
// 
//       charRes[2] = charResTwo;
// 
//       charRes[3] = 0.0;  
//       
//       //Transform back to conservative variables
//       Res[0] = (charRes[0] + 0.5*oneoverc*avRho*(charRes[2]));
//       
//       Res[1] = avU*charRes[0] + avRho*nchary*charRes[1] + 0.5*oneoverc*avU*avRho*charRes[2] + 0.5*avRho*ncharx*charRes[3];
//       
//       Res[2] = avV*charRes[0] - avRho*ncharx*charRes[1] + 0.5*oneoverc*avRho*avV*charRes[2] + 0.5*avRho*nchary*charRes[3];
// 
//       
//       Res[3] = 0.5*(avU*avU + avV*avV)*charRes[0] + avRho*(avU*nchary-avV*ncharx)*charRes[1] 
//                + 0.5*avRho*(0.5*oneoverc*(avU*avU + avV*avV) + c/(gamma-1))*charRes[2] + 0.5*avRho*(avU*ncharx+avV*nchary)*charRes[3];

     for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq) {
       rhs(stateID, iEq, m_nbEqs) = -Res[iEq];
     }
    }

    isUpdated[stateID] = true; // flagging is important!!!!!
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletNonRefEuler2DCons::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "Euler2DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler2DCons>());

  cf_assert(_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletNonRefEuler2DCons::computeCharK(std::vector<Framework::State*>& states,   RealVector& _kPlus,   RealVector& _k, RealVector& _kPast, std::vector<RealVector*>& normal, RealVector& faceLength, CFreal& Area)
  {

  CFuint _nbStatesInCell = states.size();
  RealVector _nodeArea(states.size());
  RealVector _kspace = _k;

  const CFreal kCoeff = 1./2.;
  const CFreal tCoeff = SubSystemStatusStack::getActive()->getDT()/2.;

  const RealVector& linearData = _varSet->getModel()->getPhysicalData();
    
  const CFreal U0     = linearData[EulerTerm::VX];
  const CFreal V0     = linearData[EulerTerm::VY];

    for (CFuint iState = 0; iState < _nbStatesInCell; ++iState){
      _nodeArea[iState] = faceLength[iState];
      for (CFuint iDim = 0; iDim < 2; ++iDim) {
        m_adimNormal[iDim] = (*normal[iState])[iDim];
      }
      m_adimNormal *= 1./_nodeArea[iState];

    _kspace[iState] = U0*m_adimNormal[0]+V0*m_adimNormal[1];
    _kspace[iState] *= _nodeArea[iState] * kCoeff;

    _kPast[iState]  = _kspace[iState]*tCoeff - Area/3.;
    _k[iState] =_kspace[iState]*tCoeff + Area/3.;

//     _k[iState]  = _k[iState]*tCoeff + Area;
    _kPlus[iState] = max(0.,_k[iState]);

    }

}


//////////////////////////////////////////////////////////////////////////////

CFreal StrongSubOutletNonRefEuler2DCons::distributeLDA(RealVector& m_kPlus, CFreal m_betas, CFuint boundarystate){

  m_sumKplus = m_kPlus[0];
  for (CFuint iState = 1; iState < 3; ++iState) {
    m_sumKplus  += m_kPlus[iState];
  }

  m_betas = (m_kPlus[boundarystate])/m_sumKplus;

  return m_betas;
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletNonRefEuler2DCons::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
