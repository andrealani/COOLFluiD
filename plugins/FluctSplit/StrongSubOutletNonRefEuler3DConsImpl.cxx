#include <numeric>
#include <algorithm>

#include "StrongSubOutletNonRefEuler3DConsImpl.hh"
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
#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"

#include "Framework/BlockAccumulator.hh"
#include "Framework/LSSMatrix.hh"
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

MethodCommandProvider<StrongSubOutletNonRefEuler3DConsImpl, FluctuationSplitData, FluctSplitSpaceTimeNavierStokesModule> strongSubFunctionOutletNonRefLinEuler3DConsImplProvider("StrongSubOutletNonRefEuler3DConsImpl");

//////////////////////////////////////////////////////////////////////////////

StrongSubOutletNonRefEuler3DConsImpl::StrongSubOutletNonRefEuler3DConsImpl(const std::string& name) :
  FluctuationSplitCom(name),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_updateCoeff("updateCoeff"),
  socket_nodes("nodes"),
  socket_isBState("isBState"),
  socket_volumes("volumes"),
  socket_pastStates("pastStates"),
  m_Statenormals(0), 
  _varSet(),
  _bcNormals(),
  _jacobElem(5,5),
  _jacob(5,5),
  _jacobAll(),
  bndNod2Elm(0),  
  _in(5),
  _ira(5),  
  m_nbEqs(),
  nbFreedom()
{

}

//////////////////////////////////////////////////////////////////////////////

StrongSubOutletNonRefEuler3DConsImpl::~StrongSubOutletNonRefEuler3DConsImpl()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSubOutletNonRefEuler3DConsImpl::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_nodes);
  result.push_back(&socket_isBState);  
  result.push_back(&socket_volumes);
    result.push_back(&socket_pastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletNonRefEuler3DConsImpl::setup()
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
  
/// store the number of freedoms
  
  nbFreedom.resize(bndNbNodes);
  
  RealVector ijIDs;
  
  const LSSIdxMapping& idxMapping =
    getMethodData().getLinearSystemSolver()[0]->getLocalToGlobalMapping();
    
  // prepares to loop over cells by getting the GeometricEntityPool
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = MeshDataStack::getActive()->getTrs("InnerCells");
  
  for (CFuint iState = 0; iState < bndNbNodes; ++iState) {
    
    vector <CFuint> cells = bndNod2Elm[iState];
    CFuint nbElemsOfState = cells.size();
    
    ijIDs.resize(nbElemsOfState*4);
    RealVector selectedIDs;
    selectedIDs.resize(nbElemsOfState*4);
    
    for(CFuint h=0; h<nbElemsOfState*4; h++) {
      ijIDs[h] = 0;
      selectedIDs [h] = 0;
    }
    
      for (CFuint elem=0; elem<nbElemsOfState; ++elem) {
 
        // build the GeometricEntity
        CFuint cellID = cells[elem];
        geoData.idx = cellID;
        GeometricEntity& cell = *geoBuilder->buildGE();
        vector<State*> *const statesInCell = cell.getStates();
        const CFuint nbStatesInCell = statesInCell->size();

        for (CFuint elemState=0; elemState<nbStatesInCell; ++elemState) {
	  State *const currState = (*statesInCell)[elemState];
	  CFuint currID=currState->getLocalID();
	  ijIDs[elem*nbStatesInCell+elemState] = currID;
	}
     }
   
   CFuint freedom=1;
  
   selectedIDs[0] = ijIDs[0];
   for (CFuint numfree=1; numfree<ijIDs.size(); numfree++) {
     bool found = false;
     for (CFuint num=0; num<freedom; num++) {
	if (ijIDs[numfree]==selectedIDs[num]) {
	   found = true;
	}
     }
     if (!found) {
	selectedIDs[freedom] = ijIDs[numfree];
	freedom += 1;
	}
     }
   
  nbFreedom[iState] = freedom;
  }
   

   
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletNonRefEuler3DConsImpl::executeOnTrs()
{
  Common::SafePtr<LinearSystemSolver> lss =
  getMethodData().getLinearSystemSolver()[0];

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();
  jacobMatrix->flushAssembly(); 
  
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

  // block accumulator 1*1
  auto_ptr<BlockAccumulator> acc(lss->createBlockAccumulator(1, 1, m_nbEqs));
  
  const LSSIdxMapping& idxMapping =
    getMethodData().getLinearSystemSolver()[0]->getLocalToGlobalMapping();

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
      
      
       CFuint NNodes = 4*nbElemsOfState;
	    
       vector <RealMatrix> jacobians;
       jacobians.resize(NNodes);
	
	for(CFuint ii=0; ii<NNodes; ii++) {
	  jacobians[ii].resize(m_nbEqs,m_nbEqs);
	  for(CFuint jj=0; jj<m_nbEqs; jj++)
	    for(CFuint kk=0; kk<m_nbEqs; kk++)
	    (jacobians[ii])(jj,kk) = 0.0;
	}

	CFreal Dii[4];
        CFreal D23[4];

	for(CFuint ii=0; ii<4; ii++) {	
	  Dii[ii]=0.;
	  D23[ii]=0.;
	}
	
      CFuint ijIDs[NNodes];    	
  
    CFuint elemCount = 0;   
       
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
	
	for (CFuint elemState=0; elemState<nbStatesInCell; ++elemState) {
	  State *const currState = (*statesInCell)[elemState];
	  CFuint currID=currState->getLocalID();
	  CFuint currGlobalID = idxMapping.getColID(currID)*m_nbEqs;
	  ijIDs[elemCount*nbStatesInCell+elemState] = currID;
	}

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

	for (CFuint iStat = 0; iStat < nbStatesInCell; ++iStat) {
	  Dii[iStat] = m_betas*m_k[iStat];
// 	  D23[iStat] = m_betas*0.25*c/Area*(nchary*(*m_Statenormals[iStat])[0]-ncharx*(*m_Statenormals[iStat])[1]);
	  //this is just an approximation, but works fine
	  D23[iStat] = 0.0;
	}
	// First approximation
        for (CFuint iJac = 0; iJac < nbStatesInCell; ++iJac) {
	  if (iJac == boundaryState) {
	   jacobians[elemCount*nbStatesInCell+iJac](0,0) = Dii[iJac];
	   jacobians[elemCount*nbStatesInCell+iJac](1,1) = Dii[iJac];
	   jacobians[elemCount*nbStatesInCell+iJac](2,2) = Dii[iJac];   
	   jacobians[elemCount*nbStatesInCell+iJac](3,3) =  Dii[iJac];
	   jacobians[elemCount*nbStatesInCell+iJac](4,4) =  Dii[iJac];
	  }
	  else {}
	}
	
/****************** End of matrix distribution in characteristic *******************/

      elemCount += 1;
      
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
     
     
// comupute the JACOBIANS ***************************************************        


      CFuint nbEqs = m_nbEqs;

//merge all the contributions
       CFint reducedset = nbFreedom[iState];
       vector <RealMatrix> jacobian_merged;
       jacobian_merged.resize(reducedset);
       CFuint jacobianIDs[reducedset];
	
	for(CFuint ii=0; ii<reducedset; ii++) {
	  (jacobian_merged[ii]).resize(nbEqs,nbEqs);
	  for(CFuint jj=0; jj<nbEqs; jj++)
	    for(CFuint kk=0; kk<nbEqs; kk++)
	    (jacobian_merged[ii])(jj,kk) = 0.0;
	}


     CFuint setted = 1;
// set the first item     
     jacobian_merged[0]=jacobians[0];
     jacobianIDs[0]=ijIDs[0];
// merge contributions
     for (CFuint Inclnode = 1; Inclnode < NNodes; ++Inclnode) {
        const CFuint entryID = ijIDs[Inclnode];
        bool found = false;
	for(CFuint ii=0; ii<setted; ii++) {
	  if(entryID==jacobianIDs[ii]) {
	    found = true;
            jacobian_merged[ii]+=jacobians[Inclnode];
	  }
	}
        if(found==false) {
            jacobian_merged[setted]=jacobians[Inclnode];
	    jacobianIDs[setted]=ijIDs[Inclnode];
	    setted+=1;
	}
     }
   
     acc->setRowIndex(0, stateID);
   
     for (CFuint Inclnode = 0; Inclnode < reducedset; ++Inclnode) {
        const CFuint entryID = jacobianIDs[Inclnode];

	acc->setColIndex(0, entryID);

        _jacob =  jacobian_merged[Inclnode];

	for (CFuint iVar = 0; iVar < nbEqs; ++iVar) {
	  for (CFuint jVar = 0; jVar < nbEqs; ++jVar) {
	    acc->setValue(0,0, iVar, jVar,
			  _jacob(iVar,jVar));
          }
         }
        jacobMatrix->setValues(*acc);
       }
       
       jacobians.resize(0);      
       
      }

    isUpdated[stateID] = true; // flagging is important!!!!!
    }
    
  }
  
  // this should be an intermediate lightweight assembly !!!
  // it is needed because here you SET values while elsewhere
  // you ADD values
  jacobMatrix->flushAssembly();
   
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletNonRefEuler3DConsImpl::configure ( Config::ConfigArgs& args )
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

void StrongSubOutletNonRefEuler3DConsImpl::computeCharK(std::vector<Framework::State*>& states,   RealVector& _kPlus,   RealVector& _k, RealVector& _kPast, std::vector<RealVector*>& normal, RealVector& faceLength, CFreal& Area)
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

CFreal StrongSubOutletNonRefEuler3DConsImpl::distributeLDA(RealVector& m_kPlus, CFreal m_betas, CFuint boundarystate){

  m_sumKplus = m_kPlus[0];
  for (CFuint iState = 1; iState < 4; ++iState) {
    m_sumKplus  += m_kPlus[iState];
  }

  m_betas = (m_kPlus[boundarystate])/m_sumKplus;

  return m_betas;
}


//////////////////////////////////////////////////////////////////////////////
void StrongSubOutletNonRefEuler3DConsImpl::unsetup()
{
  _bcNormals.resize(0);
  _jacob.resize(0,0);
  _jacobElem.resize(0,0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
