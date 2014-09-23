
#include "Common/BadValueException.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"

#include "FluctSplit/STU_HOCRD_SplitStrategy.hh"
#include "FluctSplit/SpaceTime_Splitter.hh"
#include "FluctSplit/FluctSplitSpaceTime.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<STU_HOCRD_SplitStrategy,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitSpaceTimeModule>
STU_HOCRDFluctSplitStrategyProvider("STU_HOCRD");

//////////////////////////////////////////////////////////////////////////////

void STU_HOCRD_SplitStrategy::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("isFirstIteP1","Choose that the first time iteration is done as a P1 in space");


}
//////////////////////////////////////////////////////////////////////////////
STU_HOCRD_SplitStrategy::STU_HOCRD_SplitStrategy(const std::string& name) :
  FluctuationSplitStrategy(name),
  CellVolume(),
  socket_updateCoeff("updateCoeff"),
  socket_interStates("interStates"),
  m_interStates(),
  m_splitter(CFNULL),
  m_phi_inter(),
  m_phi_present(),
  qdstates(),
  qdExtraVars(),
  m_statesBkp(),
  facenormal(),
  subfacetable(0,0),
  qd0a(),
  qd1a(),
  wqda(),
  m_updateVar(CFNULL),
  faceflux(),
  qdnodes(),
  subelemfacedir(0,0),
  subelemtable(0,0),
  subelemtable_state(0,0),
  substates(),
  temp_states(),
  m_phi_time(),
  qd0t(),
  qd1t(),
  qd2t(),
  wqdt(),
  m_u(),
  m_flux_time(),
  m_flux_space(),
  subresidual(),
  socket_pastStates("pastStates"),
  m_pastStates(),
  m_phi_past(),
  m_interResiduals(0),
  m_pastResiduals(0)
{
  addConfigOptionsTo(this);

  m_firstitP1 = false;
  setParameter("isFirstIteP1",&m_firstitP1);
}

//////////////////////////////////////////////////////////////////////////////

STU_HOCRD_SplitStrategy::~STU_HOCRD_SplitStrategy()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void STU_HOCRD_SplitStrategy::setup()
{
  CFAUTOTRACE;

  // first call parent method
  FluctuationSplitStrategy::setup();

  const CFuint maxNbStatesInCell = MeshDataStack::getActive()->Statistics().getMaxNbStatesInCell();
  DistributionData& ddata = getMethodData().getDistributionData();
  m_interStates.resize(maxNbStatesInCell);
  temp_states.resize(maxNbStatesInCell);
  m_pastStates.resize(maxNbStatesInCell);
  // Resizing all the states of the past and temporary states 
  for (CFuint i = 0; i < maxNbStatesInCell; ++i) {
    m_interStates[i] = new State();
    temp_states[i] = new State();
    m_pastStates[i] = new State();
   }

  // Resize the backup states
  m_statesBkp.resize(maxNbStatesInCell);

  // get the splitter
  m_splitter = getMethodData().getSplitter().d_castTo<SpaceTime_Splitter>();

  // get the update coeff socket
  DataHandle< CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  m_splitter->setUpdateCoeff(updateCoeff);

  // resize the fluctruation for the present  and intermediate layer 
  // of the sub-elements (advective part)
  m_phi_inter.resize(4); // 4 sub-elements;
  m_phi_present.resize(4); // 4 sub-elements;
  m_phi_past.resize(4); // 4 sub-elements;

  m_nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for (CFuint i = 0; i< 4; ++ i)
    {
      m_phi_inter[i].resize(m_nbEqs);
      m_phi_present[i].resize(m_nbEqs);
      m_phi_past[i].resize(m_nbEqs);

    }
  
 // physical data evaluated in the quadrature points
  m_pdata.resize(4);
  for (CFuint  i = 0; i < 4; ++i) {
    PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->
      resizePhysicalData(m_pdata[i]);
  }  
  
  qdstates.resize(3); // 3 quadrature points per face
  qdstates[0] = new State();
  qdstates[1] = new State();
  qdstates[2] = new State();

  qdExtraVars.resize(3); // 3 quadrature points per face
  qdExtraVars[0] = new RealVector();
  qdExtraVars[1] = new RealVector();
  qdExtraVars[2] = new RealVector();
  
  qdnodes.resize(3); // 3 quadrature points per face
  qdnodes[0] = new Node();
  qdstates[0]->setSpaceCoordinates(qdnodes[0]);
  qdnodes[1] = new Node();
  qdstates[1]->setSpaceCoordinates(qdnodes[1]);
  qdnodes[2] = new Node();
  qdstates[2]->setSpaceCoordinates(qdnodes[2]);

  facenormal.resize(2); //working only in 2D

  // sub face table : contain the node that contain each face
  subfacetable.resize(9,2); // 9 sub faces with 2 states each

  subfacetable(0,0) = 3;
  subfacetable(0,1) = 5;

  subfacetable(1,0) = 5;
  subfacetable(1,1) = 0;

  subfacetable(2,0) = 0;
  subfacetable(2,1) = 3;

  subfacetable(3,0) = 1;
  subfacetable(3,1) = 4;

  subfacetable(4,0) = 4;
  subfacetable(4,1) = 3;

  subfacetable(5,0) = 3;
  subfacetable(5,1) = 1;

  subfacetable(6,0) = 4;
  subfacetable(6,1) = 2;

  subfacetable(7,0) = 2;
  subfacetable(7,1) = 5;

  subfacetable(8,0) = 5;
  subfacetable(8,1) = 4;

  // sub elemt table : contain the faces of each sub-element
  subelemtable.resize(4,3); // 4 sub elems with 3 faces each
  subelemtable_state.resize(4,3); // 4 sub elems with 3 states each

  subelemfacedir.resize(4,3);

  subelemtable(0,0) = 0;
  subelemtable(0,1) = 1;
  subelemtable(0,2) = 2;

  subelemtable(1,0) = 3;
  subelemtable(1,1) = 4;
  subelemtable(1,2) = 5;

  subelemtable(2,0) = 6;
  subelemtable(2,1) = 7;
  subelemtable(2,2) = 8;

  subelemtable(3,0) = 0; // faces have negative orientation
  subelemtable(3,1) = 4;
  subelemtable(3,2) = 8;

// sub elemt table : contain the orientation of the faces of each sub-element
  subelemfacedir(0,0) = 1.;
  subelemfacedir(0,1) = 1.;
  subelemfacedir(0,2) = 1.;

  subelemfacedir(1,0) = 1.;
  subelemfacedir(1,1) = 1.;
  subelemfacedir(1,2) = 1.;

  subelemfacedir(2,0) = 1.;
  subelemfacedir(2,1) = 1.;
  subelemfacedir(2,2) = 1.;

  subelemfacedir(3,0) = -1.; // faces have negative orientation
  subelemfacedir(3,1) = -1.;
  subelemfacedir(3,2) = -1.;

  subelemtable_state(0,0) = 0;
  subelemtable_state(0,1) = 3;
  subelemtable_state(0,2) = 5;

  subelemtable_state(1,0) = 3;
  subelemtable_state(1,1) = 1;
  subelemtable_state(1,2) = 4;

  subelemtable_state(2,0) = 5;
  subelemtable_state(2,1) = 4;
  subelemtable_state(2,2) = 2;

  subelemtable_state(3,0) = 4; 
  subelemtable_state(3,1) = 5;
  subelemtable_state(3,2) = 3;


  const CFreal s  = std::sqrt( 0.6 );
  const CFreal a0 = ( 1.0 - s )*0.5;
  const CFreal a1 = ( 1.0 + s )*0.5;
  
  

  qd0a.resize(3); // quadrature points per face
  qd1a.resize(3); // quadrature points per face

  qd0a[0] = a0;  qd1a[0] = a1;
  qd0a[1] = a1;  qd1a[1] = a0;
  qd0a[2] = .5;  qd1a[2] = .5;

  wqda.resize(3); // 3 quadrature points per face
  wqda[0] = 5.0/18.0;
  wqda[1] = 5.0/18.0;
  wqda[2] = 8.0/18.0;

  m_updateVar = getMethodData().getUpdateVar();

  faceflux.resize(subfacetable.nbRows()); // one flux per sub face
  for (CFuint i = 0; i < faceflux.size(); ++i)
    faceflux[i].resize(m_nbEqs);

  (ddata.phi_time).resize(m_nbEqs);
  // sub elemt table
  substates.resize(3);   // 3 states in each sub element
  subresidual.resize(3); // 3 residuals in each sub element
  for(CFuint i = 0; i < 3; ++i)
    subresidual[i].resize(m_nbEqs);

  m_phi_time.resize(m_nbEqs);
  m_flux_time.resize(m_nbEqs);
  m_flux_space.resize(m_nbEqs);

  qd0t.resize(4); // 4 quadrature points per triangle
  qd1t.resize(4); // 4 quadrature points per triangle
  qd2t.resize(4); // 4 quadrature points per triangle

  qd0t[0] = 1.0/3.0;  qd1t[0] = 1.0/3.0;  qd2t[0] = 1.0/3.0;
  qd0t[1] = 0.6;  qd1t[1] = 0.2;  qd2t[1] = 0.2;  
  qd0t[2] = 0.2;  qd1t[2] = 0.6;  qd2t[2] = 0.2;   
  qd0t[3] = 0.2;  qd1t[3] = 0.2;  qd2t[3] = 0.6; 

  wqdt.resize(4); // 4 quadrature points per triangle
  // the weight is devided by four because we use it to integreate on a sub_element
  // which area is volume/4
  wqdt[0] = -9.0/64.0;
  wqdt[1] = 25.0/192.0;
  wqdt[2] = 25.0/192.0;
  wqdt[3] = 25.0/192.0;
  m_u.resize(m_nbEqs);

  ddata.isHO = true;


  Common::SafePtr<TopologicalRegionSet> innerCells = MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = innerCells->getLocalNbGeoEnts();
  CFuint nbStatesInCell;
  CFuint nbsubElm = 4;
// Resizing the storage for the full past Residuals
  m_interResiduals.resize(nbCells);
  m_pastResiduals.resize(nbCells);
  // Resize the RealVector
  _stdTrsGeoBuilder.getDataGE().trs = innerCells;

  for(CFuint iGeoEnt = 0; iGeoEnt < nbCells; ++iGeoEnt) {

    CFLogDebugMax("Cell " << iGeoEnt << "\n");

    // build the GeometricEntity
    _stdTrsGeoBuilder.getDataGE().idx = iGeoEnt;
    GeometricEntity& cell = *_stdTrsGeoBuilder.buildGE();

    nbStatesInCell = cell.nbStates();
    m_interResiduals[iGeoEnt].resize(m_nbEqs*nbsubElm);
    m_pastResiduals[iGeoEnt].resize(m_nbEqs*nbsubElm);
    //release the GeometricEntity
    _stdTrsGeoBuilder.releaseGE();
  }

  
}

//////////////////////////////////////////////////////////////////////////////

void STU_HOCRD_SplitStrategy::computeFluctuation(vector<RealVector>& residual)
{
  // The first iteration needs a special function because at this point
  // only two layer are available althought we need three.

 CFuint nbIter = SubSystemStatusStack::getActive()->getNbIter();
 if(nbIter == 1) {
   if(!m_firstitP1){
    
       docomputeFirstFluctuation(residual);
 
   }
   else
     {

        docomputeFirstFluctuationP1(residual);

     }
   }
 else{ docomputeFluctuation(residual);

 }
}
// //////////////////////////////////////////////////////////////////////////////

void STU_HOCRD_SplitStrategy::docomputeFirstFluctuation(vector<RealVector>& residual)
{
  ///@todo Nv: There is no moving mesh implemented
  DistributionData& ddata = getMethodData().getDistributionData();
  vector<State*>& states = *ddata.states;
  const CFuint nbStatesInCell = states.size();
//  const CFuint nbStatesInSubCell = 3;
  ddata.isfirstP1 = false;

  m_splitter->setDT(SubSystemStatusStack::getActive()->getDT());
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();
  const CFreal timeStep = dt/dim;
 
  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [ddata.cellID]);
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  CellVolume = volumes[ddata.cellID];
  // backup the states of this cell
  setCurrentCell();

  
  // Set the cell volume such that it can be accessed by the splitter
  m_splitter->setCellVolume(volumes[ddata.cellID]);

  // Since there is no moving mesh here, we set the past cell volume to the same as 
  // the non moving one
  m_splitter->setPastCellVolume(volumes[ddata.cellID]);

  // reset the residual because we will accumulate the sub element contributions
  for (CFuint i = 0; i < nbStatesInCell; ++i)
  {
    residual[i] = 0.0;
  }  
//get the interstates in this cell (intermediate layer)
  DataHandle<State*> interStatesStorage = socket_interStates.getDataHandle();
  for (CFuint i = 0; i < nbStatesInCell; ++i) {
    const CFuint stateID = (*ddata.states)[i]->getLocalID();
   State & temp =  *m_interStates[i];
  //We need the local ID when we do Linearized Euler
    //to access the meanflow 
    temp.clone(*(*ddata.states)[i]);
    *m_interStates[i] = *interStatesStorage[stateID];
  }
  if (SubSystemStatusStack::getActive()->isFirstStep()){
  // We compute the advective fluctuation on the present layer and the
  // intermediate one and we store it
  ddata.time = SubSystemStatusStack::getActive()->getCurrentTime();
  computeHOSTFluxIntegral(m_interStates,m_phi_inter);
  for (CFuint iSubElm = 0; iSubElm < 4; ++iSubElm){
      for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
        m_interResiduals[ddata.cellID][(iSubElm*m_nbEqs)+iEq] = m_phi_inter[iSubElm][iEq];
    }
  }

}
if (!SubSystemStatusStack::getActive()->isFirstStep()){
//If it is not the first iteration then m_phi_inter has already been computed 
  for (CFuint iSubElm = 0; iSubElm < 4; ++iSubElm){
      for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
         m_phi_inter[iSubElm][iEq] = m_interResiduals[ddata.cellID][(iSubElm*m_nbEqs)+iEq];
    }
  }
}
  ddata.time = SubSystemStatusStack::getActive()->getCurrentTime() + dt;
  computeHOSTFluxIntegral(states,m_phi_present);
  ddata.time = SubSystemStatusStack::getActive()->getCurrentTime();
   /*****         Triangle 0-3-5          *****/

  substates[0] = states[0];
  substates[1] = states[3];
  substates[2] = states[5];

  computeHOtime(0,3,5);

  // linearize the states in the cell
  // includes the transformation from update to linearization
  // variables to evaluate the jacobians in the average state
  // and then do the transformation from update to consistent variables
  
  ddata.tStates = computeConsistentStates(&substates);
  ddata.subStates = &substates;

  for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
    m_flux_space[iEq] = timeStep*(m_phi_inter[0][iEq] + m_phi_present[0][iEq])+m_phi_time[iEq];
  }

  //Need to scale the cell normals to the half
  cellnormals.scale(0.5);
  m_splitter->computeK(substates, &cellnormals);
  //Need to unscale after
  cellnormals.unscale();

  //Transform back the residuals to the distributed variables
  SafePtr<RealVector> phi = &ddata.phi;
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_flux_space);


  // Distribute the residual to the states of the present layer
  m_splitter->distribute(subresidual);

  residual[0] += subresidual[0];
  residual[3] += subresidual[1];
  residual[5] += subresidual[2];

   /*****         Triangle 3-1-4          *****/
  substates[0] = states[3];
  substates[1] = states[1];
  substates[2] = states[4];

  computeHOtime(3,1,4);
  // linearize the states in the cell
  // includes the transformation from update to linearization
  // variables to evaluate the jacobians in the average state
  // and then do the transformation from update to consistent variables
  
  ddata.tStates = computeConsistentStates(&substates);
  ddata.subStates = &substates;

  for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
    m_flux_space[iEq] = timeStep*(m_phi_inter[1][iEq] + m_phi_present[1][iEq])+m_phi_time[iEq];
  }

  //Need to scale the cell normals to the half
  cellnormals.scale(0.5);
  m_splitter->computeK(substates, &cellnormals);
  //Need to unscale after
  cellnormals.unscale();

  //Transform back the residuals to the distributed variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_flux_space);

  // Distribute the residual to the states of the present layer
  m_splitter->distribute(subresidual);

  residual[3] += subresidual[0];
  residual[1] += subresidual[1];
  residual[4] += subresidual[2];
  
 /*****         Triangle 5-4-2          *****/
  substates[0] = states[5];
  substates[1] = states[4];
  substates[2] = states[2];

  computeHOtime(5,4,2);

  // linearize the states in the cell
  // includes the transformation from update to linearization
  // variables to evaluate the jacobians in the average state
  // and then do the transformation from update to consistent variables
  
  ddata.tStates = computeConsistentStates(&substates);
  ddata.subStates = &substates;

  for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
    m_flux_space[iEq] = timeStep*(m_phi_inter[2][iEq] + m_phi_present[2][iEq])+m_phi_time[iEq];
  }

  //Need to scale the cell normals to the half
  cellnormals.scale(0.5);
  m_splitter->computeK(substates, &cellnormals);
  //Need to unscale after
  cellnormals.unscale();

  //Transform back the residuals to the distributed variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_flux_space);

  // Distribute the residual to the states of the present layer
  m_splitter->distribute(subresidual);

  residual[5] += subresidual[0];
  residual[4] += subresidual[1];
  residual[2] += subresidual[2];

 /*****         Triangle 4-5-3          *****/
  substates[0] = states[4];
  substates[1] = states[5];
  substates[2] = states[3];

  computeHOtime(4,5,3);

  // linearize the states in the cell
  // includes the transformation from update to linearization
  // variables to evaluate the jacobians in the average state
  // and then do the transformation from update to consistent variables
  
  ddata.tStates = computeConsistentStates(&substates);
  ddata.subStates = &substates;

  for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
    m_flux_space[iEq] = timeStep*(m_phi_inter[3][iEq] + m_phi_present[3][iEq])+m_phi_time[iEq];
  }
  //Need to scale the cell normals to the half
  cellnormals.scale(-0.5);
  m_splitter->computeK(substates, &cellnormals);
  //Need to unscale after
  cellnormals.unscale();
  //Transform back the residuals to the distributed variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_flux_space);
  // Distribute the residual to the states of the present layer
  m_splitter->distribute(subresidual);
  residual[4] += subresidual[0];
  residual[5] += subresidual[1];
  residual[3] += subresidual[2];
}

// //////////////////////////////////////////////////////////////////////////////

void STU_HOCRD_SplitStrategy::docomputeFirstFluctuationP1(vector<RealVector>& residual)
{
  ///@todo Nv: There is no moving mesh implemented
  DistributionData& ddata = getMethodData().getDistributionData();
  vector<State*>& states = *ddata.states;
  const CFuint nbStatesInCell = states.size();
//  const CFuint nbStatesInSubCell = 3;
ddata.isfirstP1 = true;
  m_splitter->setDT(SubSystemStatusStack::getActive()->getDT());
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();
  const CFreal timeStep = dt/dim;
  
  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [ddata.cellID]);
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  CellVolume = volumes[ddata.cellID];
  // backup the states of this cell
  setCurrentCell();

  
  // Set the cell volume such that it can be accessed by the splitter
  m_splitter->setCellVolume(volumes[ddata.cellID]);
  // Since there is no moving mesh here, we set the past cell volume to the same as 
  // the non moving one
  m_splitter->setPastCellVolume(volumes[ddata.cellID]);

  // reset the residual because we will accumulate the sub element contributions
  for (CFuint i = 0; i < nbStatesInCell; ++i)
  {
    residual[i] = 0.0;
  }

  
//get the interstates in this cell (intermediate layer)
   DataHandle<State*> interStatesStorage = socket_interStates.getDataHandle();
   for (CFuint i = 0; i < nbStatesInCell; ++i) {
     const CFuint stateID = (*ddata.states)[i]->getLocalID();
     State & temp =  *m_interStates[i];
     //We need the local ID when we do Linearized Euler
     //to access the meanflow 
     temp.clone(*(*ddata.states)[i]);
     *m_interStates[i] = *interStatesStorage[stateID];
   }

  if (SubSystemStatusStack::getActive()->isFirstStep()){
  // We compute the advective fluctuation on the present layer and the
  // intermediate one and we store it
  ddata.time = SubSystemStatusStack::getActive()->getCurrentTime();
  computeHOSTFluxIntegralP1(m_interStates,m_phi_inter);
  for (CFuint iSubElm = 0; iSubElm < 4; ++iSubElm){
      for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
        m_interResiduals[ddata.cellID][(iSubElm*m_nbEqs)+iEq] = m_phi_inter[iSubElm][iEq];
    }
  }
}
if (!SubSystemStatusStack::getActive()->isFirstStep()){
//If it is not the first iteration then m_phi_inter has already been computed 
  for (CFuint iSubElm = 0; iSubElm < 4; ++iSubElm){
      for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
         m_phi_inter[iSubElm][iEq] = m_interResiduals[ddata.cellID][(iSubElm*m_nbEqs)+iEq];
    }
  }
}
ddata.time = SubSystemStatusStack::getActive()->getCurrentTime()+dt;
  computeHOSTFluxIntegralP1(states,m_phi_present);
ddata.time = SubSystemStatusStack::getActive()->getCurrentTime();
   /*****         Triangle 0-3-5          *****/

  substates[0] = states[0];
  substates[1] = states[3];
  substates[2] = states[5];

  computeHOtimeP1(0,3,5);

  // linearize the states in the cell
  // includes the transformation from update to linearization
  // variables to evaluate the jacobians in the average state
  // and then do the transformation from update to consistent variables
  
  ddata.tStates = computeConsistentStates(&substates);
  ddata.subStates = &substates;

  for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
    m_flux_space[iEq] = timeStep*(m_phi_inter[0][iEq] + m_phi_present[0][iEq])+m_phi_time[iEq];
  }

  //Need to scale the cell normals to the half
  cellnormals.scale(0.5);
  m_splitter->computeK(substates, &cellnormals);
  //Need to unscale after
  cellnormals.unscale();

  //Transform back the residuals to the distributed variables
  SafePtr<RealVector> phi = &ddata.phi;
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_flux_space);

  // Distribute the residual to the states of the present layer
  m_splitter->distribute(subresidual);

  residual[0] += subresidual[0];
  residual[3] += subresidual[1];
  residual[5] += subresidual[2];

   /*****         Triangle 3-1-4          *****/

  substates[0] = states[3];
  substates[1] = states[1];
  substates[2] = states[4];

  computeHOtimeP1(3,1,4);

  // linearize the states in the cell
  // includes the transformation from update to linearization
  // variables to evaluate the jacobians in the average state
  // and then do the transformation from update to consistent variables
  
  ddata.tStates = computeConsistentStates(&substates);
  ddata.subStates = &substates;

  for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
    m_flux_space[iEq] = timeStep*(m_phi_inter[1][iEq] + m_phi_present[1][iEq])+m_phi_time[iEq];
  }

  //Need to scale the cell normals to the half
  cellnormals.scale(0.5);
  m_splitter->computeK(substates, &cellnormals);
  //Need to unscale after
  cellnormals.unscale();

  //Transform back the residuals to the distributed variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_flux_space);

  // Distribute the residual to the states of the present layer
  m_splitter->distribute(subresidual);

  residual[3] += subresidual[0];
  residual[1] += subresidual[1];
  residual[4] += subresidual[2];
  
 /*****         Triangle 5-4-2          *****/

  substates[0] = states[5];
  substates[1] = states[4];
  substates[2] = states[2];

  computeHOtimeP1(5,4,2);

  // linearize the states in the cell
  // includes the transformation from update to linearization
  // variables to evaluate the jacobians in the average state
  // and then do the transformation from update to consistent variables
  
  ddata.tStates = computeConsistentStates(&substates);
  ddata.subStates = &substates;

  for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
    m_flux_space[iEq] = timeStep*(m_phi_inter[2][iEq] + m_phi_present[2][iEq])+ m_phi_time[iEq];
  }

  //Need to scale the cell normals to the half
  cellnormals.scale(0.5);
  m_splitter->computeK(substates, &cellnormals);
  //Need to unscale after
  cellnormals.unscale();

  //Transform back the residuals to the distributed variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_flux_space);

  // Distribute the residual to the states of the present layer
  m_splitter->distribute(subresidual);

  residual[5] += subresidual[0];
  residual[4] += subresidual[1];
  residual[2] += subresidual[2];

 /*****         Triangle 4-5-3          *****/

  substates[0] = states[4];
  substates[1] = states[5];
  substates[2] = states[3];

  computeHOtimeP1(4,5,3);

  // linearize the states in the cell
  // includes the transformation from update to linearization
  // variables to evaluate the jacobians in the average state
  // and then do the transformation from update to consistent variables
  
  ddata.tStates = computeConsistentStates(&substates);
  ddata.subStates = &substates;

  for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
    m_flux_space[iEq] = timeStep*(m_phi_inter[3][iEq] + m_phi_present[3][iEq])+m_phi_time[iEq];
  }

  //Need to scale the cell normals to the half
  cellnormals.scale(-0.5);
  m_splitter->computeK(substates, &cellnormals);
  //Need to unscale after
  cellnormals.unscale();

  //Transform back the residuals to the distributed variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_flux_space);

  // Distribute the residual to the states of the present layer
  m_splitter->distribute(subresidual);

  residual[4] += subresidual[0];
  residual[5] += subresidual[1];
  residual[3] += subresidual[2];
}
//////////////////////////////////////////////////////////////////////////////
void STU_HOCRD_SplitStrategy::docomputeFluctuation(std::vector<RealVector>& residual){
///@todo Nv: There is no moving mesh implemented
  DistributionData& ddata = getMethodData().getDistributionData();
  vector<State*>& states = *ddata.states;
  const CFuint nbStatesInCell = states.size();

  m_splitter->setDT(SubSystemStatusStack::getActive()->getDT());
//  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal dt1 = SubSystemStatusStack::getActive()->getDT();
  const CFreal dt0 =SubSystemStatusStack::getActive()->getPreviousDT();
  const CFreal q = dt1/dt0;
  const CFreal Jpast  = -dt1*(q*q)/(6.0*(1.0 + q)); 
  const CFreal Jinter  = dt1*(3.0 + q)/6.0;
  const CFreal Jpresent  = dt1*(3.0 + 2.0*q)/(6.0*(1.0 + q)); 

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [ddata.cellID]);
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  CellVolume = volumes[ddata.cellID];
  // backup the states of this cell
  setCurrentCell();

  // Set the cell volume such that it can be accessed by the splitter
  m_splitter->setCellVolume(volumes[ddata.cellID]);

  // reset the residual because we will accumulate the sub element contributions
  for (CFuint i = 0; i < nbStatesInCell; ++i)
  {
    residual[i] = 0.0;
  }

  //get the interstates in this cell (intermediate layer)
   DataHandle<State*> interStatesStorage = socket_interStates.getDataHandle();
   for (CFuint i = 0; i < nbStatesInCell; ++i) {
     const CFuint stateID = (*ddata.states)[i]->getLocalID();
     State & temp =  *m_interStates[i];
    //We need the local ID when we do Linearized Euler
     //to access the meanflow 
     temp.clone(*(*ddata.states)[i]);
     *m_interStates[i] = *interStatesStorage[stateID];
   }

  //get the paststates in this cell (past layer)
  DataHandle<State*> pastStatesStorage = socket_pastStates.getDataHandle();
  for (CFuint i = 0; i < nbStatesInCell; ++i) {
    const CFuint stateID = (*ddata.states)[i]->getLocalID();
    State & temp =  *m_pastStates[i];
    //We need the local ID when we do Linearized Euler
    //to access the meanflow 
    temp.clone(*(*ddata.states)[i]);
    *m_pastStates[i] = *pastStatesStorage[stateID];
  }

  // We compute the advective fluctuation on the present layer and the
  // intermediate one

if (SubSystemStatusStack::getActive()->isFirstStep()){
  // We compute the advective fluctuation on the present layer and the
  // intermediate one and we store it
  ddata.time = SubSystemStatusStack::getActive()->getCurrentTime();
  computeHOSTFluxIntegral(m_interStates,m_phi_inter);
  ddata.time = SubSystemStatusStack::getActive()->getCurrentTime() - dt0;
  computeHOSTFluxIntegral(m_pastStates,m_phi_past);
  ddata.time = SubSystemStatusStack::getActive()->getCurrentTime();
  for (CFuint iSubElm = 0; iSubElm < 4; ++iSubElm){
      for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
        m_interResiduals[ddata.cellID][(iSubElm*m_nbEqs)+iEq] = m_phi_inter[iSubElm][iEq];
        m_pastResiduals[ddata.cellID][(iSubElm*m_nbEqs)+iEq] = m_phi_past[iSubElm][iEq];
    }
  }
}

if (!SubSystemStatusStack::getActive()->isFirstStep()){
//If it is not the first iteration then m_phi_inter has already been computed 
  for (CFuint iSubElm = 0; iSubElm < 4; ++iSubElm){
      for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
         m_phi_inter[iSubElm][iEq] = m_interResiduals[ddata.cellID][(iSubElm*m_nbEqs)+iEq];
         m_phi_past[iSubElm][iEq] = m_pastResiduals[ddata.cellID][(iSubElm*m_nbEqs)+iEq];
    }
  }
}
  ddata.time = SubSystemStatusStack::getActive()->getCurrentTime() + dt1;
  computeHOSTFluxIntegral(states,m_phi_present);
  ddata.time = SubSystemStatusStack::getActive()->getCurrentTime();
   /*****         Triangle 0-3-5          *****/

  substates[0] = states[0];
  substates[1] = states[3];
  substates[2] = states[5];

  computeHOtime(0,3,5);

  // linearize the states in the cell
  // includes the transformation from update to linearization
  // variables to evaluate the jacobians in the average state
  // and then do the transformation from update to consistent variables
  
  ddata.tStates = computeConsistentStates(&substates);
  ddata.subStates = &substates;

  for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
    m_flux_space[iEq] = Jpast*m_phi_past[0][iEq] + Jinter*m_phi_inter[0][iEq] + 
                        Jpresent*m_phi_present[0][iEq]+m_phi_time[iEq];
  }

  //Need to scale the cell normals to the half
  cellnormals.scale(0.5);
  m_splitter->computeK(substates, &cellnormals);
  //Need to unscale after
  cellnormals.unscale();

  //Transform back the residuals to the distributed variables
  SafePtr<RealVector> phi = &ddata.phi;
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_flux_space);

  // Distribute the residual to the states of the present layer
  m_splitter->distribute(subresidual);

  residual[0] += subresidual[0];
  residual[3] += subresidual[1];
  residual[5] += subresidual[2];

   /*****         Triangle 3-1-4          *****/

  substates[0] = states[3];
  substates[1] = states[1];
  substates[2] = states[4];

  computeHOtime(3,1,4);

  // linearize the states in the cell
  // includes the transformation from update to linearization
  // variables to evaluate the jacobians in the average state
  // and then do the transformation from update to consistent variables
  
  ddata.tStates = computeConsistentStates(&substates);
  ddata.subStates = &substates;

  for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
    m_flux_space[iEq] = Jpast*m_phi_past[1][iEq] + Jinter*m_phi_inter[1][iEq] + 
                        Jpresent*m_phi_present[1][iEq] + m_phi_time[iEq];
  }

  //Need to scale the cell normals to the half
  cellnormals.scale(0.5);
  m_splitter->computeK(substates, &cellnormals);
  //Need to unscale after
  cellnormals.unscale();

  //Transform back the residuals to the distributed variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_flux_space);

  // Distribute the residual to the states of the present layer
  m_splitter->distribute(subresidual);

  residual[3] += subresidual[0];
  residual[1] += subresidual[1];
  residual[4] += subresidual[2];
  
 /*****         Triangle 5-4-2          *****/

  substates[0] = states[5];
  substates[1] = states[4];
  substates[2] = states[2];

  computeHOtime(5,4,2);

  // linearize the states in the cell
  // includes the transformation from update to linearization
  // variables to evaluate the jacobians in the average state
  // and then do the transformation from update to consistent variables
  
  ddata.tStates = computeConsistentStates(&substates);
  ddata.subStates = &substates;

  for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
    m_flux_space[iEq] = Jpast*m_phi_past[2][iEq] + Jinter*m_phi_inter[2][iEq] + 
                        Jpresent*m_phi_present[2][iEq] + m_phi_time[iEq];
  }

  //Need to scale the cell normals to the half
  cellnormals.scale(0.5);
  m_splitter->computeK(substates, &cellnormals);
  //Need to unscale after
  cellnormals.unscale();

  //Transform back the residuals to the distributed variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_flux_space);

  // Distribute the residual to the states of the present layer
  m_splitter->distribute(subresidual);

  residual[5] += subresidual[0];
  residual[4] += subresidual[1];
  residual[2] += subresidual[2];

 /*****         Triangle 4-5-3          *****/

  substates[0] = states[4];
  substates[1] = states[5];
  substates[2] = states[3];

  computeHOtime(4,5,3);

  // linearize the states in the cell
  // includes the transformation from update to linearization
  // variables to evaluate the jacobians in the average state
  // and then do the transformation from update to consistent variables
  
  ddata.tStates = computeConsistentStates(&substates);
  ddata.subStates = &substates;

  for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq){
    m_flux_space[iEq] = Jpast*m_phi_past[3][iEq] + Jinter*m_phi_inter[3][iEq] + 
                        Jpresent*m_phi_present[3][iEq] + m_phi_time[iEq];
  }

  //Need to scale the cell normals to the half
  cellnormals.scale(-0.5);
  m_splitter->computeK(substates, &cellnormals);
  //Need to unscale after
  cellnormals.unscale();

  //Transform back the residuals to the distributed variables
  *phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_flux_space);

  // Distribute the residual to the states of the present layer
  m_splitter->distribute(subresidual);

  residual[4] += subresidual[0];
  residual[5] += subresidual[1];
  residual[3] += subresidual[2];

}

//////////////////////////////////////////////////////////////////////////////

void STU_HOCRD_SplitStrategy::computeHOSTFluxIntegral(std::vector<Framework::State*>& states, std::vector<RealVector>& m_phi)
{
  
  cf_assert(qdstates.size() == 3); // only triags so three quadrature points per face
  DistributionData& ddata = getMethodData().getDistributionData();
  CFuint nbCellStates = states.size();
  DataHandle<InwardNormalsData*> normals = socket_normals.getDataHandle();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
ddata.isfirstP1 = false;

  const State& state0 = *(states[0]);
  const State& state1 = *(states[1]);
  const State& state2 = *(states[2]);
  const State& state3 = *(states[3]);
  const State& state4 = *(states[4]);
  const State& state5 = *(states[5]);

  vector<Node*>& nodes = *ddata.cell->getNodes();
  cf_assert(nodes.size()  == 3); // P1 triangles for geometry space

  const CFreal x1 = (*nodes[0])[XX];
  const CFreal x2 = (*nodes[1])[XX];
  const CFreal x3 = (*nodes[2])[XX];

  const CFreal y1 = (*nodes[0])[YY];
  const CFreal y2 = (*nodes[1])[YY];
  const CFreal y3 = (*nodes[2])[YY];
  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [ddata.cellID]);

  const CFreal nx1 = cellnormals.getNodalNormComp(0,XX);
  const CFreal nx2 = cellnormals.getNodalNormComp(1,XX);
  const CFreal nx3 = cellnormals.getNodalNormComp(2,XX);

  const CFreal ny1 = cellnormals.getNodalNormComp(0,YY);
  const CFreal ny2 = cellnormals.getNodalNormComp(1,YY);
  const CFreal ny3 = cellnormals.getNodalNormComp(2,YY);
  const CFreal inv_volume = 1.0 / CellVolume;

  CFuint nbFace = subfacetable.nbRows();
  //Computation  of the fluctuation of each faces of the cell
  for (CFuint iFace = 0; iFace < nbFace; ++iFace)
    {
      facenormal[XX] = 0.5*cellnormals.getNodalNormComp(iFace%3,XX);
      facenormal[YY] = 0.5*cellnormals.getNodalNormComp(iFace%3,YY);

      for (CFuint iQd = 0; iQd < 3; ++iQd)
        {
          const Node& node0 = states[subfacetable(iFace,0)]->getCoordinates();
          const Node& node1 = states[subfacetable(iFace,1)]->getCoordinates();
          const CFreal x = qd0a[iQd] * node0[XX] + qd1a[iQd] * node1[XX];
          const CFreal y = qd0a[iQd] * node0[YY] + qd1a[iQd] * node1[YY];
          const CFreal L1 = 1.0 + 0.5*( ( x - x1 )*nx1 + ( y - y1 )*ny1 ) * inv_volume ;
          const CFreal L2 = 1.0 + 0.5*( ( x - x2 )*nx2 + ( y - y2 )*ny2 ) * inv_volume ;
          const CFreal L3 = 1.0 + 0.5*( ( x - x3 )*nx3 + ( y - y3 )*ny3 ) * inv_volume ;

          (*qdstates[iQd]) = (L1*( 2.0*L1 - 1.0 ) * state0) +
                             (L2*( 2.0*L2 - 1.0 ) * state1) +
                             (L3*( 2.0*L3 - 1.0 ) * state2) +
                             (4.0*L1*L2           * state3) +
                             (4.0*L3*L2           * state4) +
                             (4.0*L1*L3           * state5);

	 (*qdnodes[iQd])[XX] = x;
  	 (*qdnodes[iQd])[YY] = y;

    }
    computeStatesData(3, m_updateVar, qdstates, m_pdata, qdExtraVars); // three quadrature points per face
    faceflux[iFace] = 0.;
    for (CFuint iQd = 0; iQd < 3; ++iQd)
      {
        faceflux[iFace] += wqda[iQd] * m_updateVar->getFlux()(m_pdata[iQd],facenormal);
      }
    } 
  CFuint subelemtab_nbRows = subelemtable.nbRows();
  CFuint subelemtab_nbCols = subelemtable.nbCols();

  cf_assert(m_phi.size() == subelemtab_nbRows);


  for (CFuint iState = 0; iState < nbCellStates; ++iState) {
    ddata.cell->setState(iState, states[iState]);
  }

  // recomposition of the fluctuation of each sub cell (taking care of the orientation)
  for (CFuint iRow = 0; iRow < subelemtab_nbRows; ++iRow)
  {
    RealVector& phi = (m_phi[iRow]);
    phi = 0.;
    for (CFuint jCol = 0; jCol < subelemtab_nbCols; ++jCol)
    {
      phi -= subelemfacedir(iRow,jCol) * faceflux[subelemtable(iRow,jCol)];
    }

   for (CFuint jCol = 0; jCol < subelemtab_nbCols; ++jCol)
    {
    substates[jCol] = states[subelemtable_state(iRow,jCol)];
    }   
    ddata.subStates = &substates;
  
   const CFuint nbST = getMethodData().getSourceTermSplitter()->size();
    for (CFuint i = 0; i < nbST; ++i)
    {
      ddata.sourceTermID = i;
      // computes the source term and places it in ddata.phiS
      // To check : that the state used is indeed the past state 
      // For monopole, dipole and quadripole it is not making any difference
      // because they only depend on the time, but for others could be a problem
      getMethodData().getSourceTermSplitter(i)->computeSourceTerm(cellnormals);

      phi -= ddata.phiS;
    }
    
  
  }

  for (CFuint iState = 0; iState < nbCellStates; ++iState) {
    ddata.cell->setState(iState, m_statesBkp[iState]);
  }

}

//////////////////////////////////////////////////////////////////////////////

void STU_HOCRD_SplitStrategy::computeHOtime(CFuint i1, CFuint i2, CFuint i3){
  DistributionData& ddata = getMethodData().getDistributionData();
  vector<State*>& states = *ddata.states;

  // Coordinate of the vertex of the element
  const Node& node0_vert = states[0]->getCoordinates();
  const Node& node1_vert = states[1]->getCoordinates();
  const Node& node2_vert = states[2]->getCoordinates();
  
  const CFreal x1_vert = node0_vert[XX];
  const CFreal x2_vert = node1_vert[XX];
  const CFreal x3_vert = node2_vert[XX];

  const CFreal y1_vert = node0_vert[YY];
  const CFreal y2_vert = node1_vert[YY];
  const CFreal y3_vert = node2_vert[YY];

  //Coordinates of the node defining the surface where the integral is computed
  //  This correspond to vertex of the sub-element on which we integrate
  Node& node0 = states[i1]->getCoordinates();
  Node& node1 = states[i2]->getCoordinates();
  Node& node2 = states[i3]->getCoordinates();

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()[ddata.cellID]);

  const CFreal nx1 = cellnormals.getNodalNormComp(0,XX);
  const CFreal nx2 = cellnormals.getNodalNormComp(1,XX);
  const CFreal nx3 = cellnormals.getNodalNormComp(2,XX);

  const CFreal ny1 = cellnormals.getNodalNormComp(0,YY);
  const CFreal ny2 = cellnormals.getNodalNormComp(1,YY);
  const CFreal ny3 = cellnormals.getNodalNormComp(2,YY);

  const CFreal inv_volume = 1.0/CellVolume;

   *temp_states[0] = *states[0] - *m_interStates[0];
   *temp_states[1] = *states[1] - *m_interStates[1];
   *temp_states[2] = *states[2] - *m_interStates[2];
   *temp_states[3] = *states[3] - *m_interStates[3];
   *temp_states[4] = *states[4] - *m_interStates[4];
   *temp_states[5] = *states[5] - *m_interStates[5];
   
  m_phi_time = 0.0;
  for (CFuint iQd = 0; iQd < 4; ++iQd) {
     //point od quadrature
     const CFreal x = qd0t[iQd] * node0[XX] + qd1t[iQd] * node1[XX] + qd2t[iQd] * node2[XX] ;
     const CFreal y = qd0t[iQd] * node0[YY] + qd1t[iQd] * node1[YY] + qd2t[iQd] * node2[YY] ;

     // Linear basis function
     CFreal L1 = 1.0 + 0.5*( ( x - x1_vert )*nx1 + ( y - y1_vert )*ny1 )*inv_volume ;
     CFreal L2 = 1.0 + 0.5*( ( x - x2_vert )*nx2 + ( y - y2_vert )*ny2 )*inv_volume ;
     CFreal L3 = 1.0 + 0.5*( ( x - x3_vert )*nx3 + ( y - y3_vert )*ny3 )*inv_volume ;

         m_u = (L1*(2.0*L1 - 1.0)*(*(temp_states[0])) +
	      L2*(2.0*L2 - 1.0)*(*(temp_states[1])) +
	      L3*(2.0*L3 - 1.0)*(*(temp_states[2])) +
	      4.0*L1*L2*(*(temp_states[3])) +
	      4.0*L2*L3*(*(temp_states[4])) +
	      4.0*L3*L1*(*(temp_states[5])));
  
        m_phi_time += m_u*CellVolume*wqdt[iQd];
        


  }

}


//////////////////////////////////////////////////////////////////////////////

void STU_HOCRD_SplitStrategy::computeHOtimeP1(CFuint i1, CFuint i2, CFuint i3){
  DistributionData& ddata = getMethodData().getDistributionData();
  vector<State*>& states = *ddata.states;


   *temp_states[0] = *states[i1] - *m_interStates[i1];
   *temp_states[1] = *states[i2] - *m_interStates[i2];
   *temp_states[2] = *states[i3] - *m_interStates[i3];
    
  m_phi_time = (CellVolume/12.0)*(*temp_states[0] + *temp_states[1] + *temp_states[2]) ;

}

//////////////////////////////////////////////////////////////////////////////

void STU_HOCRD_SplitStrategy::computeHOSTFluxIntegralP1(std::vector<Framework::State*>& states, std::vector<RealVector>& m_phi)
{
  
  cf_assert(qdstates.size() == 3); // only triags so three quadrature points per face
  DistributionData& ddata = getMethodData().getDistributionData();

  vector<Node*>& nodes = *ddata.cell->getNodes();

  cf_assert(nodes.size()  == 3); // P1 triangles for geometry space

  InwardNormalsData& cellnormals = (*socket_normals.getDataHandle()
				    [ddata.cellID]);

 
  CFuint nbFace = subfacetable.nbRows();
  //Computation  of the fluctuation of each faces of the cell
  for (CFuint iFace = 0; iFace < nbFace; ++iFace)
    {
      facenormal[XX] = 0.5*cellnormals.getNodalNormComp(iFace%3,XX);
      facenormal[YY] = 0.5*cellnormals.getNodalNormComp(iFace%3,YY);


      for (CFuint iQd = 0; iQd < 3; ++iQd)
        {
          const Node& node0 = states[subfacetable(iFace,0)]->getCoordinates();
          const Node& node1 = states[subfacetable(iFace,1)]->getCoordinates();

          const CFreal x = qd0a[iQd] * node0[XX] + qd1a[iQd] * node1[XX];
          const CFreal y = qd0a[iQd] * node0[YY] + qd1a[iQd] * node1[YY];
         // const CFreal L1 = (x - node0[XX])*;
          //const CFreal L2 = 1.0 + 0.5*( ( x - node1[XX] )*nx2 + ( y - node1[YY] )*ny2 ) * inv_volume ;

          const State& state0 = *(states[subfacetable(iFace,0)]);
          const State& state1 = *(states[subfacetable(iFace,1)]);

          (*qdstates[iQd]) = qd0a[iQd]* state0 + qd1a[iQd]* state1 ;

	 (*qdnodes[iQd])[XX] = x;
  	 (*qdnodes[iQd])[YY] = y;

    }

    computeStatesData(3, m_updateVar, qdstates, m_pdata, qdExtraVars); // three quadrature points per face
    
    faceflux[iFace] = 0.;

    for (CFuint iQd = 0; iQd < 3; ++iQd)
      {

        faceflux[iFace] += wqda[iQd] * m_updateVar->getFlux()(m_pdata[iQd],facenormal);

      }
    } 

  CFuint subelemtab_nbRows = subelemtable.nbRows();
  CFuint subelemtab_nbCols = subelemtable.nbCols();

  cf_assert(m_phi.size() == subelemtab_nbRows);

  // recomposition of the fluctuation of each sub cell (taking care of the orientation)
  for (CFuint iRow = 0; iRow < subelemtab_nbRows; ++iRow)
  {
    RealVector& phi = (m_phi[iRow]);
    phi = 0.;
    for (CFuint jCol = 0; jCol < subelemtab_nbCols; ++jCol)
    {
      phi -= subelemfacedir(iRow,jCol) * faceflux[subelemtable(iRow,jCol)];
    }

 for (CFuint jCol = 0; jCol < subelemtab_nbCols; ++jCol)
    {
    substates[jCol] = states[subelemtable_state(iRow,jCol)];
    }
    ddata.subStates = &substates;
   const CFuint nbST = getMethodData().getSourceTermSplitter()->size();
    for (CFuint i = 0; i < nbST; ++i)
    {
      ddata.sourceTermID = i;
      // computes the source term and places it in ddata.phiS
      // To check : that the state used is indeed the past state
      // For monopole, dipole and quadripole it is not making any difference
      // because they only depend on the time, but for others could be a problem
      getMethodData().getSourceTermSplitter(i)->computeSourceTerm(cellnormals);
      phi -= ddata.phiS;
     }
     
  }
}


//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
STU_HOCRD_SplitStrategy::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = FluctuationSplitStrategy::needsSockets();
  result.push_back(&socket_updateCoeff);
   result.push_back(&socket_interStates);
    result.push_back(&socket_pastStates);
  return result;  
}

//////////////////////////////////////////////////////////////////////////////

void STU_HOCRD_SplitStrategy::setCurrentCell(){

  DistributionData& ddata = getMethodData().getDistributionData();  
  // back up the update states and sets them in the linearizer
  // AL: don't touch this! update states have to be backed up
  // because linearizing states are temporarily stored 
  // in the cell and used to contour integrate
  const CFuint nbCellStates = ddata.states->size();

  for (CFuint iState = 0; iState < nbCellStates; ++iState) {

    m_statesBkp[iState] = (*ddata.states)[iState]; 
  } 
  getMethodData().getLinearizer()->setUpdateStates(&m_statesBkp);
   
  ddata.tStates = computeConsistentStates(ddata.states);
}

//////////////////////////////////////////////////////////////////////////////

void STU_HOCRD_SplitStrategy::unsetup()
{
  // AL: this could be changed ... 
  // where are deleted the qdstates ???
  
  for (CFuint i = 0; i < qdExtraVars.size(); ++i) {
    deletePtr(qdExtraVars[i]);
  }


  FluctuationSplitStrategy::unsetup();
}
 
//////////////////////////////////////////////////////////////////////////////

   } // namespace FluctSplit



} // namespace COOLFluiD
