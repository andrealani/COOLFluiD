#include "Framework/SubSystemStatus.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/BadValueException.hh"
#include "Common/PE.hh"
#include "Framework/PathAppender.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolume/DerivativeComputer.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "AeroCoef/AeroCoefFVM.hh"
#include "AeroCoef/NavierStokesSkinFrictionHeatFluxCC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NavierStokesSkinFrictionHeatFluxCC,
		      DataProcessingData,
		      AeroCoefFVMModule>
navierStokesSkinFrictionHeatFluxCCProvider("NavierStokesSkinFrictionHeatFluxCC");

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxCC::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("StantonNumberID","ID identifying definition of Stanton number");
}
    
//////////////////////////////////////////////////////////////////////////////

NavierStokesSkinFrictionHeatFluxCC::NavierStokesSkinFrictionHeatFluxCC(const std::string& name) :
  AeroForcesFVMCC(name),
  socket_wallDistance("wallDistance",false),
  _derivComputer(CFNULL),
  _diffVar(CFNULL),
  _qradFluxWall(CFNULL),
  _hasRadiationCoupling(false),
  _states(),
  _values(),
  _gradients(),
  _rhoWall(0.),
  _muWall(0.),
  _yPlus(0.),
  _tau(0.),
  _heatFluxRad(0.),
  _tau3D(DIM_3D,DIM_3D),
  m_Cf3D(DIM_3D,DIM_3D)
{
  addConfigOptionsTo(this);
  
  _stantonNumID = 0;
  setParameter("StantonNumberID",&_stantonNumID);
}

//////////////////////////////////////////////////////////////////////////////

NavierStokesSkinFrictionHeatFluxCC::~NavierStokesSkinFrictionHeatFluxCC()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
NavierStokesSkinFrictionHeatFluxCC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = 
    AeroForcesFVMCC::needsSockets();
  result.push_back(&socket_wallDistance);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxCC::setup()
{
  CFAUTOTRACE;
  
  AeroForcesFVMCC::setup();
  
  _derivComputer = m_fvmccData->getDerivativeComputer();
  
  _diffVar = m_fvmccData->getDiffusiveVar().d_castTo<NavierStokesVarSet>();
  cf_assert(_diffVar.isNotNull());
  cf_assert(m_updateVarSet.isNotNull());
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint maxNbNodesInControlVolume =
    _derivComputer->getMaxNbVerticesInControlVolume();
  
  _states.resize(maxNbNodesInControlVolume);
  _values.resize(nbEqs, maxNbNodesInControlVolume);
  
  _gradients.resize(nbEqs);
  for (CFuint i = 0; i< nbEqs; ++i) {
    _gradients[i] = new RealVector(PhysicalModelStack::getActive()->getDim());
  }

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  if (m_TID == 0) {
    m_TID = (dim == DIM_3D) ? 4 : 3;
  }
  
  if (m_UID == 0) m_UID = 1;
  if (m_VID == 0) m_VID = 2;
  
  if (m_WID == 0) {
    m_WID = (dim == DIM_3D) ? 3 : 0;
  }
  
  m_varNames.clear();
  for (CFuint i = 0; i < dim; ++i) {
    const std::string xdim = "x" + Common::StringOps::to_str(i);
    this->m_varNames.push_back(xdim);
  }
  m_varNames.push_back("P");
  m_varNames.push_back("T");
  m_varNames.push_back("rho");
  m_varNames.push_back("Cp");
  m_varNames.push_back("heatF");
  m_varNames.push_back("Stanton");
  m_varNames.push_back("yplus");
  m_varNames.push_back("Cf");
  m_varNames.push_back("muWall"); 
  m_varNames.push_back("heatFRad");
  
  cf_always_assert(this->m_varNames.size() == 10 + dim); 
  
  // check if the radiative heat is stored
  const string qradName = MeshDataStack::getActive()->getPrimaryNamespace() + "_qradFluxWall";
  _hasRadiationCoupling = MeshDataStack::getActive()->getDataStorage()->checkData(qradName);
  if (_hasRadiationCoupling) {
    _qradFluxWall = MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(qradName);
  }  
}

//////////////////////////////////////////////////////////////////////////////
    
void NavierStokesSkinFrictionHeatFluxCC::unsetup()
{
  for (CFuint i = 0; i< _gradients.size(); ++i) {
    deletePtr(_gradients[i]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxCC::computeWall()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle< RealVector> nstates = socket_nstates.getDataHandle();
  DataHandle< CFreal> normals = socket_normals.getDataHandle();
  
  cf_always_assert(PhysicalModelStack::getActive()->getDim() == DIM_2D);
  const bool isPerturb = false;
  const CFuint dummyVarID = 0;
  
  // set the state values (pointers) corresponding to the
  // vertices of the control volume
  _derivComputer->computeControlVolume(_states, m_currFace);
  
  // compute speed components and temperature in the given states
  // if you are updating in conservative variables your nodal values
  // MUST be already in primitive variables (there is inconsistency here !!!)
  const CFuint nbCVStates = _derivComputer->getNbVerticesInControlVolume(m_currFace);
  _diffVar->setGradientVars(_states, _values, nbCVStates);
  
  // compute control volume around the face and gradients
  _derivComputer->computeGradients(m_currFace, _values, _gradients);
  
  // compute the average values [p u v w T]
  _derivComputer->computeAverageValues(m_currFace, _states, (*m_avState));
  
  // this is needed for LTE
  _diffVar->setComposition((*m_avState), isPerturb, dummyVarID);
  
  // Compute dynamic viscosity and density
  // We are computing at the wall where the wall distance is NULL!
  _diffVar->setWallDistance(0.);
  _muWall = _diffVar->getDynViscosity((*m_avState), _gradients);
  _rhoWall = _diffVar->getDensity((*m_avState));
  
  // Compute the friction at the wall
  computeTauWall();
  
  // Compute y+ value
  if(socket_wallDistance.isConnected())  {
    computeYplus();
  }
  
  updateWriteData();
    
  // Do some extra computation on the face if needed
  computeExtraValues();
}
    
//////////////////////////////////////////////////////////////////////////////
    
void NavierStokesSkinFrictionHeatFluxCC::computeTauWall()
{
  // const CFreal Uref = m_updateVarSet->getModel()->getVelRef();
  //const CFreal refLength = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  //const CFreal mu_ref = _diffVar->getModel().getDynViscosityDim
  //(m_updateVarSet->getModel()->getPressRef(), m_updateVarSet->getModel()->getTempRef());
  
  // this will not work adimensional
  // const CFreal tauRef = 1.0; // mu_ref*Uref/refLength;
  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
  const bool adim = Framework::PhysicalModelStack::getActive()->getImplementor()->isAdimensional();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  _tau = 0.0;
  m_Cf = 0.0;
  
  if (dim == DIM_2D) {
    
    _tau = _muWall*
      (m_unitNormal[YY]*MathFunctions::innerProd(*_gradients[m_UID],m_unitNormal) -
       m_unitNormal[XX]*MathFunctions::innerProd(*_gradients[m_VID],m_unitNormal));
    
    if(adim){
      //check this
      const CFreal mInf = m_updateVarSet->getModel()->getMachInf();
      const CFreal Re = _diffVar->getModel().getReynolds();
      m_Cf = _tau / (0.5*mInf*sqrt(gamma)*Re);
    }
    else{
      // friction coefficient 
      m_Cf = _tau / (0.5*m_rhoInf*m_uInf*m_uInf);
    } 
    
    m_frictionForces[XX] =  m_Cf*m_normal[YY];
    m_frictionForces[YY] = -m_Cf*m_normal[XX];
  }
  else{
    cf_assert(dim == DIM_3D);
    
    CFreal divU = (2.0/3.0)*((*_gradients[m_UID])[0] + (*_gradients[m_VID])[1] + (*_gradients[m_WID])[2]) ;
    
    _tau3D(XX,XX) = _muWall*(2.0*(*_gradients[m_UID])[XX] - divU) ;
    _tau3D(YY,YY) = _muWall*(2.0*(*_gradients[m_VID])[YY] - divU) ;
    _tau3D(ZZ,ZZ) = _muWall*(2.0*(*_gradients[m_WID])[ZZ] - divU) ;
    
    _tau3D(XX,YY) = _tau3D(YY,XX) = _muWall*((*_gradients[m_UID])[YY] + (*_gradients[m_VID])[XX]) ;
    _tau3D(XX,ZZ) = _tau3D(ZZ,XX) = _muWall*((*_gradients[m_UID])[ZZ] + (*_gradients[m_WID])[XX]) ;
    _tau3D(YY,ZZ) = _tau3D(ZZ,YY) = _muWall*((*_gradients[m_VID])[ZZ] + (*_gradients[m_WID])[YY]) ;
      
    if(adim){
      // check this
      const CFreal mInf = m_updateVarSet->getModel()->getMachInf();
      const CFreal Re = _diffVar->getModel().getReynolds();
      m_Cf3D = _tau3D / (0.5*mInf*sqrt(gamma)*Re);
    }
    else{
      CFreal rhoInf = m_pInf / (m_updateVarSet->getModel()->getR() * m_TInf);
      m_Cf3D = _tau3D / (0.5*rhoInf*m_uInf*m_uInf);
    }
    
    //Compute the viscous force coefficients on this wall face.
    m_frictionForces = m_Cf3D * m_normal;
  }
  
  // forces coefficients must be still divided by the wet surface
  m_frictionForces /= m_refArea;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxCC::prepareOutputFileWall()
{ 
  const std::string nsp = getMethodData().getNamespace();
  PE::GetPE().setBarrier(nsp);
  
  // only the first processor writes the header of the output file 
  if (PE::GetPE().GetRank (nsp) == 0) {
    SafePtr<TopologicalRegionSet> currTrs = this->getCurrentTRS();
    boost::filesystem::path file = Environment::DirPaths::getInstance().getResultsDir() /
      boost::filesystem::path(this->m_nameOutputFileWall + currTrs->getName());
    file = Framework::PathAppender::getInstance().appendAllInfo  
      (file,this->m_appendIter,this->m_appendTime,false);   
    
    SelfRegistPtr<Environment::FileHandlerOutput> fhandle =
      Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    // append to the existing file 
    ofstream& fout = fhandle->open(file);
    
    fout << "TITLE = Unstructured Surface Quantities" << "\n";
    fout << "VARIABLES = "; 
    for (CFuint i = 0; i < this->m_varNames.size(); ++i) {
      fout << this->m_varNames[i] << " ";
    }
    fout << "\n";
    fout.close();
  } 
  
  PE::GetPE().setBarrier(nsp);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxCC::updateOutputFileWall()
{  
  const std::string nsp = getMethodData().getNamespace();
  PE::GetPE().setBarrier(nsp);
  
  // all processors will write their own data one after the other 
  for (CFuint i = 0; i < PE::GetPE().GetProcessorCount(nsp); ++i) {
    if (i == PE::GetPE().GetRank (nsp)) {
      if (getCurrentTRS()->getLocalNbGeoEnts() > 0) {
	SafePtr<TopologicalRegionSet> currTrs = getCurrentTRS();
	boost::filesystem::path file = Environment::DirPaths::getInstance().getResultsDir() /
	  boost::filesystem::path(m_nameOutputFileWall + currTrs->getName());
	file = Framework::PathAppender::getInstance().appendAllInfo  
	  (file,this->m_appendIter,this->m_appendTime,false);   
     
	SelfRegistPtr<Environment::FileHandlerOutput> fhandle =
	  Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
	ofstream& fout = fhandle->open(file, ios::app);
	
	Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
	  geoBuilder = m_fvmccData->getFaceTrsGeoBuilder();
	
	SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
	geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
	
	FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
	geoData.trs = currTrs;
	geoData.isBFace = true;
	
	const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
	for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
	  // build the GeometricEntity
	  geoData.idx = iFace;
	  m_currFace = geoBuilder->buildGE();
	  m_fvmccData->getCurrentFace() = m_currFace;
	  
	  // only faces whose internal State is parallel updatable will write
	  // their data to avoid redudance due to overlap 
	  if (m_currFace->getState(0)->isParUpdatable()) {
	    // compute the face normal
	    const vector<Node*>& faceNodes = *m_currFace->getNodes();
	    const CFuint nbFaceNodes = faceNodes.size();
	    
	    // compute the face mid point
	    m_coord = 0.0;
	    for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode) {
	      m_coord += *faceNodes[iNode];
	    }
	    m_coord /= nbFaceNodes;
	    
	    fout << m_coord << " ";
	   
	    const CFuint index = m_mapTrsFaceToID.find(m_currFace->getID());
	    const CFuint nbVars = m_valuesMat.nbRows(); 
	    for (CFuint iVar = 0; iVar < nbVars; ++iVar) {
	      fout << m_valuesMat(iVar, index) << " ";
	    }
	   fout << "\n";
	  }
	  
	  geoBuilder->releaseGE();
	}
      }
    }
    
    PE::GetPE().setBarrier(nsp);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxCC::updateWriteData() // Vastalya: this is used to write wall data
{  
  const CFreal refLength = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  const CFreal heatFluxRef = m_updateVarSet->getModel()->getTempRef()/refLength;
  const CFuint index = m_mapTrsFaceToID.find(m_currFace->getID());
  
  // compute the radiative heat flux
  _heatFluxRad = 0.;
  if (_hasRadiationCoupling) {
    // AL: the ordering of the wall TRS's must be consistent with the one in RadiativeTransferMonteCarlo
    cf_assert(_qradFluxWall.size() > 0);
    cf_assert(index < _qradFluxWall.size());
    _heatFluxRad = _qradFluxWall[index]*heatFluxRef;
  } 
  
  const CFreal heatFlux = _diffVar->getHeatFlux( (*m_avState), _gradients, m_unitNormal)*heatFluxRef + _heatFluxRad;
  
  CFreal pDim = 0.;
  CFreal rhoDim = 0.;
  CFreal TDim = 0.;
  computeDimensionalPressDensTemp(pDim, rhoDim, TDim);
  
  const CFreal rhoInf = m_pInf / (m_updateVarSet->getModel()->getRdim() * m_TInf);
  
  CFreal stantonNumber = 0.0;
  switch(_stantonNumID) {
  case(0):
    stantonNumber = heatFlux/(rhoInf*pow(m_uInf,3.0));
    break;
  case(1):
    stantonNumber = heatFlux / ((m_updateVarSet->getModel()->getCp()*(m_TInf - TDim) + 
				 0.5* m_uInf*m_uInf )*rhoInf*m_uInf);
    break;
  case(2):
    stantonNumber = heatFlux / (m_updateVarSet->getModel()->getCp()*_rhoWall*m_uInf);
    break;
  default:
    //stantonNumber = heatFlux/(rhoInf*pow(m_uInf,3.0));
    stantonNumber = heatFlux/(m_updateVarSet->getModel()->getCp()*(m_TInf - TDim)*rhoInf*m_uInf);
    break;
  }
  //stantonNumber =  -6; //Vatsalya: just for trial //  it works
  CFreal Cp = (pDim - m_pInf);
  Cp /= (0.5*rhoInf*m_uInf*m_uInf);
  
  // fill in the 2D-array with all the data to be output 
  if (PhysicalModelStack::getActive()->getDim() == DIM_2D){
    updateValuesMatAndResidual(0, index, pDim);    
    updateValuesMatAndResidual(1, index, TDim);    
    updateValuesMatAndResidual(2, index, rhoDim);    
    updateValuesMatAndResidual(3, index, Cp);    
    updateValuesMatAndResidual(4, index, -heatFlux);    
    updateValuesMatAndResidual(5, index, -stantonNumber);
    updateValuesMatAndResidual(6, index, this->_yPlus);
    updateValuesMatAndResidual(7, index, this->m_Cf);
    updateValuesMatAndResidual(8, index, _muWall);
    updateValuesMatAndResidual(9, index, -_heatFluxRad);
  }
  else{
    throw Common::NotImplementedException (FromHere(),"NavierStokesSkinFrictionHeatFluxCC::updateOutputFile() is only 2D");
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxCC::computeYplus()
{
  CFAUTOTRACE;

  _yPlus = 0.;

  DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();

  // Compute the distance to the first cell: y0
  //First get the inner state
  State* innerState = m_currFace->getState(0);
  cf_assert(!innerState->isGhost());

  CFreal y0 = wallDistance[innerState->getLocalID()];

  // Compute the y+
  // y^{+} = \frac {\sqrt{\rho} * \sqrt{\tau} * y0} {\mu}

  const CFreal refSpeed = (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::V];
  const CFreal rhoRef = (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::RHO];

  _yPlus = sqrt(_rhoWall) * sqrt(fabs(_tau)) * y0 / _muWall;
  _yPlus *= sqrt(rhoRef) * refSpeed;

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxCC::computeDimensionalPressDensTemp
(CFreal& pDim, CFreal& rhoDim, CFreal& TDim)
{ 
  // output the data
  const CFuint PID = 0;
  const CFreal rhoRef = (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::RHO];
  
  pDim = m_updateVarSet->getModel()->getPressureFromState((*m_avState)[PID]) * 
    (m_updateVarSet->getModel()->getPressRef());
  TDim = (*m_avState)[m_TID] * (m_updateVarSet->getModel()->getTempRef());
  rhoDim = _rhoWall * rhoRef;
}
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////




