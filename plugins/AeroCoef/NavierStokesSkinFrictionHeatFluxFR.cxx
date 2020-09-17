#include "Framework/SubSystemStatus.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/BadValueException.hh"
#include "Common/PE.hh"
#include "Framework/PathAppender.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FiniteVolume/DerivativeComputer.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "AeroCoef/AeroCoefFR.hh"
#include "AeroCoef/NavierStokesSkinFrictionHeatFluxFR.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::FluxReconstructionMethod;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NavierStokesSkinFrictionHeatFluxFR,
		      DataProcessingData,
		      AeroCoefFRModule>
navierStokesSkinFrictionHeatFluxFRProvider("NavierStokesSkinFrictionHeatFluxFR");

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxFR::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("StantonNumberID","ID identifying definition of Stanton number");
}
    
//////////////////////////////////////////////////////////////////////////////

NavierStokesSkinFrictionHeatFluxFR::NavierStokesSkinFrictionHeatFluxFR(const std::string& name) :
  AeroForcesFR(name),
  socket_wallDistance("wallDistance",false),
  m_diffVar(CFNULL),
  m_qradFluxWall(CFNULL),
  m_hasRadiationCoupling(false),
  m_rhoWall(0.),
  m_muWall(0.),
  m_yPlus(0.),
  m_tau(0.),
  m_heatFluxRad(0.),
  m_tau3D(DIM_3D,DIM_3D),
  m_Cf3D(DIM_3D,DIM_3D)
{
  addConfigOptionsTo(this);
  
  m_stantonNumID = 1;
  setParameter("StantonNumberID",&m_stantonNumID);
}

//////////////////////////////////////////////////////////////////////////////

NavierStokesSkinFrictionHeatFluxFR::~NavierStokesSkinFrictionHeatFluxFR()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
NavierStokesSkinFrictionHeatFluxFR::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = 
    AeroForcesFR::needsSockets();
  result.push_back(&socket_wallDistance);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxFR::setup()
{
  CFAUTOTRACE;
  
  AeroForcesFR::setup();
  
//   m_derivComputer = m_fvmccData->getDerivativeComputer();
  
  m_diffVar = m_frData->getDiffusiveVar().d_castTo<NavierStokesVarSet>();
  cf_assert(m_diffVar.isNotNull());
  cf_assert(m_updateVarSet.isNotNull());
  
  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = m_frData->getFRLocalData();
  cf_assert(frLocalData.size() > 0);

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
  m_varNames.push_back("Cfcrit");
  m_varNames.push_back("Cf");
  m_varNames.push_back("muWall"); 
  m_varNames.push_back("gamma");
  
  cf_always_assert(this->m_varNames.size() == 10 + dim); 
  
  // check if the radiative heat is stored
  const string qradName = MeshDataStack::getActive()->getPrimaryNamespace() + "m_qradFluxWall";
  m_hasRadiationCoupling = MeshDataStack::getActive()->getDataStorage()->checkData(qradName);
  if (m_hasRadiationCoupling) {
    m_qradFluxWall = MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(qradName);
  }  
}

//////////////////////////////////////////////////////////////////////////////
    
void NavierStokesSkinFrictionHeatFluxFR::unsetup()
{
  AeroForcesFR::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxFR::computeWall()
{
  CFAUTOTRACE;
  
  const CFreal R = m_updateVarSet->getModel()->getR();
  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  // loop over flx pnts
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    // compute coordinates of output point
    const RealVector coord = m_currFace->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);

    // dereference current state
    const RealVector& state = *(m_cellStatesFlxPnt[iFlx]);

    // dereference current unit normal
    const RealVector& normal = m_unitNormalFlxPnts[iFlx];
    
    const bool isPerturb = false;
    const CFuint dummyVarID = 0;
  
    // this is needed for LTE
    m_diffVar->setComposition(state, isPerturb, dummyVarID);

    // pressure
    m_rhoWall = m_diffVar->getDensity(state);
//     const CFreal invRho = 1./m_rhoWall;
//     CFreal rhoK2 = 0.0;
//     for (CFuint iDim = 0; iDim < m_dim; ++iDim)
//     {
//       rhoK2 += state[iDim+1]*state[iDim+1];
//     }
//     rhoK2 *= 0.5*invRho;
//     const CFreal p = gammaMinus1*(state[m_nbrEqs-1] - rhoK2);
// 
//     // temperature
//     const CFreal T = p*invRho/R;
// 
//     // dimensional values
//     const CFreal TDim   = T   * m_updateVarSet->getModel()->getTempRef();
//     const CFreal pDim   = p   * m_updateVarSet->getModel()->getPressRef();
//     const CFreal rhoDim = m_rhoWall * (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::RHO];
// 
//     // pressure coefficient
//     m_Cp = (pDim - m_pInf) / (0.5*m_rhoInf*m_uInf*m_uInf);
// 
//     // compute dynamic viscosity
//     m_diffVar->setWallDistance(0.); // we are at the wall
//     m_muWall  = m_diffVar->getDynViscosity(state,m_cellGradFlxPnt[iFlx]);
// 
//     // compute the friction at the wall
//     const CFreal refU = m_updateVarSet->getModel()->getVelRef();
//     const CFreal refLength = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
//     const CFreal refMu = (m_diffVar->getModel().getReferencePhysicalData())[NSTerm::MU];
//     const CFreal refTau = refMu*refU/refLength;
//     const bool adim = Framework::PhysicalModelStack::getActive()->getImplementor()->isAdimensional();
// 
//     // compute shear stress and skin friction tensors
//     RealMatrix tauTensor(m_dim,m_dim);
//     RealMatrix skinFrictionTensor(m_dim,m_dim);
//     RealVector frictionForce(m_dim);
//     CFreal skinFriction = 0.0;
//     if (m_dim == DIM_2D)
//     {
//       const RealVector& gradU = *(m_cellGradFlxPnt[iFlx][1]);
//       const RealVector& gradV = *(m_cellGradFlxPnt[iFlx][2]);
//       const CFreal divU = (2.0/3.0)*(gradU[XX] + gradV[YY]);
// 
//       tauTensor(XX,XX) = m_muWall*(2.0*gradU[XX] - divU);
//       tauTensor(YY,YY) = m_muWall*(2.0*gradV[YY] - divU);
// 
//       tauTensor(XX,YY) = tauTensor(YY,XX) = m_muWall*(gradU[YY] + gradV[XX]);
// 
//       if(adim)
//       {
//         const CFreal machInf = m_updateVarSet->getModel()->getMachInf();
//         const CFreal Re = m_diffVar->getModel().getReynolds();
// 
//         /// @todo check this
//         skinFrictionTensor = tauTensor/(0.5*machInf*sqrt(gamma)/Re);
//       }
//       else
//       {
//         skinFrictionTensor = tauTensor/(0.5*m_rhoInf*m_uInf*m_uInf);
//       }
// 
//       // compute the viscous force on this wall face.
//       frictionForce[XX] = refTau*(tauTensor(XX,XX)*normal[XX] + tauTensor(XX,YY)*normal[YY]);
//       frictionForce[YY] = refTau*(tauTensor(YY,XX)*normal[XX] + tauTensor(YY,YY)*normal[YY]);
// 
//       // compute adimensional skin friction
//       RealVector skinFrictionVector(2);
//       skinFrictionVector[XX] = skinFrictionTensor(XX,XX)*normal[XX] + skinFrictionTensor(XX,YY)*normal[YY];
//       skinFrictionVector[YY] = skinFrictionTensor(YY,XX)*normal[XX] + skinFrictionTensor(YY,YY)*normal[YY];
//       if (m_flowDir[XX]*normal[YY] - m_flowDir[YY]*normal[XX] > 0.0)
//       {
//         // sign is reversed to obtain the force on the body
//         skinFriction = - skinFrictionVector[XX]*normal[YY] + skinFrictionVector[YY]*normal[XX];
//       }
//       else
//       {
//         // sign is reversed to obtain the force on the body
//         skinFriction = skinFrictionVector[XX]*normal[YY] - skinFrictionVector[YY]*normal[XX];
//       }
//     }
//     else
//     {
//       cf_assert(m_dim == DIM_3D);
// 
//       const RealVector& gradU = *(m_cellGradFlxPnt[iFlx][1]);
//       const RealVector& gradV = *(m_cellGradFlxPnt[iFlx][2]);
//       const RealVector& gradW = *(m_cellGradFlxPnt[iFlx][3]);
//       const CFreal divU = (2.0/3.0)*(gradU[XX] + gradV[YY] + gradW[ZZ]);
// 
//       tauTensor(XX,XX) = m_muWall*(2.0*gradU[XX] - divU);
//       tauTensor(YY,YY) = m_muWall*(2.0*gradV[YY] - divU);
//       tauTensor(ZZ,ZZ) = m_muWall*(2.0*gradW[ZZ] - divU);
// 
//       tauTensor(XX,YY) = tauTensor(YY,XX) = m_muWall*(gradU[YY] + gradV[XX]);
//       tauTensor(XX,ZZ) = tauTensor(ZZ,XX) = m_muWall*(gradU[ZZ] + gradW[XX]);
//       tauTensor(YY,ZZ) = tauTensor(ZZ,YY) = m_muWall*(gradV[ZZ] + gradW[YY]);
// 
//       if(adim)
//       {
//         const CFreal machInf = m_updateVarSet->getModel()->getMachInf();
//         const CFreal Re = m_diffVar->getModel().getReynolds();
// 
//         skinFrictionTensor = tauTensor/(0.5*machInf*sqrt(gamma)/Re);
//       }
//       else
//       {
//         skinFrictionTensor = tauTensor/(0.5*m_rhoInf*m_uInf*m_uInf);
//       }
// 
//       //Compute the viscous force on this wall face.
//       frictionForce[XX] = refTau*(tauTensor(XX,XX)*normal[XX] + tauTensor(XX,YY)*normal[YY] + tauTensor(XX,ZZ)*normal[ZZ]);
//       frictionForce[YY] = refTau*(tauTensor(YY,XX)*normal[XX] + tauTensor(YY,YY)*normal[YY] + tauTensor(YY,ZZ)*normal[ZZ]);
//       frictionForce[ZZ] = refTau*(tauTensor(ZZ,XX)*normal[XX] + tauTensor(ZZ,YY)*normal[YY] + tauTensor(ZZ,ZZ)*normal[ZZ]);
// 
//       // compute adimensional skin friction
//       RealVector skinFrictionVector(3);
//       skinFrictionVector[XX] = skinFrictionTensor(XX,XX)*normal[XX] + skinFrictionTensor(XX,YY)*normal[YY] + skinFrictionTensor(XX,ZZ)*normal[ZZ];
//       skinFrictionVector[YY] = skinFrictionTensor(YY,XX)*normal[XX] + skinFrictionTensor(YY,YY)*normal[YY] + skinFrictionTensor(YY,ZZ)*normal[ZZ];
//       skinFrictionVector[ZZ] = skinFrictionTensor(YY,XX)*normal[XX] + skinFrictionTensor(YY,YY)*normal[YY] + skinFrictionTensor(ZZ,ZZ)*normal[ZZ];
//       skinFriction = skinFrictionVector.norm2();// like this, the sign of this coefficient is always positive. How to keep the sign?
//     }
// 
//     // heat flux
//     const CFreal heatFluxDim = m_diffVar->getHeatFlux(state,m_cellGradFlxPnt[iFlx],normal); // check dimensionality
// 
//     // stanton number
//     const CFreal stanton = 0.0;
//    const CFreal stanton = heatFluxDim/(m_rhoInf*pow(m_uInf,3.0));

    // Compute the friction at the wall
    computeTauWall(iFlx);

    updateWriteData(iFlx);

    // Do some extra computation on the face if needed
    computeExtraValues();

//    // Compute y+ value
//    if(socket_wallDistance.isConnected())  {
//      computeYplus();
//    }
  }
}
    
//////////////////////////////////////////////////////////////////////////////
    
void NavierStokesSkinFrictionHeatFluxFR::computeTauWall(CFuint flxIdx)
{
  // const CFreal Uref = m_updateVarSet->getModel()->getVelRef();
  //const CFreal refLength = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  //const CFreal mu_ref = m_diffVar->getModel().getDynViscosityDim
  //(m_updateVarSet->getModel()->getPressRef(), m_updateVarSet->getModel()->getTempRef());
  
  m_muWall = m_diffVar->getDynViscosity(*(m_cellStatesFlxPnt[flxIdx]), m_cellGradFlxPnt[flxIdx]);
  
  // this will not work adimensional
  // const CFreal tauRef = 1.0; // mu_ref*Uref/refLength;
  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
  const bool adim = Framework::PhysicalModelStack::getActive()->getImplementor()->isAdimensional();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  m_tau = 0.0;
  m_Cf = 0.0;
  
  if (dim == DIM_2D) {
    
    m_tau = m_muWall*
      (m_unitNormalFlxPnts[flxIdx][YY]*MathFunctions::innerProd(*(m_cellGradFlxPnt[flxIdx][m_UID]),m_unitNormalFlxPnts[flxIdx]) -
       m_unitNormalFlxPnts[flxIdx][XX]*MathFunctions::innerProd(*(m_cellGradFlxPnt[flxIdx][m_VID]),m_unitNormalFlxPnts[flxIdx]));
    
    if(adim){
      //check this
      const CFreal mInf = m_updateVarSet->getModel()->getMachInf();
      const CFreal Re = m_diffVar->getModel().getReynolds();
      m_Cf = m_tau / (0.5*mInf*sqrt(gamma)*Re);
    }
    else{
      // friction coefficient 
      m_Cf = m_tau / (0.5*m_rhoInf*m_uInf*m_uInf);
    } 
    
    if (m_cellStatesFlxPnt[flxIdx]->size() == 8)
    {
      m_yPlus = (*(m_cellStatesFlxPnt[flxIdx]))[7]*sqrt(m_muWall*m_rhoWall)/(0.5*m_rhoInf*m_uInf*m_uInf);
    }
    
    m_frictionForces[XX] =  m_Cf*m_unitNormalFlxPnts[flxIdx][YY];
    m_frictionForces[YY] = -m_Cf*m_unitNormalFlxPnts[flxIdx][XX];
  }
  else{
    cf_assert(dim == DIM_3D);
    
    CFreal divU = (2.0/3.0)*((*(m_cellGradFlxPnt[flxIdx][m_UID]))[0] + (*(m_cellGradFlxPnt[flxIdx][m_VID]))[1] + (*(m_cellGradFlxPnt[flxIdx][m_WID]))[2]) ;
    
    m_tau3D(XX,XX) = m_muWall*(2.0*(*(m_cellGradFlxPnt[flxIdx][m_UID]))[XX] - divU) ;
    m_tau3D(YY,YY) = m_muWall*(2.0*(*(m_cellGradFlxPnt[flxIdx][m_VID]))[YY] - divU) ;
    m_tau3D(ZZ,ZZ) = m_muWall*(2.0*(*(m_cellGradFlxPnt[flxIdx][m_WID]))[ZZ] - divU) ;
    
    m_tau3D(XX,YY) = m_tau3D(YY,XX) = m_muWall*((*(m_cellGradFlxPnt[flxIdx][m_UID]))[YY] + (*(m_cellGradFlxPnt[flxIdx][m_VID]))[XX]) ;
    m_tau3D(XX,ZZ) = m_tau3D(ZZ,XX) = m_muWall*((*(m_cellGradFlxPnt[flxIdx][m_UID]))[ZZ] + (*(m_cellGradFlxPnt[flxIdx][m_WID]))[XX]) ;
    m_tau3D(YY,ZZ) = m_tau3D(ZZ,YY) = m_muWall*((*(m_cellGradFlxPnt[flxIdx][m_VID]))[ZZ] + (*(m_cellGradFlxPnt[flxIdx][m_WID]))[YY]) ;
      
    if(adim){
      // check this
      const CFreal mInf = m_updateVarSet->getModel()->getMachInf();
      const CFreal Re = m_diffVar->getModel().getReynolds();
      m_Cf3D = m_tau3D / (0.5*mInf*sqrt(gamma)*Re);
    }
    else{
      CFreal rhoInf = m_pInf / (m_updateVarSet->getModel()->getR() * m_TInf);
      m_Cf3D = m_tau3D / (0.5*rhoInf*m_uInf*m_uInf);
    }
    
    //Compute the viscous force coefficients on this wall face.
    m_frictionForces = m_Cf3D * m_unitNormalFlxPnts[flxIdx];
  }
  
  // forces coefficients must be still divided by the wet surface
  m_frictionForces /= m_refArea;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxFR::prepareOutputFileWall()
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

void NavierStokesSkinFrictionHeatFluxFR::updateOutputFileWall()
{  
  const std::string nsp = getMethodData().getNamespace();
  PE::GetPE().setBarrier(nsp);
  
  // all processors will write their own data one after the other 
  for (CFuint i = 0; i < PE::GetPE().GetProcessorCount(nsp); ++i) 
  {
    if (i == PE::GetPE().GetRank (nsp)) 
    {
      if (getCurrentTRS()->getLocalNbGeoEnts() > 0) 
      {
	SafePtr<TopologicalRegionSet> currTrs = getCurrentTRS();
	boost::filesystem::path file = Environment::DirPaths::getInstance().getResultsDir() /
	  boost::filesystem::path(m_nameOutputFileWall + currTrs->getName());
	file = Framework::PathAppender::getInstance().appendAllInfo  
	  (file,this->m_appendIter,this->m_appendTime,false);   
     
	SelfRegistPtr<Environment::FileHandlerOutput> fhandle =
	  Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
	ofstream& fout = fhandle->open(file, ios::app);
	
	Common::SafePtr<GeometricEntityPool<FaceToCellGEBuilder> >
	  geoBuilder = m_faceBuilder;
	
	SafePtr<FaceToCellGEBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
	//geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);
	
	// get InnerCells TopologicalRegionSet
        SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");
	
	FaceToCellGEBuilder::GeoData& geoData = geoBuilder->getDataGE();
	geoData.cellsTRS = cellTrs;
        geoData.facesTRS = currTrs;
        geoData.isBoundary = true;
	
	const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
	for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) 
	{
	  // build the GeometricEntity
	  geoData.idx = iFace;
	  m_currFace = geoBuilder->buildGE();
	  //m_fvmccData->getCurrentFace() = m_currFace;
	  
	  // GET THE NEIGHBOURING CELL
          m_intCell = m_currFace->getNeighborGeo(0);

          // GET THE STATES IN THE NEIGHBOURING CELL
          m_cellStates = m_intCell->getStates();
	  
	  // only faces whose internal State is parallel updatable will write
	  // their data to avoid redudance due to overlap 
	  if ((*m_cellStates)[0]->isParUpdatable()) 
	  {
	    // loop over flx pnts
            for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
            {
              // compute coordinates of output point
              m_coord = m_currFace->computeCoordFromMappedCoord((*m_flxLocalCoords)[iFlx]);
	      
	      fout << m_coord << " ";
	      
	      const CFuint index = m_mapTrsFaceToID.find(m_currFace->getID()*m_nbrFaceFlxPnts+iFlx);
	      const CFuint nbVars = m_valuesMat.nbRows(); 
	      for (CFuint iVar = 0; iVar < nbVars; ++iVar) 
	      {
	        fout << m_valuesMat(iVar, index) << " ";
	      }
	      fout << "\n";
            }
            
// 	    // compute the face normal
// 	    const vector<Node*>& faceNodes = *m_currFace->getNodes();
// 	    const CFuint nbFaceNodes = faceNodes.size();
// 	    
// 	    // compute the face mid point
// 	    m_coord = 0.0;
// 	    for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode) 
// 	    {
// 	      m_coord += *faceNodes[iNode];
// 	    }
// 	    m_coord /= nbFaceNodes;
// 	    
	  }  
	  geoBuilder->releaseGE();
	}
      }
    }
    
    PE::GetPE().setBarrier(nsp);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxFR::updateWriteData(CFuint flxIdx)
{  
  const CFreal refLength = PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  const CFreal heatFluxRef = m_updateVarSet->getModel()->getTempRef()/refLength;
  const CFuint index = m_mapTrsFaceToID.find(m_currFace->getID()*m_nbrFaceFlxPnts+flxIdx);
  
  // compute the radiative heat flux
  m_heatFluxRad = 0.;
  if (m_hasRadiationCoupling) {
    // AL: the ordering of the wall TRS's must be consistent with the one in RadiativeTransferMonteCarlo
    cf_assert(m_qradFluxWall.size() > 0);
    cf_assert(index < m_qradFluxWall.size());
    m_heatFluxRad = m_qradFluxWall[index]*heatFluxRef;
  } 
  
  const CFreal heatFlux = m_diffVar->getHeatFlux( *(m_cellStatesFlxPnt[flxIdx]), m_cellGradFlxPnt[flxIdx], m_unitNormalFlxPnts[flxIdx])*heatFluxRef + m_heatFluxRad;
  //CFLog(INFO, "state: " << *(m_cellStatesFlxPnt[flxIdx]) << ", Tgrad: " << (*(m_cellGradFlxPnt[flxIdx][3])) << ", ref: " << heatFluxRef << "\n");
  
  if (m_cellStatesFlxPnt[flxIdx]->size() == 8)
  {
    m_heatFluxRad = (*(m_cellStatesFlxPnt[flxIdx]))[6];
  }
  
  CFreal pDim = 0.;
  CFreal rhoDim = 0.;
  CFreal TDim = 0.;

  computeDimensionalPressDensTemp(pDim, rhoDim, TDim, flxIdx);
  
  const CFreal rhoInf = m_pInf / (m_updateVarSet->getModel()->getRdim() * m_TInf);
  //const CFreal rhoInf = 1.0;
  
  CFreal stantonNumber = 0.0;
  switch(m_stantonNumID) {
  case(0):
    stantonNumber = heatFlux/(rhoInf*pow(m_uInf,3.0));
    break;
  case(1):
    stantonNumber = heatFlux / ((m_updateVarSet->getModel()->getCp()*(m_TInf - TDim) + 
				 0.5* m_uInf*m_uInf )*rhoInf*m_uInf);
    break;
  case(2):
    stantonNumber = heatFlux / (m_updateVarSet->getModel()->getCp()*m_rhoWall*m_uInf);
    break;
  default:
    //stantonNumber = heatFlux/(rhoInf*pow(m_uInf,3.0));
    stantonNumber = heatFlux/(m_updateVarSet->getModel()->getCp()*(m_TInf - TDim)*rhoInf*m_uInf);
    break;
  }
  
  CFreal Cp = (pDim - m_pInf);

  Cp /= (0.5*rhoInf*m_uInf*m_uInf);
  
  // fill in the 2D-array with all the data to be output 
  if (PhysicalModelStack::getActive()->getDim() == DIM_2D){
    updateValuesMatAndResidual(0, index, pDim);    
    updateValuesMatAndResidual(1, index, TDim);    
    updateValuesMatAndResidual(2, index, rhoDim);    
    updateValuesMatAndResidual(3, index, Cp);    
    updateValuesMatAndResidual(4, index, heatFlux);    
    updateValuesMatAndResidual(5, index, stantonNumber);
    updateValuesMatAndResidual(6, index, this->m_yPlus);
    updateValuesMatAndResidual(7, index, this->m_Cf);
    updateValuesMatAndResidual(8, index, m_muWall);
    updateValuesMatAndResidual(9, index, m_heatFluxRad);
  }
  else{
    throw Common::NotImplementedException (FromHere(),"NavierStokesSkinFrictionHeatFluxFR::updateOutputFile() is only 2D");
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxFR::computeYplus()
{
  CFAUTOTRACE;

  m_yPlus = 0.;

//   DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();
// 
//   // Compute the distance to the first cell: y0
//   //First get the inner state
//   State* innerState = m_currFace->getState(0);
// 
//   CFreal y0 = wallDistance[innerState->getLocalID()];
// 
//   // Compute the y+
//   // y^{+} = \frac {\sqrt{\rho} * \sqrt{\tau} * y0} {\mu}
// 
//   const CFreal refSpeed = (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::V];
//   const CFreal rhoRef = (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::RHO];
// 
//   m_yPlus = sqrt(m_rhoWall) * sqrt(fabs(m_tau)) * y0 / m_muWall;
//   m_yPlus *= sqrt(rhoRef) * refSpeed;

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxFR::computeDimensionalPressDensTemp(CFreal& pDim, CFreal& rhoDim, CFreal& TDim, CFuint flxIdx)
{
  CFAUTOTRACE;
  
  const bool Puvt = m_frData->getUpdateVarStr() == "Puvt";
  
  const CFreal rhoRef = (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::RHO];
  
  // output the data
  const CFuint PID = 0;
  
  if (Puvt)
  {
    pDim = m_updateVarSet->getModel()->getPressureFromState((*(m_cellStatesFlxPnt[flxIdx]))[PID]) * 
      (m_updateVarSet->getModel()->getPressRef());
    TDim = (*(m_cellStatesFlxPnt[flxIdx]))[m_TID] * (m_updateVarSet->getModel()->getTempRef());
    rhoDim = m_rhoWall * rhoRef;
  }
  else
  {
    // get some data needed further
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    const CFuint nbEqsM1 = nbEqs - 1;
    const CFreal R = m_updateVarSet->getModel()->getR();
    const CFreal gamma = m_updateVarSet->getModel()->getGamma();
    const CFreal gammaMinus1 = gamma - 1.;
  
    const CFreal rho = (*(m_cellStatesFlxPnt[flxIdx]))[0];
    const CFreal invRho = 1./rho;
    CFreal rhoK2 = 0.0;
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      rhoK2 += (*(m_cellStatesFlxPnt[flxIdx]))[iDim+1]*(*(m_cellStatesFlxPnt[flxIdx]))[iDim+1];
    }
    rhoK2 *= 0.5*invRho;
    const CFreal p = gammaMinus1*((*(m_cellStatesFlxPnt[flxIdx]))[nbEqsM1] - rhoK2);

    // temperature
    const CFreal T = p*invRho/R;

    // dimensional values
    TDim   = T   * m_updateVarSet->getModel()->getTempRef();
    pDim   = p   * m_updateVarSet->getModel()->getPressRef();
    rhoDim = rho * (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::RHO];
  }  
}
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////




