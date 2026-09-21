#include <iomanip>

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
  m_varNames.push_back(hasTransitionLayout() ? "gamma" : "heatFRadiative");
  
  cf_always_assert(this->m_varNames.size() == 10 + dim); 

  // velocity components in the states
  m_velocityIDs.resize(m_dim);
  m_velocityIDs[XX] = m_UID;
  m_velocityIDs[YY] = m_VID;
  if (m_dim == DIM_3D)
  {
    m_velocityIDs[ZZ] = m_WID;
  }

  m_traction.resize(m_dim);
  
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
  for (CFuint iFlx = 0; iFlx < m_nbrFaceFlxPnts; ++iFlx)
  {
    const RealVector& state = *m_cellStatesFlxPnt[iFlx];

    // this is needed for LTE
    m_diffVar->setComposition(state,false,0);

    m_rhoWall = m_diffVar->getDensity(state);

    // friction at the wall
    computeTauWall(iFlx);

    updateWriteData(iFlx);

    // extra computation on the face if needed
    computeExtraValues();
  }
}

//////////////////////////////////////////////////////////////////////////////

bool NavierStokesSkinFrictionHeatFluxFR::hasTransitionLayout() const
{
  const std::string updateVarStr = m_frData->getUpdateVarStr();
  const std::string convectiveName = PhysicalModelStack::getActive()->getConvectiveName();

  return updateVarStr == "Puvt" && (convectiveName.find("GReKLogO") != std::string::npos || convectiveName.find("GReKO") != std::string::npos);
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokesSkinFrictionHeatFluxFR::computeStantonNumber(CFreal heatFlux, CFreal temperature, CFuint flxIdx)
{
  switch (m_stantonNumID)
  {
    case 0:
      return heatFlux/(m_rhoInf*std::pow(m_uInf,3.));
    case 1:
      return heatFlux/((m_updateVarSet->getModel()->getCp()*(m_TInf-temperature)+0.5*m_uInf*m_uInf)*m_rhoInf*m_uInf);
    case 2:
      return heatFlux/(m_updateVarSet->getModel()->getCp()*m_rhoWall*m_uInf);
    default:
      return heatFlux/(m_updateVarSet->getModel()->getCp()*(m_TInf-temperature)*m_rhoInf*m_uInf);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxFR::computeTauWall(CFuint flxIdx)
{
  m_muWall = m_diffVar->getDynViscosity(*m_cellStatesFlxPnt[flxIdx],m_cellGradFlxPnt[flxIdx]);

  const RealVector& normal = m_bndFaceDiffData.unitNormals[flxIdx];
  const RealVector& diffFlux = m_bndFaceDiffData.diffFluxes[flxIdx];

  const bool adim = PhysicalModelStack::getActive()->getImplementor()->isAdimensional();
  const CFreal scale = adim ? m_updateVarSet->getModel()->getPressRef() :
    1./PhysicalModelStack::getActive()->getImplementor()->getRefLength();
  const CFreal dynamicPressure = 0.5*m_rhoInf*m_uInf*m_uInf;

  // tangential viscous traction on the fluid, tau_t = (I - n n^T) F_mom
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    m_traction[iDim] = scale*diffFlux[m_velocityIDs[iDim]];
  }
  m_traction -= MathFunctions::innerProd(m_traction,normal)*normal;

  // friction force coefficients on the body, Cf = -tau_t/(q_inf refArea)
  m_frictionForces = -m_traction/(dynamicPressure*m_refArea);
  m_frictionForcesFlxPnts[flxIdx] = m_frictionForces;

  // 2D: component of the force on the body along (-ny, nx); 3D: magnitude
  m_tau = m_dim == DIM_2D ? m_traction[XX]*normal[YY]-m_traction[YY]*normal[XX] : m_traction.norm2();
  m_Cf = m_tau/dynamicPressure;

  m_yPlus = 0.;
  if (hasTransitionLayout())
  {
    m_yPlus = (*m_cellStatesFlxPnt[flxIdx])[m_dim+5]*sqrt(m_muWall*m_rhoWall)/dynamicPressure;
  }
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
    fout << "# n points out of fluid; heatF positive into wall; Cf/forces act on body (2D Cf tangent=(-ny,nx)). Physical diffusive flux only (no convective or LLAV flux).\n";
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
        fout << std::setprecision(17);
	
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
            for (CFuint iFlx = 0; iFlx < getNbrFaceFlxPnts(m_currFace->getID()); ++iFlx)
            {
              // compute coordinates of output point, with the flux points of the face type in 3D
              Common::SafePtr< std::vector< RealVector > > flxLocalCoords = m_flxLocalCoords;
              if (m_dim == 3)
              {
                const CFuint faceType = m_currFace->getShape() == CFGeoShape::TRIAG ? 0 : 1;
                flxLocalCoords = &(*m_frData->getFRLocalData()[0]->getFaceFlxPntsLocalCoordsPerType())[faceType];
              }
              m_coord = m_currFace->computeCoordFromMappedCoord((*flxLocalCoords)[iFlx]);
	      
	      fout << m_coord << " ";
	      
	      const CFuint index = m_mapTrsFaceToID.find(m_currFace->getID()*m_nbrFaceFlxPntsMax+iFlx);
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
  const bool adim = PhysicalModelStack::getActive()->getImplementor()->isAdimensional();
  const CFreal velocityRef = m_updateVarSet->getModel()->getVelRef();
  const CFreal heatFluxRef = adim ?
    m_updateVarSet->getModel()->getReferencePhysicalData()[EulerTerm::RHO]*velocityRef*velocityRef*velocityRef : 1./refLength;
  const CFuint index = m_mapTrsFaceToID.find(m_currFace->getID()*m_nbrFaceFlxPntsMax+flxIdx);
  
  // compute the radiative heat flux
  m_heatFluxRad = 0.;
  if (m_hasRadiationCoupling) {
    // AL: the ordering of the wall TRS's must be consistent with the one in RadiativeTransferMonteCarlo
    cf_assert(m_qradFluxWall.size() > 0);
    cf_assert(index < m_qradFluxWall.size());
    m_heatFluxRad = m_qradFluxWall[index]*heatFluxRef;
  } 
  
  // heat flux into the wall, heatF = u_b.F_mom - F_E, with F_E the diffusive
  // energy flux (conduction, species and modal enthalpy, viscous work)
  m_updateVarSet->computePhysicalData(*m_cellStatesFlxPnt[flxIdx],m_dataState);
  const RealVector& diffFlux = m_bndFaceDiffData.diffFluxes[flxIdx];

  CFreal viscousWork = 0.;
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    viscousWork += m_dataState[EulerTerm::VX+iDim]*diffFlux[m_velocityIDs[iDim]];
  }

  const CFreal heatFlux = (viscousWork-diffFlux[m_TID])*heatFluxRef + m_heatFluxRad;

  // last output column: gamma of the transition model or the radiative heat flux
  const CFreal extraOutput = hasTransitionLayout() ? (*m_cellStatesFlxPnt[flxIdx])[m_dim+4] : m_heatFluxRad;

  CFreal pDim = 0.;
  CFreal rhoDim = 0.;
  CFreal TDim = 0.;

  computeDimensionalPressDensTemp(pDim, rhoDim, TDim, flxIdx);
  
  const CFreal stantonNumber = computeStantonNumber(heatFlux,TDim,flxIdx);

  CFreal Cp = (pDim - m_pInf);

  Cp /= (0.5*m_rhoInf*m_uInf*m_uInf);
  
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
    updateValuesMatAndResidual(9, index, extraOutput);
  }
  else{
    updateValuesMatAndResidual(0, index, pDim);    
    updateValuesMatAndResidual(1, index, TDim);    
    updateValuesMatAndResidual(2, index, rhoDim);    
    updateValuesMatAndResidual(3, index, Cp);    
    updateValuesMatAndResidual(4, index, heatFlux);    
    updateValuesMatAndResidual(5, index, stantonNumber);
    updateValuesMatAndResidual(6, index, this->m_yPlus);
    updateValuesMatAndResidual(7, index, this->m_Cf);
    updateValuesMatAndResidual(8, index, m_muWall);
    updateValuesMatAndResidual(9, index, extraOutput);
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
  
  const bool Puvt = m_frData->getUpdateVarStr() == "Puvt" || m_frData->getUpdateVarStr() == "Pvt";
  
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




