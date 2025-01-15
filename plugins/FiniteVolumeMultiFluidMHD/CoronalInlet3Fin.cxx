#include "Framework/MethodCommandProvider.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "FiniteVolumeMultiFluidMHD/CoronalInlet3Fin.hh"
#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CoronalInlet3Fin, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> 
coronalInlet3FinFVMCCProvider("CoronalInlet3FinFVMCC");

//////////////////////////////////////////////////////////////////////////////

void CoronalInlet3Fin::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFuint> >("ProjectionIDs","IDs corresponding to projection fields.");
  options.addConfigOption< vector<CFuint> >("VarIDs","IDs of the variables from which values are read by file");
  options.addConfigOption< bool >("m_BfromFile","");
  options.addConfigOption< bool >("rotate","");
  options.addConfigOption< CFreal >("Vin","Boundary velocity of ions and neutrals");
  options.addConfigOption< CFreal >("rhoin","Boundary density of ions");
  options.addConfigOption< CFreal >("Tin","Boundary temperature (Ti = Tn)");
  options.addConfigOption< CFreal >("Cin","Assumed boundary concentration of ions");
  options.addConfigOption< CFreal >("MptoMe","Assumed mass ratio between protons and electrons");
}
      
//////////////////////////////////////////////////////////////////////////////

CoronalInlet3Fin::CoronalInlet3Fin(const std::string& name) :
  SuperInlet(name),
  m_mapGeoToTrs(),
  m_mapTrs2Twall(),
  m_BfromFile(false)
{
  addConfigOptionsTo(this);

  
  _projectionIDs = vector<CFuint>();
  setParameter("ProjectionIDs",&_projectionIDs);
  m_varIDs = vector<CFuint>();
  setParameter("VarIDs",&m_varIDs);

  _m_BfromFile = false;
  _rotate = false;
  _Vin = 1935.7;
  _rhoin = 1.67e-13;
  _Tin = 1.8e6; 
  _Cin = 1e-6;
  _mptome = 1836.15267343; 

  setParameter("m_BfromFile",&_m_BfromFile);
  setParameter("rotate",&_rotate);
  setParameter("Vin",&_Vin);
  setParameter("rhoin",&_rhoin);
  setParameter("Tin",&_Tin);
  setParameter("Cin",&_Cin);
  setParameter("MptoMe",&_mptome);

}

//////////////////////////////////////////////////////////////////////////////

CoronalInlet3Fin::~CoronalInlet3Fin()
{
}

//////////////////////////////////////////////////////////////////////////////

void CoronalInlet3Fin::setup()
{
  SuperInlet::setup();

  if (m_varIDs.size() > 0) {m_BfromFile = true;}
  
  m_mapGeoToTrs =
    MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");
  
  vector< SafePtr<TopologicalRegionSet> >& trsList = this->getTrsList();
  
  for (CFuint iTrs = 0; iTrs < trsList.size(); ++iTrs) {
    SafePtr<TopologicalRegionSet> trs = trsList[iTrs];
    const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();
    RealVector* tWall = new RealVector(0.0,nbTrsFaces);
    m_mapTrs2Twall.insert(&*trs, tWall);
  }
  m_mapTrs2Twall.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

void CoronalInlet3Fin::unsetup()
{
  SuperInlet::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void CoronalInlet3Fin::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);


  RealVector B_PFSS_dimless(3);
  RealVector B_PFSS(3);

  if (m_initialSolutionMap.size() > 0) {
    SafePtr<MapGeoToTrsAndIdx> mapGeoToTrs =
      MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");
    const CFuint faceIdx = mapGeoToTrs->getIdxInTrs(face->getID());
    const string name = getCurrentTRS()->getName();
    SafePtr<RealVector> initialValues = m_initialSolutionMap.find(name);
    const CFuint nbVars = m_initialSolutionIDs.size();
    const CFuint startID = faceIdx*nbVars;
    for (CFuint i = 0; i < nbVars; ++i) {
      const CFuint varID = m_initialSolutionIDs[i];
      const CFuint idx = startID+i;
      cf_assert(idx < initialValues->size());

      B_PFSS_dimless[i] = (*initialValues)[idx]/2.2e-4; // save the Poisson PFSS solution in a vector
      B_PFSS[i] = (*initialValues)[idx];
      (*ghostState)[varID] = 2.*(*initialValues)[idx] - (*innerState)[varID];
      CFLog(DEBUG_MIN, "CoronalInlet3Fin::setGhostState() => [" << varID << "] => " << (*initialValues)[idx] << " | " << (*innerState)[varID] << "\n");
    }
  }

  CFreal BrFromFile = 0.;
  m_BfromFile = _m_BfromFile;
  if (m_BfromFile) { 
    vector<Node*>& nodesInFace = *face->getNodes();
    const CFuint nbNodesInFace = nodesInFace.size();
    SafePtr<TopologicalRegionSet> trs = m_mapGeoToTrs->getTrs(face->getID());
    cf_assert(trs.isNotNull());
    
    // build the mapTrs2NodalValues storage
    SafePtr<NodalStatesExtrapolator<CellCenterFVMData> > nse =
      this->getMethodData().getNodalStatesExtrapolator();
    SafePtr<vector<NodalStatesExtrapolator<CellCenterFVMData>::MapTrs2NodalValues*> >
      mapTrs2NodalValues = nse->getMapTrs2NodalValues();
    
    cf_assert(m_varIDs.size() == 1);
    cf_assert(m_varIDs[0] == 0);
    
    const CFuint varID = 0;
    cf_assert(varID < mapTrs2NodalValues->size());
    
    RealVector& bArray = *(*mapTrs2NodalValues)[varID]->find(&*trs);
    CFMap<CFuint,CFuint>* mapNodeIDs = nse->getMapTrs2NodeIDs()->find(&*trs);
    
    BrFromFile = 0.0;
    for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
      const CFuint localNodeID = nodesInFace[iNode]->getLocalID();
      // CFLog(INFO, "Barray = " << bArray[mapNodeIDs->find(localNodeID)] << "\n");
      BrFromFile += bArray[mapNodeIDs->find(localNodeID)];
    }
    BrFromFile /= nbNodesInFace;

  }

  CFreal RSun = 6.9551e8; 
  CFreal density_code = _rhoin; // the default value is 1.67e-13;
  CFreal velocity_code = 2.2e-4/std::sqrt(1.256e-6*1.67e-13);
  CFreal pressure_code = std::pow(2.2e-4,2)/1.256e-6;
  CFreal B_code = 2.2e-4;

  CFreal latI = 0.;
  const CFreal PI = MathTools::MathConsts::CFrealPi();
  const CFreal xG = ghostState->getCoordinates()[XX];
  const CFreal xG_dimless = ghostState->getCoordinates()[XX]/RSun;
  const CFreal yG = ghostState->getCoordinates()[YY];
  const CFreal yG_dimless = ghostState->getCoordinates()[YY]/RSun;
  const CFreal zG = ghostState->getCoordinates()[ZZ];
  const CFreal zG_dimless = ghostState->getCoordinates()[ZZ]/RSun;
  const CFreal xI = innerState->getCoordinates()[XX];
  const CFreal xI_dimless = innerState->getCoordinates()[XX]/RSun;
  const CFreal yI = innerState->getCoordinates()[YY];
  const CFreal yI_dimless = innerState->getCoordinates()[YY]/RSun;
  const CFreal zI = innerState->getCoordinates()[ZZ];
  const CFreal zI_dimless = innerState->getCoordinates()[ZZ]/RSun;
  const CFreal xBoundary = (xG + xI)/2.0;
  const CFreal xBoundary_dimless = (xG_dimless + xI_dimless)/2.0;
  const CFreal yBoundary = (yG + yI)/2.0;
  const CFreal yBoundary_dimless = (yG_dimless + yI_dimless)/2.0;
  const CFreal zBoundary = (zG + zI)/2.0;
  const CFreal zBoundary_dimless = (zG_dimless + zI_dimless)/2.0;
  const CFreal rBoundary = std::sqrt(xBoundary*xBoundary + yBoundary*yBoundary + zBoundary*zBoundary);
  const CFreal rBoundary_dimless = std::sqrt(xBoundary_dimless*xBoundary_dimless + yBoundary_dimless*yBoundary_dimless + zBoundary_dimless*zBoundary_dimless);
  const CFreal rhoBoundary = std::sqrt(xBoundary*xBoundary + yBoundary*yBoundary);
  const CFreal rhoBoundary_dimless = std::sqrt(xBoundary_dimless*xBoundary_dimless + yBoundary_dimless*yBoundary_dimless);
  const CFreal thetaBoundary = std::atan2(std::sqrt(xBoundary*xBoundary + yBoundary*yBoundary),zBoundary);
  const CFreal thetaBoundary_dimless = std::atan2(std::sqrt(xBoundary_dimless*xBoundary_dimless + yBoundary_dimless*yBoundary_dimless),zBoundary_dimless);
  const CFreal phiBoundary = std::atan2(yBoundary,xBoundary);
  const CFreal phiBoundary_dimless = std::atan2(yBoundary_dimless,xBoundary_dimless);


  const CFreal rG = std::sqrt(xG*xG + yG*yG + zG*zG);
  const CFreal rG_dimless = std::sqrt(xG_dimless*xG_dimless + yG_dimless*yG_dimless + zG_dimless*zG_dimless);
  const CFreal rhoG = std::sqrt(xG*xG + yG*yG);
  const CFreal rhoG_dimless = std::sqrt(xG_dimless*xG_dimless + yG_dimless*yG_dimless);
  const CFreal thetaG = std::atan2(std::sqrt(xG_dimless*xG_dimless + yG_dimless*yG_dimless),zG_dimless);
  const CFreal phiG = std::atan2(yG_dimless,xG_dimless);
  const CFreal rI = std::sqrt(xI*xI + yI*yI + zI*zI);
  const CFreal rI_dimless = std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless + zI_dimless*zI_dimless);
  const CFreal thetaI = std::atan2(std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless),zI_dimless);
  const CFreal rhoI = std::sqrt(xI*xI + yI*yI);
  const CFreal rhoI_dimless = std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless);
  const CFreal phiI = std::atan2(yI_dimless,xI_dimless);
  const CFreal rB = (rG + rI) / 2.0;



  // DENSITY

  // The mass of ions is set, the mass of electrons is determined from those of protons an the mass ratio
  CFreal m_proton = 1.6726219e-27; 
  CFreal m_electron = m_proton / _mptome; 

  // The inner boundary density of neutrals is determined from the assumed concentration set by the user
  CFreal concentration = _Cin;

  // As the inner boundary, we also assume charge neutrality
  (*ghostState)[8] = 2.0*density_code / _mptome - (*innerState)[8];
  (*ghostState)[9] = 2.0*density_code - (*innerState)[9];
  (*ghostState)[10] = 2.0*density_code * concentration - (*innerState)[10];


  // MAGNETIC FIELD 
  // The principle here is similar as what is applied in ideal-MHD
  CFreal BrBoundary_dimless; 

  // Get only the Br from the PFSS
  BrBoundary_dimless = xBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[0] + yBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[1] + zBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[2];

  // Get the Cartesian dimensionless states
  CFreal BxI_dimless = (*innerState)[0]/B_code;
  CFreal ByI_dimless = (*innerState)[1]/B_code;
  CFreal BzI_dimless = (*innerState)[2]/B_code;

  // Transform them to spherical
  CFreal BrI_dimless = xI_dimless/rI_dimless*BxI_dimless + yI_dimless/rI_dimless*ByI_dimless + zI_dimless/rI_dimless*BzI_dimless;
  CFreal BthetaI_dimless = xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*BxI_dimless + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*ByI_dimless - rhoI_dimless/rI_dimless*BzI_dimless;
  CFreal BphiI_dimless = -yI_dimless/rhoI_dimless*BxI_dimless + xI_dimless/rhoI_dimless*ByI_dimless;

  // The Bphi ghost is set to Bphi inner --> zero gradient
  CFreal BphiG_dimless = BphiI_dimless;
  CFreal BphiBoundary_dimless = (BphiG_dimless + BphiI_dimless)/2.0;

  // With Btheta, we do the same
  CFreal BthetaG_dimless = BthetaI_dimless;
  CFreal BthetaBoundary_dimless = (BthetaI_dimless + BthetaG_dimless)/2.0;

  // Transform back to Cartesian
  CFreal BxBoundary_dimless = xBoundary_dimless/rBoundary_dimless*BrBoundary_dimless - yBoundary_dimless/rhoBoundary_dimless*BphiBoundary_dimless + xBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*BthetaBoundary_dimless;
  CFreal ByBoundary_dimless = yBoundary_dimless/rBoundary_dimless*BrBoundary_dimless + xBoundary_dimless/rhoBoundary_dimless*BphiBoundary_dimless + yBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*BthetaBoundary_dimless;
  CFreal BzBoundary_dimless = zBoundary_dimless/rBoundary_dimless*BrBoundary_dimless - rhoBoundary_dimless/rBoundary_dimless*BthetaBoundary_dimless;

  // Assign
  (*ghostState)[0] = (2*BxBoundary_dimless - (*innerState)[0]/B_code)*B_code;
  (*ghostState)[1] = (2*ByBoundary_dimless - (*innerState)[1]/B_code)*B_code;
  (*ghostState)[2] = (2*BzBoundary_dimless - (*innerState)[2]/B_code)*B_code;

  CFreal BxI = (*innerState)[0];
  CFreal ByI = (*innerState)[1];
  CFreal BzI = (*innerState)[2];

  CFreal BxBoundary= BxBoundary_dimless * B_code; 
  CFreal ByBoundary= ByBoundary_dimless * B_code;
  CFreal BzBoundary= BzBoundary_dimless * B_code;


  

  // VELOCITY
  // The principle here is again similar to ideal-MHD --> we prescribe a velocity and correct it to make it have the same direction as the B field


  /*CFreal VxBoundaryn;
  CFreal VyBoundaryn;
  CFreal VzBoundaryn;

  CFreal VxBoundaryi;
  CFreal VyBoundaryi;
  CFreal VzBoundaryi;

  // We prescribe an initial radial outflow
  CFreal VrBoundary = 1935.0; //m/s
  CFreal VthetaBoundary = 0.0;
  CFreal VphiBoundary = 0.0;

  CFreal VxBoundary= xBoundary/rBoundary*VrBoundary- yBoundary/rhoBoundary*VphiBoundary+ xBoundary*zBoundary/(rhoBoundary*rBoundary)*VthetaBoundary;
  CFreal VyBoundary= yBoundary/rBoundary*VrBoundary+ xBoundary/rhoBoundary*VphiBoundary+ yBoundary*zBoundary/(rhoBoundary*rBoundary)*VthetaBoundary;
  CFreal VzBoundary= zBoundary/rBoundary*VrBoundary- rhoBoundary/rBoundary*VthetaBoundary;

  VxBoundaryi = VxBoundary;
  VyBoundaryi = VyBoundary;
  VzBoundaryi = VzBoundary;
  VxBoundaryn = VxBoundary;
  VyBoundaryn = VyBoundary;
  VzBoundaryn = VzBoundary;

  (*ghostState)[11] = (2.*VxBoundaryi- (*innerState)[11]);
  (*ghostState)[12] = (2.*VyBoundaryi- (*innerState)[12]);
  (*ghostState)[13] = (2.*VzBoundaryi- (*innerState)[13]);
  (*ghostState)[14] = (2.*VxBoundaryi- (*innerState)[14]);
  (*ghostState)[15] = (2.*VyBoundaryi- (*innerState)[15]);
  (*ghostState)[16] = (2.*VzBoundaryi- (*innerState)[16]);
  (*ghostState)[17] = (2.*VxBoundaryn- (*innerState)[17]);
  (*ghostState)[18] = (2.*VyBoundaryn- (*innerState)[18]);
  (*ghostState)[19] = (2.*VzBoundaryn- (*innerState)[19]);*/

  // We compue the unit vector of the magnetic field and from there determine the direction of the velocity field
  CFreal BxB = BxBoundary;
  CFreal ByB = ByBoundary;
  CFreal BzB = BzBoundary;
  CFreal BBmag = std::sqrt(BxB*BxB+ ByB*ByB+ BzB*BzB);
  CFreal VBmag = _Vin; 
  CFreal BxuB = BxB/BBmag;
  CFreal ByuB = ByB/BBmag;
  CFreal BzuB = BzB/BBmag;
  CFreal VxBB = BxuB * VBmag;
  CFreal VyBB = ByuB * VBmag;
  CFreal VzBB = BzuB * VBmag;

  CFreal VrBB = xBoundary/rBoundary*VxBB+ yBoundary/rBoundary*VyBB + zBoundary/rBoundary*VzBB;
  CFreal VthetaBB = xBoundary*zBoundary/(rhoBoundary*rBoundary)*VxBB+ yBoundary*zBoundary/(rhoBoundary*rBoundary)*VyBB- rhoBoundary/rBoundary*VzBB;
  CFreal VphiBB = -yBoundary/rhoBoundary*VxBB + xBoundary/rhoBoundary*VyBB;

  // If the flow points inwards, we flip the vector
  if (VrBB < 0.0){
    VxBB = -VxBB;
    VyBB = -VyBB;
    VzBB = -VzBB;
  }

  // We assign to the species (all the same)
  (*ghostState)[11] = (2.*VxBB- (*innerState)[11]);
  (*ghostState)[12] = (2.*VyBB- (*innerState)[12]);
  (*ghostState)[13] = (2.*VzBB- (*innerState)[13]);
  (*ghostState)[14] = (2.*VxBB- (*innerState)[14]);
  (*ghostState)[15] = (2.*VyBB- (*innerState)[15]);
  (*ghostState)[16] = (2.*VzBB- (*innerState)[16]);
  (*ghostState)[17] = (2.*VxBB- (*innerState)[17]);
  (*ghostState)[18] = (2.*VyBB- (*innerState)[18]);
  (*ghostState)[19] = (2.*VzBB- (*innerState)[19]);


  // TEMPERATURE
  CFreal T_Sun = 1.5e6; 
  CFreal mu_cor = 1.27; 
  CFreal mH = 1.67e-27;   
  CFreal BRef = 2.2e-4;   
  CFreal mu0 = 1.2566e-6;
  CFreal rhoRef = 1.67e-13;
  CFreal kB = 1.38e-23;

  // The temperature is set such that pressure is similar to what we had in ideal MHD
  CFreal pBoundary = 0.25801090625/2.0 * _rhoin/1.67e-13;
  CFreal TBoundary = pBoundary/ (2.0 * density_code / mH) / kB;

  // For now, we assume that all the fluids have the same temperature at the inner boundary
  (*ghostState)[20] = (2.* _Tin - (*innerState)[20]);
  (*ghostState)[21] = (2.* _Tin - (*innerState)[21]); 
  (*ghostState)[22] = (2.* _Tin - (*innerState)[22]); 


  // ELECTRIC FIELD
  // We assume it to be zero at the boundary (perfect conductor) 
  (*ghostState)[3] = -(*innerState)[3]; 
  (*ghostState)[4] = -(*innerState)[4]; 
  (*ghostState)[5] = -(*innerState)[5]; 


  // DIVERGENCE CLEANING
  (*ghostState)[6] = -(*innerState)[6]; // for a perfect conductor, psi(boundary) --> 0.0
  (*ghostState)[7] = (*innerState)[7]; // for a perfect conductor, the gradient of phi(boundary) is 0  --> d phi / d n = 0

}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
