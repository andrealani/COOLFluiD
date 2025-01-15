#include "Framework/MethodCommandProvider.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "FiniteVolumeMultiFluidMHD/CoronalInlet2Fin.hh"
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

MethodCommandProvider<CoronalInlet2Fin, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> 
coronalInlet2FinFVMCCProvider("CoronalInlet2FinFVMCC");

//////////////////////////////////////////////////////////////////////////////

void CoronalInlet2Fin::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFuint> >("ProjectionIDs","IDs corresponding to projection fields.");
  options.addConfigOption< vector<CFuint> >("VarIDs","IDs of the variables from which values are read by file");

  options.addConfigOption< bool >("m_BfromFile","");
  options.addConfigOption< bool >("rotate","");
  options.addConfigOption< bool >("B_theta","");
  options.addConfigOption< CFreal >("Vin","");
  options.addConfigOption< CFreal >("rhoin","");
  options.addConfigOption< CFreal >("Tin","");
  options.addConfigOption< CFreal >("Cin","");

}
      
//////////////////////////////////////////////////////////////////////////////

CoronalInlet2Fin::CoronalInlet2Fin(const std::string& name) :
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
  _B_theta = true;
  _Vin = 1935.7;
  _rhoin = 1.67e-13;
  _Tin = 1.8e6; 
  _Cin = 1e-6;
  setParameter("m_BfromFile",&_m_BfromFile);
  setParameter("rotate",&_rotate);
  setParameter("B_theta",&_B_theta);
  setParameter("Vin",&_Vin);
  setParameter("rhoin",&_rhoin);
  setParameter("Tin",&_Tin);
  setParameter("Cin",&_Cin);

}

//////////////////////////////////////////////////////////////////////////////

CoronalInlet2Fin::~CoronalInlet2Fin()
{
}

//////////////////////////////////////////////////////////////////////////////

void CoronalInlet2Fin::setup()
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

void CoronalInlet2Fin::unsetup()
{
  SuperInlet::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void CoronalInlet2Fin::setGhostState(GeometricEntity *const face)
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
      CFLog(DEBUG_MIN, "CoronalInlet2Fin::setGhostState() => [" << varID << "] => " << (*initialValues)[idx] << " | " << (*innerState)[varID] << "\n");
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

  CFreal RSun = 6.9551e8; // m

  CFreal density_code = _rhoin; //1.67e-13;
  CFreal velocity_code = 2.2e-4/std::sqrt(1.256e-6*1.67e-13);
  CFreal pressure_code = std::pow(2.2e-4,2)/1.256e-6;
  CFreal B_code = 2.2e-4;

  CFreal RSS = 1.4953465e10; // m
  CFreal B0dip = 0.00022; // Tesla
  CFreal latI = 0.;
  const CFreal PI = MathTools::MathConsts::CFrealPi();
  const CFreal xG = ghostState->getCoordinates()[XX];// RSun;
  const CFreal xG_dimless = ghostState->getCoordinates()[XX]/RSun;
  const CFreal yG = ghostState->getCoordinates()[YY];//*RSun;
  const CFreal yG_dimless = ghostState->getCoordinates()[YY]/RSun;
  const CFreal zG = ghostState->getCoordinates()[ZZ];//*RSun;
  const CFreal zG_dimless = ghostState->getCoordinates()[ZZ]/RSun;
  const CFreal xI = innerState->getCoordinates()[XX];//*RSun;
  const CFreal xI_dimless = innerState->getCoordinates()[XX]/RSun;
  const CFreal yI = innerState->getCoordinates()[YY];//*RSun;
  const CFreal yI_dimless = innerState->getCoordinates()[YY]/RSun;
  const CFreal zI = innerState->getCoordinates()[ZZ];//*RSun;
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
  const CFreal rI = std::sqrt(xI*xI + yI*yI + zI*zI);
  CFreal rB = (rG + rI) / 2.0;
  const CFreal rhoG = std::sqrt(xG*xG + yG*yG);
  const CFreal rhoG_dimless = std::sqrt(xG_dimless*xG_dimless + yG_dimless*yG_dimless);
  const CFreal thetaG = std::atan2(std::sqrt(xG_dimless*xG_dimless + yG_dimless*yG_dimless),zG_dimless);
  const CFreal phiG = std::atan2(yG_dimless,xG_dimless);
  const CFreal rI_dimless = std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless + zI_dimless*zI_dimless);
  const CFreal thetaI = std::atan2(std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless),zI_dimless);
  const CFreal rhoI = std::sqrt(xI*xI + yI*yI);
  const CFreal rhoI_dimless = std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless);
  const CFreal phiI = std::atan2(yI_dimless,xI_dimless);


  


  // DENSITY

  CFreal m_proton = 1.6726219e-27; 
  CFreal m_electron = 9.10938356e-31; 
  CFreal me_to_mp = m_electron / m_proton;


  // We prescribe the density by prescribing the ion density and deriving the neutral density from the assumed concentration
  CFreal concentration = _Cin; 
  (*ghostState)[8] = 2.0*density_code - (*innerState)[8];
  (*ghostState)[9] = 2.0*density_code * concentration - (*innerState)[9];



  // MAGNETIC FIELD 

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

  /* Uncomment this for uniform radial outflow
  CFreal VxBoundaryn;
  CFreal VyBoundaryn;
  CFreal VzBoundaryn;

  CFreal VxBoundaryi;
  CFreal VyBoundaryi;
  CFreal VzBoundaryi;

  // For now, we prescribe a small uniform outflow in spherical coordinates
  CFreal VrBoundary = 1935.0; 
  CFreal VthetaBoundary = 0.0;
  CFreal VphiBoundary = 0.0;

  // This is translated to Cartesian coordinates
  CFreal VxBoundary= xBoundary/rBoundary*VrBoundary- yBoundary/rhoBoundary*VphiBoundary+ xBoundary*zBoundary/(rhoBoundary*rBoundary)*VthetaBoundary;
  CFreal VyBoundary= yBoundary/rBoundary*VrBoundary+ xBoundary/rhoBoundary*VphiBoundary+ yBoundary*zBoundary/(rhoBoundary*rBoundary)*VthetaBoundary;
  CFreal VzBoundary= zBoundary/rBoundary*VrBoundary- rhoBoundary/rBoundary*VthetaBoundary;

  // And then attributed to both fluids
  VxBoundaryi = VxBoundary;
  VyBoundaryi = VyBoundary;
  VzBoundaryi = VzBoundary;
  VxBoundaryn = VxBoundary;
  VyBoundaryn = VyBoundary;
  VzBoundaryn = VzBoundary;

  (*ghostState)[10] = (2.*VxBoundaryi- (*innerState)[10]);
  (*ghostState)[11] = (2.*VyBoundaryi- (*innerState)[11]);
  (*ghostState)[12] = (2.*VzBoundaryi- (*innerState)[12]);
  (*ghostState)[13] = (2.*VxBoundaryn- (*innerState)[13]);
  (*ghostState)[14] = (2.*VyBoundaryn- (*innerState)[14]);
  (*ghostState)[15] = (2.*VzBoundaryn- (*innerState)[15]);
  */
 
  // Just like in the case of ideal MHD, we determine the outflow direction from the direction of the magnetic field
  CFreal BxB = BxBoundary;
  CFreal ByB = ByBoundary;
  CFreal BzB = BzBoundary;
  CFreal BBmag = std::sqrt(BxB*BxB+ ByB*ByB+ BzB*BzB);
  CFreal VBmag = _Vin ; 
  CFreal BxuB = BxB/BBmag;
  CFreal ByuB = ByB/BBmag;
  CFreal BzuB = BzB/BBmag;
  CFreal VxBB = BxuB * VBmag;
  CFreal VyBB = ByuB * VBmag;
  CFreal VzBB = BzuB * VBmag;

  CFreal VrBB = xBoundary/rBoundary*VxBB+ yBoundary/rBoundary*VyBB + zBoundary/rBoundary*VzBB;
  CFreal VthetaBB = xBoundary*zBoundary/(rhoBoundary*rBoundary)*VxBB+ yBoundary*zBoundary/(rhoBoundary*rBoundary)*VyBB- rhoBoundary/rBoundary*VzBB;
  CFreal VphiBB = -yBoundary/rhoBoundary*VxBB + xBoundary/rhoBoundary*VyBB;

  if (VrBB < 0.0){
    VxBB = -VxBB;
    VyBB = -VyBB;
    VzBB = -VzBB;
  }

  // The same outflow is the attributed to both fluids (strong coupling assumed at the inner boundary)
  (*ghostState)[10] = (2.*VxBB- (*innerState)[10]);
  (*ghostState)[11] = (2.*VyBB- (*innerState)[11]);
  (*ghostState)[12] = (2.*VzBB- (*innerState)[12]);
  (*ghostState)[13] = (2.*VxBB- (*innerState)[13]);
  (*ghostState)[14] = (2.*VyBB- (*innerState)[14]);
  (*ghostState)[15] = (2.*VzBB- (*innerState)[15]);



  // TEMPERATURE

  CFreal T_Sun = 1.5e6;   // K
  CFreal mu_cor = 1.27;   // Mean molecular weight
  CFreal mH = 1.67e-27;   // Mass hydrogen
  CFreal BRef = 2.2e-4;   // T
  CFreal mu0 = 1.2566e-6;
  CFreal rhoRef = 1.67e-13;
  CFreal kB = 1.38e-23;

  // We set the temperature so that pressure is similar to what we had in ideal MHD
  CFreal pBoundary = 0.25801090625/2.0 * _rhoin/1.67e-13;
  CFreal TBoundary = pBoundary/ (2.0 * density_code / mH) / kB;

  // We assume that the temperature at the inner boundary is the same for ions and neutrals
  (*ghostState)[16] = (2.* _Tin - (*innerState)[16]); 
  (*ghostState)[17] = (2.* _Tin - (*innerState)[17]); 




 
  // ELECTRIC FIELD
 
  // for it to be zero at the boundary 
  (*ghostState)[3] = -(*innerState)[3]; 
  (*ghostState)[4] = -(*innerState)[4]; 
  (*ghostState)[5] = -(*innerState)[5]; 



  // DIVERGENCE CLEANING

  (*ghostState)[6] = -(*innerState)[6]; 
  (*ghostState)[7] = -(*innerState)[7]; 



}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
