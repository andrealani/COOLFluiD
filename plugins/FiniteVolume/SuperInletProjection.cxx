#include "FiniteVolume/FiniteVolume.hh"
#include "SuperInletProjection.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SuperInletProjection, CellCenterFVMData, FiniteVolumeModule> 
superInletProjectionFVMCCProvider("SuperInletProjectionFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SuperInletProjection::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFuint> >("ProjectionIDs","IDs corresponding to projection fields.");
  options.addConfigOption< vector<CFuint> >("VarIDs","IDs of the variables from which values are read by file");
  options.addConfigOption< CFreal >("pBC","pressure at the boundary");
  options.addConfigOption< CFreal >("rhoBC","density at the boundary");
  options.addConfigOption< CFreal >("VrBC","radial velocity at the boundary");
  options.addConfigOption< CFint >("rotation","rotation, 0 or 1");
}
      
//////////////////////////////////////////////////////////////////////////////

SuperInletProjection::SuperInletProjection(const std::string& name) :
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


  _pBC = 0.108; //*8.0;

  //_pBC = vector<CFuint>();
  setParameter("pBC",&_pBC);

  _rhoBC = 1.0;

  //_rhoBC = vector<CFuint>();
  setParameter("rhoBC",&_rhoBC);

  _VrBC = 1935.07; //848.15;

  //_VrBC = vector<CFuint>();
  setParameter("VrBC",&_VrBC);

  _rotation = 0; 

  //_VrBC = vector<CFuint>();
  setParameter("rotation",&_rotation);
}

//////////////////////////////////////////////////////////////////////////////

SuperInletProjection::~SuperInletProjection()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletProjection::setup()
{
  SuperInlet::setup();

  if (m_varIDs.size() > 0) {m_BfromFile = true;}
  
  m_mapGeoToTrs =
    MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");
  
  // build the m_mapTrs2Twall storage
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

void SuperInletProjection::unsetup()
{
  SuperInlet::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletProjection::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // during initialization phase we store the initial solution values to be used as BC
  // This BC (the PFSS solution fixed on the inlet) was suggested by Jon Linker and
  // is evaluated by default. When using Jens' or Dana's BCs for the magnetic field
  // the B-field values in the ghost cells are overwritten.

  RealVector B_PFSS_dimless(3);

  if (m_initialSolutionMap.size() > 0) {
    /// map faces to corresponding TRS and index inside that TRS
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

      B_PFSS_dimless[i] = (*initialValues)[idx]; // save the Poisson PFSS solution in a vector

      (*ghostState)[varID] = 2.*(*initialValues)[idx] - (*innerState)[varID];
      CFLog(DEBUG_MIN, "SuperInletProjection::setGhostState() => [" << varID << "] => " << (*initialValues)[idx] << " | " << (*innerState)[varID] << "\n");
    }
  }

  // this interpolates Br directly from the magnetogram file
  // static CFreal maxBr = -1e-14;
  CFreal BrFromFile = 1; //0.;
  //std::cout << m_BfromFile << " m_BfromFile " << "\n";

  m_BfromFile = false;
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

    // maxBr = std::max(maxBr, BrFromFile);
    // CFLog(INFO, "Br from file #### varID = " << varID << ", maxBr = " << maxBr << "\n");

  }

  CFreal latG = 0.;
  CFreal RSun = 6.9551e8; // m
  CFreal RSS = 1.4953465e10; // m
  CFreal B0dip = 0.00022; // Tesla
  CFreal latI = 0.;
  const CFreal PI = MathTools::MathConsts::CFrealPi();
  const CFreal xG = ghostState->getCoordinates()[XX]*RSun;
  const CFreal xG_dimless = ghostState->getCoordinates()[XX];
  const CFreal yG = ghostState->getCoordinates()[YY]*RSun;
  const CFreal yG_dimless = ghostState->getCoordinates()[YY];
  const CFreal zG = ghostState->getCoordinates()[ZZ]*RSun;
  const CFreal zG_dimless = ghostState->getCoordinates()[ZZ];
  const CFreal xI = innerState->getCoordinates()[XX]*RSun;
  const CFreal xI_dimless = innerState->getCoordinates()[XX];
  const CFreal yI = innerState->getCoordinates()[YY]*RSun;
  const CFreal yI_dimless = innerState->getCoordinates()[YY];
  const CFreal zI = innerState->getCoordinates()[ZZ]*RSun;
  const CFreal zI_dimless = innerState->getCoordinates()[ZZ];
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



  //===== D E N S I T Y   B O U N D A R Y   C O N D I T I O N =================
  

  CFreal densityBoundary_dimless = _rhoBC;
  CFreal densityG_dimless = 2.*densityBoundary_dimless - (*innerState)[0];
  (*ghostState)[0] = densityG_dimless;


  //===== M A G N E T I C  F I E L D   B O U N D A R Y   C O N D I T I O N ===============
 

  CFreal BrBoundary_dimless; // only initialization, it is overwritten just below

  BrBoundary_dimless = xBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[0] + yBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[1] + zBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[2];

  CFreal BthetaBoundary_dimless = xBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*B_PFSS_dimless[0] + yBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*B_PFSS_dimless[1] - rhoBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[2];

  CFreal BxI_dimless = (*innerState)[4];
  CFreal ByI_dimless = (*innerState)[5];
  CFreal BzI_dimless = (*innerState)[6];
  CFreal BrI_dimless = xI_dimless/rI_dimless*BxI_dimless + yI_dimless/rI_dimless*ByI_dimless + zI_dimless/rI_dimless*BzI_dimless;
  CFreal BthetaI_dimless = xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*BxI_dimless + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*ByI_dimless - rhoI_dimless/rI_dimless*BzI_dimless;
  CFreal BphiI_dimless = -yI_dimless/rhoI_dimless*BxI_dimless + xI_dimless/rhoI_dimless*ByI_dimless;
  
  CFreal BphiG_dimless = BphiI_dimless;
  CFreal BphiBoundary_dimless = (BphiG_dimless + BphiI_dimless)/2.0;
  BthetaBoundary_dimless = (BthetaI_dimless + BthetaI_dimless)/2.0;



  CFreal BxBoundary_dimless = xBoundary_dimless/rBoundary_dimless*BrBoundary_dimless - yBoundary_dimless/rhoBoundary_dimless*BphiBoundary_dimless + xBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*BthetaBoundary_dimless;
  CFreal ByBoundary_dimless = yBoundary_dimless/rBoundary_dimless*BrBoundary_dimless + xBoundary_dimless/rhoBoundary_dimless*BphiBoundary_dimless + yBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*BthetaBoundary_dimless;
  CFreal BzBoundary_dimless = zBoundary_dimless/rBoundary_dimless*BrBoundary_dimless - rhoBoundary_dimless/rBoundary_dimless*BthetaBoundary_dimless;


  (*ghostState)[4] = 2*BxBoundary_dimless - (*innerState)[4];
  (*ghostState)[5] = 2*ByBoundary_dimless - (*innerState)[5];
  (*ghostState)[6] = 2*BzBoundary_dimless - (*innerState)[6];




  //===== V E L O C I T Y   B O U N D A R Y   C O N D I T I O N ===============

  CFreal BxB = BxBoundary_dimless;
  CFreal ByB = ByBoundary_dimless;
  CFreal BzB = BzBoundary_dimless; 
  CFreal BBmag = std::sqrt(BxB*BxB+ ByB*ByB+ BzB*BzB);
  CFreal VBmag = _VrBC/(2.2e-4/sqrt(1.2566e-6*1.67e-13)); //1935.07/(2.2e-4/sqrt(1.2566e-6*1.67e-13));
  CFreal BxuB = BxB/BBmag;
  CFreal ByuB = ByB/BBmag;
  CFreal BzuB = BzB/BBmag;
  CFreal VxBB = BxuB * VBmag;
  CFreal VyBB = ByuB * VBmag;
  CFreal VzBB = BzuB * VBmag;

  //Compute the spherical components from the parallel velocity 
  CFreal VrBB = xBoundary_dimless/rBoundary_dimless*VxBB+ yBoundary_dimless/rBoundary_dimless*VyBB + zBoundary_dimless/rBoundary_dimless*VzBB;
  CFreal VthetaBB = xBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*VxBB+ yBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*VyBB- rhoBoundary_dimless/rBoundary_dimless*VzBB;
  CFreal VphiBB = -yBoundary_dimless/rhoBoundary_dimless*VxBB + xBoundary_dimless/rhoBoundary_dimless*VyBB;

  //If the radial component would indicate inflow, switch the orientation
  if (VrBB < 0.0){
    VxBB = -VxBB;
    VyBB = -VyBB;
    VzBB = -VzBB;
  }

  //With the fixed-orientation components, recalculate the spherical components
  VrBB = xBoundary_dimless/rBoundary_dimless*VxBB+ yBoundary_dimless/rBoundary_dimless*VyBB + zBoundary_dimless/rBoundary_dimless*VzBB;
  VthetaBB = xBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*VxBB+ yBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*VyBB- rhoBoundary_dimless/rBoundary_dimless*VzBB;
  VphiBB = -yBoundary_dimless/rhoBoundary_dimless*VxBB + xBoundary_dimless/rhoBoundary_dimless*VyBB;

  //Add the rotation motion
 
  if (_rotation == 1){
    VphiBB = VphiBB + rBoundary_dimless*std::sin(thetaBoundary_dimless)*3.86e-3;} 

  //Recompute Cartesian components
  CFreal VxBBR = xBoundary_dimless/rBoundary_dimless*VrBB - yBoundary_dimless/rhoBoundary_dimless*VphiBB + xBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*VthetaBB;
  CFreal VyBBR = yBoundary_dimless/rBoundary_dimless*VrBB + xBoundary_dimless/rhoBoundary_dimless*VphiBB + yBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*VthetaBB;
  CFreal VzBBR = zBoundary_dimless/rBoundary_dimless*VrBB - rhoBoundary_dimless/rBoundary_dimless*VthetaBB;  


  (*ghostState)[1] = 2.0*VxBBR - (*innerState)[1];
  (*ghostState)[2] = 2.0*VyBBR - (*innerState)[2];
  (*ghostState)[3] = 2.0*VzBBR - (*innerState)[3];


 
 


  //===== P R E S S U R E   B O U N D A R Y   C O N D I T I O N ===============
  // Temperature kept constant at 1.5e6 K
  // Pressure at inner boundary:

  CFreal T_Sun = 1.5e6;   // K
  CFreal mu_cor = 1.27;   // Mean molecular weight
  CFreal mH = 1.67e-27;   // Mass hydrogen
  CFreal BRef = 2.2e-4;   // T
  CFreal mu0 = 1.2566e-6;
  CFreal rhoRef = 1.67e-13;
  CFreal kB = 1.38e-23;


  CFreal PressureBoundary_dimless = _pBC;


  (*ghostState)[7] = 2.*PressureBoundary_dimless - (*innerState)[7];


  


  //===== P S I   B O U N D A R Y   C O N D I T I O N S =======================
    (*ghostState)[8] = (*innerState)[8];
  

}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
