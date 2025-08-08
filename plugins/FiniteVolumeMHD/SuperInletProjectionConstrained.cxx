#include "FiniteVolume/FiniteVolume.hh"
#include "SuperInletProjectionConstrained.hh"
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

MethodCommandProvider<SuperInletProjectionConstrained, CellCenterFVMData, FiniteVolumeModule> 
superInletProjectionConstrainedFVMCCProvider("SuperInletProjectionConstrainedFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SuperInletProjectionConstrained::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFuint> >("ProjectionIDs","IDs corresponding to projection fields.");
  options.addConfigOption< vector<CFuint> >("VarIDs","IDs of the variables from which values are read by file");
  options.addConfigOption< CFreal >("pBC","pressure at the boundary");
  options.addConfigOption< CFreal >("rhoBC","density at the boundary");
  options.addConfigOption< CFreal >("VrBC","radial velocity at the boundary");
  options.addConfigOption< CFint >("rotation","rotation, 0 or 1");
}
      
//////////////////////////////////////////////////////////////////////////////

SuperInletProjectionConstrained::SuperInletProjectionConstrained(const std::string& name) :
  SuperInlet(name),
  m_mapGeoToTrs(),
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

  _VrBC = 1340.; //848.15;

  //_VrBC = vector<CFuint>();
  setParameter("VrBC",&_VrBC);

  _rotation = 0; 

  //_VrBC = vector<CFuint>();
  setParameter("rotation",&_rotation);
}

//////////////////////////////////////////////////////////////////////////////

SuperInletProjectionConstrained::~SuperInletProjectionConstrained()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletProjectionConstrained::setup()
{
  SuperInlet::setup();

  if (m_varIDs.size() > 0) {m_BfromFile = true;}
  
  m_mapGeoToTrs =
    MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletProjectionConstrained::unsetup()
{
  SuperInlet::unsetup();
}


void SuperInletProjectionConstrained::setGhostState(GeometricEntity *const face)
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
	//CFLog(INFO, "faceIdx="<< faceIdx<<"\n");
    const string name = getCurrentTRS()->getName();
	//CFLog(INFO, "name="<< name<<"\n"); //Inlet
    SafePtr<RealVector> initialValues = m_initialSolutionMap.find(name);
    const CFuint nbVars = m_initialSolutionIDs.size();
	//CFLog(INFO, "nbVars="<<nbVars<<"\n");
    const CFuint startID = faceIdx*nbVars;
	//CFLog(INFO, "startID="<<startID<<"\n");
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
  //m_BfromFile = false;
  CFreal BrFromFile = 0.;
  if (m_initialSolutionMap.size() > 0 && m_BfromFile == true) {
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
	//CFLog(INFO, "nbNodesInFace="<<nbNodesInFace<<"\n");
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


  // Determine the latitude of the ghost cell:
  // Only needed when using the dead-zone
  if (thetaG > -PI && thetaG < -PI*0.5) {
     latG = std::abs(thetaG) - PI*0.5;
  } else if (thetaG > -PI*0.5 && thetaG < 0) {
     latG = PI*0.5 - std::abs(thetaG);
  } else if (thetaG > 0. && thetaG < PI*0.5) {
     latG = PI*0.5 - thetaG;
  } else if (thetaG > PI*0.5 && thetaG < PI) {
     latG = thetaG - PI*0.5;
  } 

  //else {
  //  CFLog(INFO, "Error: value of theta for the point in question outside expected range\n");
  //}




  CFreal T_Sun = 1.5e6;   // K
  CFreal mu_cor = 1.27;   // Mean molecular weight
  CFreal mH = 1.67e-27;   // Mass hydrogen
  CFreal Bref = 2.2e-4;   // T
  CFreal mu0 = 1.2566e-6;
  CFreal rhoref = 1.67e-13;
  CFreal kB = 1.38e-23;
  CFreal pref = 0.03851;
  CFreal vref = 4.8e5;
  CFreal Veps = 2.0e-6;

  //===== M A G N E T I C  F I E L D   B O U N D A R Y   C O N D I T I O N ===============
 

//  CFreal BrBoundary_dimless; // only initialization, it is overwritten just below
//  BrBoundary_dimless = xBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[0] + yBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[1] + zBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[2];
  
  CFreal BrBoundary_dimless = BrFromFile;
  if (!m_BfromFile) {
   BrBoundary_dimless = xBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[0] + yBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[1] + zBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[2];
  }
// CFLog(INFO, "m_BfromFile[" << m_BfromFile << "], Br from File vs BrBoundary_dimless = " << BrFromFile << " ? " << BrBoundary_dimless << " \n"); 


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
  //-------->> Mark 2025.03.30----------------------------
  //BrBoundary_dimless = (BrBoundary_dimless + BrI_dimless) / 2.0;
  //---------- Mark 2025.03.30<<--------------------------

  //BphiBoundary_dimless=0.0;
  //BthetaBoundary_dimless=0.0;


  CFreal BxBoundary_dimless = xBoundary_dimless/rBoundary_dimless*BrBoundary_dimless - yBoundary_dimless/rhoBoundary_dimless*BphiBoundary_dimless + xBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*BthetaBoundary_dimless;
  CFreal ByBoundary_dimless = yBoundary_dimless/rBoundary_dimless*BrBoundary_dimless + xBoundary_dimless/rhoBoundary_dimless*BphiBoundary_dimless + yBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*BthetaBoundary_dimless;
  CFreal BzBoundary_dimless = zBoundary_dimless/rBoundary_dimless*BrBoundary_dimless - rhoBoundary_dimless/rBoundary_dimless*BthetaBoundary_dimless;


  (*ghostState)[4] = 2*BxBoundary_dimless - (*innerState)[4];
  (*ghostState)[5] = 2*ByBoundary_dimless - (*innerState)[5];
  (*ghostState)[6] = 2*BzBoundary_dimless - (*innerState)[6];




  //===== V E L O C I T Y   B O U N D A R Y   C O N D I T I O N ===============

  // Read the inner state velocity need for extrapolation radially inwards:
  CFreal VxI_dimless = (*innerState)[1];
  CFreal VyI_dimless = (*innerState)[2];
  CFreal VzI_dimless = (*innerState)[3];
  CFreal VrI_dimless = xI_dimless/rI_dimless*VxI_dimless + yI_dimless/rI_dimless*VyI_dimless + zI_dimless/rI_dimless*VzI_dimless;
  CFreal VthetaI_dimless = xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*VxI_dimless + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*VyI_dimless - rhoI_dimless/rI_dimless*VzI_dimless;
  CFreal VphiI_dimless = -yI_dimless/rhoI_dimless*VxI_dimless + xI_dimless/rhoI_dimless*VyI_dimless;

  CFreal VrBoundaryUpdate_dimless = VrI_dimless; 
  CFreal VthetaBoundaryUpdate_dimless = VthetaI_dimless; 
  CFreal VphiBoundaryUpdate_dimless = VphiI_dimless;

  if (_rotation == 1) {
    VphiBoundaryUpdate_dimless = 2.*rBoundary_dimless*std::sin(thetaBoundary_dimless)*3.86e-3;
  }
  else if (_rotation == -1) {
    VphiBoundaryUpdate_dimless = -2.*rBoundary_dimless*std::sin(thetaBoundary_dimless)*3.86e-3;
  }
  else {  
    VphiBoundaryUpdate_dimless = 0.0;  
  }


  CFreal VxBoundaryUpdate_dimless = xBoundary_dimless/rBoundary_dimless*VrBoundaryUpdate_dimless - yBoundary_dimless/rhoBoundary_dimless*VphiBoundaryUpdate_dimless + xBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*VthetaBoundaryUpdate_dimless;
  CFreal VyBoundaryUpdate_dimless = yBoundary_dimless/rBoundary_dimless*VrBoundaryUpdate_dimless + xBoundary_dimless/rhoBoundary_dimless*VphiBoundaryUpdate_dimless + yBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*VthetaBoundaryUpdate_dimless;
  CFreal VzBoundaryUpdate_dimless = zBoundary_dimless/rBoundary_dimless*VrBoundaryUpdate_dimless - rhoBoundary_dimless/rBoundary_dimless*VthetaBoundaryUpdate_dimless;

  CFreal Bmag = std::sqrt(BxBoundary_dimless*BxBoundary_dimless + ByBoundary_dimless*ByBoundary_dimless + BzBoundary_dimless*BzBoundary_dimless);
  CFreal Vmag = std::sqrt(VxBoundaryUpdate_dimless*VxBoundaryUpdate_dimless + VyBoundaryUpdate_dimless*VyBoundaryUpdate_dimless + VzBoundaryUpdate_dimless*VzBoundaryUpdate_dimless);

  Vmag = std::sqrt(VxI_dimless*VxI_dimless + VyI_dimless*VyI_dimless + VzI_dimless*VzI_dimless);
  CFreal Vmagold=Vmag;
  //Vmag = VrI_dimless;
  //>> OK 
  CFreal Vmagmax = 1.0e3;
  //<<
  //if (Vmag * vref > Vmagmax) {
  //  Vmag = Vmagmax / vref;
  //}

  CFreal smooth_factor = 0.0;
  if(Vmag>0.0){
  CFreal dist = (Vmag * vref - Vmagmax*1.0)/1.0e0;
  smooth_factor = 0.5 + 0.5 * std::tanh(dist * 3.141592654);
  CFreal Vmagnew = smooth_factor * Vmagmax + Vmag * vref * (1.0-smooth_factor);
  Vmag = Vmagnew / vref;
  }
  else{
  CFreal dist = (-Vmag * vref - Vmagmax*0.1)/2.0e0;
  smooth_factor = 0.5 + 0.5 * std::tanh(dist * 3.141592654);
  CFreal Vmagnew = smooth_factor * (-Vmagmax) + Vmag * vref * (1.0-smooth_factor);  
  Vmag = Vmagnew / vref;
  }

  CFreal Bxu = BxBoundary_dimless/Bmag;
  CFreal Byu = ByBoundary_dimless/Bmag;
  CFreal Bzu = BzBoundary_dimless/Bmag;

  CFreal sign = 1.0;
  CFreal Vdotg = BxBoundary_dimless*VxBoundaryUpdate_dimless + ByBoundary_dimless*VyBoundaryUpdate_dimless + BzBoundary_dimless*VzBoundaryUpdate_dimless;

  if (Vdotg < 0.0) {
    sign = -1.0;
  }

  //>> mark 2024.06.17 comment for time-evolving coronal simulation
  //CFreal VxB = Bxu * Vmag * sign;
  //CFreal VyB = Byu * Vmag * sign;
  //CFreal VzB = Bzu * Vmag * sign;
  //<

  //>> mark 2024.06.17 comment for quasisteady state coronal simulation  
  CFreal VxB = VxBoundaryUpdate_dimless/Vmagold*Vmag;
  CFreal VyB = VyBoundaryUpdate_dimless/Vmagold*Vmag;
  CFreal VzB = VzBoundaryUpdate_dimless/Vmagold*Vmag;
  //<<

  //CFreal Vr=xBoundary_dimless/rBoundary_dimless*VxB+yBoundary_dimless/rBoundary_dimless*VyB+zBoundary_dimless/rBoundary_dimless*VzB;
  //if(Vr<0.0){
  //	  VxB=0.0;
  //	  VyB=0.0;
  //  VzB=0.0;
  //}



  (*ghostState)[1] = 2.*VxB -(*innerState)[1]; 
  (*ghostState)[2] = 2.*VyB -(*innerState)[2]; 
  (*ghostState)[3] = 2.*VzB -(*innerState)[3]; 



  //vA = B / SQRT(rho mu)
  //vA^2 = B^2 / (rho mu)
  // rho mu = B^2 / vA^2 --> rho = 1/mu B^2 / vA^2


  CFreal densityBoundary_dimless = _rhoBC; 
  CFreal PressureBoundary_dimless = _pBC;
  //---->> Mark 2025.03.28---------------------------------------------
  //CFreal Temperature_ratio = PressureBoundary_dimless / densityBoundary_dimless;
  //------ Mark 2025.03.28<<-------------------------------------------

  CFreal maxvA = 3.0e6;   //3.0e8; //3.0e6; //6.0e6; //
  CFreal vA = Bmag*Bref / pow((densityBoundary_dimless * rhoref * mu0),0.5);

  smooth_factor = 0.0;
  CFreal dist = (vA - maxvA)/2.0e3;
  smooth_factor = 0.5 + 0.5 * std::tanh(dist * 3.141592654);

  CFreal densityG = smooth_factor * (Bmag*Bref)*(Bmag*Bref) / (maxvA*maxvA) / mu0 + densityBoundary_dimless * rhoref * (1.0-smooth_factor); 
  (*ghostState)[0] = 2.0*(densityG / rhoref) - (*innerState)[0];

  //pressurefactor = (densityG/ 1.67e-13) / _rhoBC;


  //if (vA > maxvA) {
  //  CFreal rho = (Bmag*Bref)*(Bmag*Bref) / (maxvA*maxvA) / mu0;
  //  densityBoundary_dimless = rho / rhoref;
  //}
  //CFreal densityG_dimless = 2.*densityBoundary_dimless - (*innerState)[0];
  //(*ghostState)[0] = densityBoundary_dimless; // densityG_dimless;



  //  CFreal densityG = smooth_factor * BBmagdim*BBmagdim/vAmax/vAmax/1.257e-6 + _rhoBC * 1.67e-13 * (1.0-smooth_factor); // /1.67e-13;
  //  (*ghostState)[0] = 2.0*(densityG / 1.67e-13) - (*innerState)[0];

  //  pressurefactor = (densityG/ 1.67e-13) / _rhoBC;





  CFreal minPlasmaBeta = 0.001; // 0.02;
  CFreal maxPlasmaBeta = 1.0;
  CFreal pmag = (Bmag*Bref)*(Bmag*Bref)/2./mu0;
  
  //dist = VrBoundaryUpdate_dimless/Vmagold*Vmag/Veps;
  //smooth_factor = 0.5 + 0.5 * std::tanh(dist * 3.141592654);
  //CFreal densityBoundary_dimless = _rhoBC*smooth_factor+(*innerState)[0]*(1.0-smooth_factor);   
  //CFreal PressureBoundary_dimless = _pBC; 
  CFreal pth = PressureBoundary_dimless * pref;

  CFreal plasmaBeta = pth/pmag;

  smooth_factor = 0.0;
  dist = (minPlasmaBeta - plasmaBeta)/2.0e-6; //if current beta is 0.019, then dist is positive, smooth factor is close to 1
  smooth_factor = 0.5 + 0.5 * std::tanh(dist * 3.141592654);  //if the current beta is  0.3, then dist is negative and so smooth factor is 0

  CFreal pthG = smooth_factor * pmag * minPlasmaBeta + pth * (1.0-smooth_factor);
  (*ghostState)[7] = 2.0*(pthG / pref) - (*innerState)[7];



  /*if (pth / pmag > minPlasmaBeta) {
    pth = pmag * minPlasmaBeta;
    PressureBoundary_dimless = pth / pref;}
  if (pth / pmag < minPlasmaBeta) {
    pth = minPlasmaBeta * pmag;
    PressureBoundary_dimless =  pth / pref;}

  (*ghostState)[7] = PressureBoundary_dimless; // 2.*PressureBoundary_dimless - (*innerState)[7];
  */

  
    (*ghostState)[8] = (*innerState)[8];
  

}

//////////////////////////////////////////////////////////////////////////////
//Mark 2023/12/14
void SuperInletProjectionConstrained::preProcess()
{
  SuperInlet::preProcess();
  CFLog(INFO, "SuperInletProjectionConstrained::preProcess() => START\n");  
  // AL: only for unsteady runs (DT > 0.)  
  if (SubSystemStatusStack::getActive()->getDT() > 0.) {
   if (SubSystemStatusStack::getActive()->getSubIter() == 0) {
    SafePtr<NodalStatesExtrapolator<CellCenterFVMData> > nse =
      this->getMethodData().getNodalStatesExtrapolator();
    nse->extrapolateVarsFromFileInTime();
   }
  }
  CFLog(INFO, "SuperInletProjectionConstrained::preProcess() => END\n");
}
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
