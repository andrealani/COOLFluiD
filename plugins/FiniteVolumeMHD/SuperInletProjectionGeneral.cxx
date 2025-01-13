#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "FiniteVolumeMHD/SuperInletProjectionGeneral.hh"
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

MethodCommandProvider<SuperInletProjectionGeneral, CellCenterFVMData, FiniteVolumeMHDModule> 
superInletProjectionGeneralFVMCCProvider("SuperInletProjectionGeneralFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SuperInletProjectionGeneral::defineConfigOptions(Config::OptionList& options)
{
  // Here we have definitions of the flags one can set via the CFcase file    
  options.addConfigOption< vector<CFuint> >("ProjectionIDs","IDs corresponding to projection fields.");
  options.addConfigOption< vector<CFuint> >("VarIDs","IDs of the variables from which values are read by file");
  options.addConfigOption< CFreal >("pBC","pressure at the boundary");
  options.addConfigOption< CFreal >("rhoBC","density at the boundary");
  options.addConfigOption< CFreal >("VrBC","radial velocity at the boundary");
  options.addConfigOption< CFint >("rotation","rotation, 0 or 1");
  options.addConfigOption< CFreal >("betamin","limiting beta constraint");
  options.addConfigOption< CFreal >("vAmax","limiting vA constraint");
  options.addConfigOption< bool >("Nonhomogeneous","homogeneous or nonhomogenous BC prescription");
}
      
//////////////////////////////////////////////////////////////////////////////

SuperInletProjectionGeneral::SuperInletProjectionGeneral(const std::string& name) :
  SuperInlet(name),
  m_mapGeoToTrs(),
  m_BfromFile(false)
{
  addConfigOptionsTo(this);
  
  _projectionIDs = vector<CFuint>();
  setParameter("ProjectionIDs",&_projectionIDs);
  
  m_varIDs = vector<CFuint>();
  setParameter("VarIDs",&m_varIDs);
  
  // Adimensional boundary pressure
  _pBC = 0.108; 
  setParameter("pBC",&_pBC);
  
  // Adimensional boundary density
  _rhoBC = 1.0;
  setParameter("rhoBC",&_rhoBC);
  
  // Dimensional boundary radial velocity
  _VrBC = 1340.; 
  setParameter("VrBC",&_VrBC);

  // Rotation yes or no
  _rotation = 0; 
  setParameter("rotation",&_rotation);
  
  // Minimum allowed plasma beta (only matters if nonhomogeneous BC)
  _betamin = 0.02;
  setParameter("betamin",&_betamin);

  // Maximum allowed Alfven speed (only matters if nonhomogeneous BC)
  _vAmax = 2e6;
  setParameter("vAmax",&_vAmax);
  
  // Homogenenous or non-homogeneous boundary condition for p and rho
  _nonhomogeneous = false;
  setParameter("Nonhomogeneous",&_nonhomogeneous);
}

//////////////////////////////////////////////////////////////////////////////

SuperInletProjectionGeneral::~SuperInletProjectionGeneral()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletProjectionGeneral::setup()
{
  SuperInlet::setup();

  if (m_varIDs.size() > 0) {m_BfromFile = true;}
  
  m_mapGeoToTrs =
    MeshDataStack::getActive()->getMapGeoToTrs("MapFacesToTrs");
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletProjectionGeneral::unsetup()
{
  SuperInlet::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletProjectionGeneral::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

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
      CFLog(DEBUG_MIN, "SuperInletProjectionGeneral::setGhostState() => [" << varID << "] => " << (*initialValues)[idx] << " | " << (*innerState)[varID] << "\n");
    }
  }

  // this interpolates Br directly from the magnetogram file
  // static CFreal maxBr = -1e-14;
  CFreal BrFromFile = 1; //0.;
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
  const CFreal xB_dimless = xBoundary/RSun;
  const CFreal yB_dimless = yBoundary/RSun;
  const CFreal zB_dimless = zBoundary/RSun;
  const CFreal rB_dimless = rBoundary_dimless;
  const CFreal rhoB_dimless = rhoBoundary_dimless;

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


  CFreal T_Sun = 1.5e6;   
  CFreal mu_cor = 1.27;  
  CFreal mH = 1.67e-27; 
  CFreal Bref = 2.2e-4;
  CFreal mu0 = 1.2566e-6;
  CFreal rhoref = 1.67e-13;
  CFreal kB = 1.38e-23;
  CFreal pref = 0.03851;
  CFreal vref = 4.8e5;


  //===== M A G N E T I C  F I E L D   B O U N D A R Y   C O N D I T I O N ===============
 

  // Determine the radial magnetic field on the inner boundary from the PFSS
  CFreal BrBoundary_dimless = xBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[0] + yBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[1] + zBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[2];

  // We also determine the azimuthal magnetic field from the PFSS, but we overwrite it later on
  CFreal BthetaBoundary_dimless = xBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*B_PFSS_dimless[0] + yBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*B_PFSS_dimless[1] - rhoBoundary_dimless/rBoundary_dimless*B_PFSS_dimless[2];


  // We determine the current inner state of the B-field in Cartesian coordinates and then transform these into spherical coordinates
  CFreal BxI_dimless = (*innerState)[4];
  CFreal ByI_dimless = (*innerState)[5];
  CFreal BzI_dimless = (*innerState)[6];
  CFreal BrI_dimless = xI_dimless/rI_dimless*BxI_dimless + yI_dimless/rI_dimless*ByI_dimless + zI_dimless/rI_dimless*BzI_dimless;
  CFreal BthetaI_dimless = xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*BxI_dimless + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*ByI_dimless - rhoI_dimless/rI_dimless*BzI_dimless;
  CFreal BphiI_dimless = -yI_dimless/rhoI_dimless*BxI_dimless + xI_dimless/rhoI_dimless*ByI_dimless;
  
  // We allows the phi and theta components of the magnetic field to evolve freely (set the boundary value to the current inner state value) 
  CFreal BphiBoundary_dimless = BphiI_dimless; 
  BthetaBoundary_dimless = BthetaI_dimless; 

  // With the Br determined from the initial Br in PFSS and free evolving Bphi and Btheta, we transform these back to Cartesian coordinates and impose them onto the boundary 
  CFreal BxBoundary_dimless = xBoundary_dimless/rBoundary_dimless*BrBoundary_dimless - yBoundary_dimless/rhoBoundary_dimless*BphiBoundary_dimless + xBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*BthetaBoundary_dimless;
  CFreal ByBoundary_dimless = yBoundary_dimless/rBoundary_dimless*BrBoundary_dimless + xBoundary_dimless/rhoBoundary_dimless*BphiBoundary_dimless + yBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*BthetaBoundary_dimless;
  CFreal BzBoundary_dimless = zBoundary_dimless/rBoundary_dimless*BrBoundary_dimless - rhoBoundary_dimless/rBoundary_dimless*BthetaBoundary_dimless;
  (*ghostState)[4] = 2*BxBoundary_dimless - (*innerState)[4];
  (*ghostState)[5] = 2*ByBoundary_dimless - (*innerState)[5];
  (*ghostState)[6] = 2*BzBoundary_dimless - (*innerState)[6];




  //===== V E L O C I T Y   B O U N D A R Y   C O N D I T I O N ===============

  // For the velocity boundary condition, we ensure that the plasma escapes with a certain speed while following the magnetic field lines

  // We read the current inner state velocity need for extrapolation radially inwards
  CFreal VxI_dimless = (*innerState)[1];
  CFreal VyI_dimless = (*innerState)[2];
  CFreal VzI_dimless = (*innerState)[3];

  // We define the magnetic field strength in the current location (needed later for vector normalisation)
  CFreal Bmag = std::sqrt(BxBoundary_dimless*BxBoundary_dimless + ByBoundary_dimless*ByBoundary_dimless + BzBoundary_dimless*BzBoundary_dimless);

  // And we also define the current velocity magnitude
  CFreal Vmag = std::sqrt(VxI_dimless*VxI_dimless + VyI_dimless*VyI_dimless + VzI_dimless*VzI_dimless);

  CFreal smooth_factor = 0.0;
  CFreal dist = 0.0;


  // Now we have two options. If we have homogenenous boundary conditions, the speed of the outflow is solely determined by _VrBC defined by the user. 
  // If it is nonhomogeneous, it may be more appropriate to give more freedom to the plasma and only set a maximum speed allowed. 
  if (_nonhomogeneous) {

    CFreal Vmagmax = 1e5;
    smooth_factor = 0.0;
    dist = (Vmag * vref - Vmagmax)/2e4;
    smooth_factor = 0.5 + 0.5 * std::tanh(dist * 3.141592654);

    CFreal Vmagnew = smooth_factor * Vmagmax + Vmag * vref * (1.0-smooth_factor);
    Vmag = Vmagnew / vref;
    }

  else {
    Vmag = _VrBC / vref;
    }


  // Next, we compute the magnetic field unit vector
  CFreal Bxu = BxBoundary_dimless/Bmag;
  CFreal Byu = ByBoundary_dimless/Bmag;
  CFreal Bzu = BzBoundary_dimless/Bmag;

  CFreal sign = 1.0;
  CFreal Vdotg = BxBoundary_dimless*VxI_dimless + ByBoundary_dimless*VyI_dimless + BzBoundary_dimless*VzI_dimless;

  if (Vdotg < 0.0) {
    sign = -1.0;
  }

  // We Ensure that the prescribed velocity vector follows the magnetic field lines
  CFreal VxB = Bxu * Vmag * sign;
  CFreal VyB = Byu * Vmag * sign;
  CFreal VzB = Bzu * Vmag * sign;

  // From this, we transfer the Cartesian results to spherical in order to add rotation if necessary
  CFreal VrB = xB_dimless/rB_dimless*VxB + yB_dimless/rB_dimless*VyB + zB_dimless/rB_dimless*VzB;
  CFreal VthetaB = xB_dimless*zB_dimless/(rhoB_dimless*rB_dimless)*VxB + yB_dimless*zB_dimless/(rhoB_dimless*rB_dimless)*VyB - rhoB_dimless/rB_dimless*VzB;
  CFreal VphiB = -yB_dimless/rhoB_dimless*VxB + xB_dimless/rhoB_dimless*VyB;


  // We add rotation to the resulting flow
  if (_rotation == 1) {
    VphiB = VphiB +  2.*rBoundary_dimless*std::sin(thetaBoundary_dimless)*3.86e-3;
  }
  else if (_rotation == -1) {
    VphiB = VphiB - 2.*rBoundary_dimless*std::sin(thetaBoundary_dimless)*3.86e-3;
  }
  else {
    VphiB = VphiB;
  }

  // Finally, we transfer back to Cartesian and prescribe these values on the boundary surface
  VxB = xBoundary_dimless/rBoundary_dimless*VrB- yBoundary_dimless/rhoBoundary_dimless*VphiB+ xBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*VthetaB;
  VyB = yBoundary_dimless/rBoundary_dimless*VrB+ xBoundary_dimless/rhoBoundary_dimless*VphiB+ yBoundary_dimless*zBoundary_dimless/(rhoBoundary_dimless*rBoundary_dimless)*VthetaB;
  VzB = zBoundary_dimless/rBoundary_dimless*VrB- rhoBoundary_dimless/rBoundary_dimless*VthetaB;


  (*ghostState)[1] = 2.*VxB -(*innerState)[1]; 
  (*ghostState)[2] = 2.*VyB -(*innerState)[2]; 
  (*ghostState)[3] = 2.*VzB -(*innerState)[3]; 



  //===== D E N S I T Y   A N D   P R E S S U R E   B O U N D A R Y   C O N D I T I O N ===============


  // The way that the density and pressure values are prescribed depend on the nonhomogeneous flag.
  if (_nonhomogeneous) {
    // We first prescribe the homogeneous density value
    CFreal densityBoundary_dimless = _rhoBC;
    CFreal maxvA = _vAmax;
    // Then we compute the current Alfven speed
    CFreal vA = Bmag*Bref / pow((densityBoundary_dimless * rhoref * mu0),0.5);

    // If the current Alfven speed is smaller than the maximum set by the user, no adjustment is needed. Otherwise, the maximum set value is taken, with a double-sided hyperbolic profile in the middle.
    // This is achieved vie the following,
    smooth_factor = 0.0;
    dist = (vA - maxvA)/(0.1*maxvA);
    // where the smooth factor has a value of 0.0 if the Alfven speed is much smaller than the maximum set value and 1.0 if it exceeds it. The length over which the factor goes from 0 to 1 is defined by the distance above.
    smooth_factor = 0.5 + 0.5 * std::tanh(dist * 3.141592654);

    // The final density is then computed as a combination of the current value and the max value, with the weight set according to the smooth factor. 
    CFreal densityG = smooth_factor * (Bmag*Bref)*(Bmag*Bref) / (maxvA*maxvA) / mu0 + _rhoBC * rhoref * (1.0-smooth_factor); 
    (*ghostState)[0] = 2.0*(densityG / rhoref) - (*innerState)[0];



    // The same procedure is repeated for plasma beta and plasma pressure. The only difference here is that the user sets a minimum value, not maximum.
    CFreal minPlasmaBeta = _betamin; 
    CFreal pmag = (Bmag*Bref)*(Bmag*Bref)/2./mu0;
 
    CFreal PressureBoundary_dimless = _pBC; 
    CFreal pth = PressureBoundary_dimless * pref;

    CFreal plasmaBeta = pth/pmag;
  
    smooth_factor = 0.0;
    dist = (minPlasmaBeta - plasmaBeta)/(0.1*_betamin); 
    smooth_factor = 0.5 + 0.5 * std::tanh(dist * 3.141592654);  

    CFreal pthG = smooth_factor * pmag * minPlasmaBeta + pth * (1.0-smooth_factor);
    (*ghostState)[7] = 2.0*(pthG / pref) - (*innerState)[7]; }

  else {
     // Otherwise, in case of a homogeneous BC, everything is set according to the user uniformly. 
    (*ghostState)[0] = 2.0*_rhoBC - (*innerState)[0];
    (*ghostState)[7] = 2.0*_pBC - (*innerState)[7]; }




  //===== P S I    B O U N D A R Y   C O N D I T I O N ===============


  // In principle, the most physically appropriate boundary condition is (*ghostState)[8] = - (*innerState)[8]; which 
  // sets psi to zero at the boundary. In case of strong gradients near the inner boundary, however, this boundary 
  // condition may be too restrictive and lead to divergence. In that case, setting psi to 0 in the ghost cell can also
  // be done. This might result in odd speeds resolved right at the inner boundary, but otherwise the flowfield in the 
  // rest of the domain looks just like it would with (*ghostState)[8] = - (*innerState)[8] (tested). Do check the 
  // values of psi in the domain in the solution if the less restrictive BC is used, however. 
 

  (*ghostState)[8] = 0.; //-(*innerState)[8];
  

}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
