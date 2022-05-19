#include "Framework/MethodStrategyProvider.hh"

#include "MHD/MHD3DProjectionVarSet.hh"

#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "FluxReconstructionMHD/BCSuperOutletProj3DMHD.hh"

#include "Common/NotImplementedException.hh"

#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::MHD;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<
    BCSuperOutletProj3DMHD,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionMHDModule >
  BCSuperOutletProj3DMHDProvider("SuperOutletProj3DMHD");

//////////////////////////////////////////////////////////////////////////////

BCSuperOutletProj3DMHD::BCSuperOutletProj3DMHD(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData()
{
  CFAUTOTRACE;
  
  addConfigOptionsTo(this);

  m_refPhi = 0.0;
  setParameter("refPhi",&m_refPhi);
  
  m_imposeBDecay = false;
  setParameter("ImposeBDecay",&m_imposeBDecay);
}

//////////////////////////////////////////////////////////////////////////////

BCSuperOutletProj3DMHD::~BCSuperOutletProj3DMHD()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutletProj3DMHD::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("refPhi","Reference phi value imposed at the supersonic outlet.");
  options.addConfigOption< bool >("ImposeBDecay","Impose r2 decay of Br on outlet.");
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutletProj3DMHD::computeGhostStates(const vector< State* >& intStates,
                                                  vector< State* >& ghostStates,
                                                  const std::vector< RealVector >& normals,
                                                  const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    if (m_imposeBDecay) m_tempStates[iState] = (*intStates[iState]);
              
    // dereference states
    State& intState   = (*intStates  [iState]);
    State& ghostState = (*ghostStates[iState]);
    
    cf_assert(intState.size() == 9);
    cf_assert(ghostState.size() == 9);
        
    for (CFuint iEq = 0; iEq < 8; ++iEq)
    {
      ghostState[iEq] = intState[iEq];
    }
    
    //phi
    ghostState[8] = -intState[8];
    
    
//    const CFreal RSun = 6.9551e8; // m
//    const CFreal density_code = 1.67e-13;
//    const CFreal velocity_code = 2.2e-4/std::sqrt(1.256e-6*1.67e-13);
//    const CFreal pressure_code = std::pow(2.2e-4,2)/1.256e-6;
//    const CFreal B_code = 2.2e-4;
//
//    const CFreal xI = coords[iState][XX]*RSun;
//    const CFreal xI_dimless = coords[iState][XX];
//    const CFreal yI = coords[iState][YY]*RSun;
//    const CFreal yI_dimless = coords[iState][YY];
//    const CFreal zI = coords[iState][ZZ]*RSun;
//    const CFreal zI_dimless = coords[iState][ZZ];
//    const CFreal rI = (coords[iState]).norm2()*RSun;
//    const CFreal rI_dimless = (coords[iState]).norm2();
//
//    // to be checked if imposed minimum is needed (suspected problems at poles otherwise)
//    const CFreal rhoI = std::max(std::sqrt(xI*xI + yI*yI),1.0e-4*RSun);
//    const CFreal rhoI_dimless = std::max(std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless),1.0e-4);
//
//    const CFreal densityI_dimless = intState[0];
//    const CFreal densityI = densityI_dimless*density_code;
//  
//    const CFreal VxI_dimless = intState[1];
//    const CFreal VxI = VxI_dimless*velocity_code;
//    
//    const CFreal VyI_dimless = intState[2];
//    const CFreal VyI = VyI_dimless*velocity_code;
//  
//    const CFreal VzI_dimless = intState[3];
//    const CFreal VzI = VzI_dimless*velocity_code;
//  
//    const CFreal BxI_dimless = intState[4];
//    const CFreal BxI = BxI_dimless*B_code;
//  
//    const CFreal ByI_dimless = intState[5];
//    const CFreal ByI = ByI_dimless*B_code;
//  
//    const CFreal BzI_dimless = intState[6];
//    const CFreal BzI = BzI_dimless*B_code;
//  
//    const CFreal PrI_dimless = intState[7];
//    const CFreal PrI = PrI_dimless*pressure_code;
//
//    const CFreal BrI_dimless = xI_dimless/rI_dimless*BxI_dimless + yI_dimless/rI_dimless*ByI_dimless + zI_dimless/rI_dimless*BzI_dimless;
//    const CFreal BthetaI_dimless = xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*BxI_dimless + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*ByI_dimless - rhoI_dimless/rI_dimless*BzI_dimless;
//    const CFreal BphiI_dimless = -yI_dimless/rhoI_dimless*BxI_dimless + xI_dimless/rhoI_dimless*ByI_dimless;
//
//    const CFreal VrI_dimless = xI_dimless/rI_dimless*VxI_dimless + yI_dimless/rI_dimless*VyI_dimless + zI_dimless/rI_dimless*VzI_dimless;
//    const CFreal VthetaI_dimless = xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*VxI_dimless + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*VyI_dimless - rhoI_dimless/rI_dimless*VzI_dimless;
//    const CFreal VphiI_dimless = -yI_dimless/rhoI_dimless*VxI_dimless + xI_dimless/rhoI_dimless*VyI_dimless;
//
//    // Continuous BCs for r^2*Br, Btheta, Bphi:
//    const CFreal BrG_dimless = BrI_dimless;//rI_dimless*rI_dimless/(rG_dimless*rG_dimless)*BrI_dimless;
//    const CFreal BthetaG_dimless = BthetaI_dimless;
//    //Jens:
//    //CFreal BphiG_dimless = (rI_dimless/rG_dimless)*BphiI_dimless;
//    const CFreal BphiG_dimless = BphiI_dimless;
//  
//    const CFreal densityG_dimless = densityI_dimless;//rI_dimless*rI_dimless*densityI_dimless/(rG_dimless*rG_dimless);
//    
//    ghostState[0] = densityG_dimless;
//
//    // Continuous BCs for r^2*density*Vr, density*Vtheta, r*Vphi:
//    const CFreal VrG_dimless = VrI_dimless;//rI_dimless*rI_dimless*VrI_dimless*densityI_dimless/(rG_dimless*rG_dimless*densityG_dimless);
//    //CFreal VrG_dimless = rI_dimless*rI_dimless*VrI_dimless/(rG_dimless*rG_dimless);
//    const CFreal VthetaG_dimless = VthetaI_dimless;//densityI_dimless*VthetaI_dimless/densityG_dimless;
//    const CFreal VphiG_dimless = VphiI_dimless;//(rI_dimless/rG_dimless)*VphiI_dimless;
//
//    // UPDATE: BARBARA'S VELOCITY OUTLET CONDITIONS:
//    //CFreal VrG_dimless = VrI_dimless;
//    //CFreal VrG_dimless = rI_dimless*rI_dimless*VrI_dimless/(rG_dimless*rG_dimless);
//    //CFreal VthetaG_dimless = VthetaI_dimless;
//    //CFreal VphiG_dimless = VphiI_dimless;
//
//
//    // Convert all spherical vector components back to cartesian ones:
//    const CFreal BxG_dimless = xG_dimless/rG_dimless*BrG_dimless + xG_dimless*zG_dimless/(rhoG_dimless*rG_dimless)*BthetaG_dimless - yG_dimless/rhoG_dimless*BphiG_dimless;
//    const CFreal ByG_dimless = yG_dimless/rG_dimless*BrG_dimless + yG_dimless*zG_dimless/(rhoG_dimless*rG_dimless)*BthetaG_dimless + xG_dimless/rhoG_dimless*BphiG_dimless;
//    const CFreal BzG_dimless = zG_dimless/rG_dimless*BrG_dimless - rhoG_dimless/rG_dimless*BthetaG_dimless;
//    const CFreal VxG_dimless = xG_dimless/rG_dimless*VrG_dimless + xG_dimless*zG_dimless/(rhoG_dimless*rG_dimless)*VthetaG_dimless - yG_dimless/rhoG_dimless*VphiG_dimless;
//    const CFreal VyG_dimless = yG_dimless/rG_dimless*VrG_dimless + yG_dimless*zG_dimless/(rhoG_dimless*rG_dimless)*VthetaG_dimless + xG_dimless/rhoG_dimless*VphiG_dimless;
//    const CFreal VzG_dimless = zG_dimless/rG_dimless*VrG_dimless - rhoG_dimless/rG_dimless*VthetaG_dimless;
//
//    const CFreal BthetaB_dimless = 0.0;
//    const CFreal BphiB_dimless = 0.0;
//    const CFreal BrB_dimless = rI_dimless*rI_dimless/(rB_dimless*rB_dimless)*BrI_dimless;
//
//    const CFreal BxB_dimless = xB_dimless/rB_dimless*BrB_dimless + xB_dimless*zB_dimless/(rhoB_dimless*rB_dimless)*BthetaB_dimless - yB_dimless/rhoB_dimless*BphiB_dimless;
//    const CFreal ByB_dimless = yB_dimless/rB_dimless*BrB_dimless + yB_dimless*zB_dimless/(rhoB_dimless*rB_dimless)*BthetaB_dimless + xB_dimless/rhoB_dimless*BphiB_dimless;
//    const CFreal BzB_dimless = zB_dimless/rB_dimless*BrB_dimless - rhoB_dimless/rB_dimless*BthetaB_dimless;
//    
//    const CFreal VthetaB_dimless = 0.0;
//    const CFreal VphiB_dimless = 0.0;
//    const CFreal VrB_dimless = rI_dimless*rI_dimless/(rB_dimless*rB_dimless)*VrI_dimless;
//
//    const CFreal VxB_dimless = xB_dimless/rB_dimless*VrB_dimless + xB_dimless*zB_dimless/(rhoB_dimless*rB_dimless)*VthetaB_dimless - yB_dimless/rhoB_dimless*VphiB_dimless;
//    const CFreal VyB_dimless = yB_dimless/rB_dimless*VrB_dimless + yB_dimless*zB_dimless/(rhoB_dimless*rB_dimless)*VthetaB_dimless + xB_dimless/rhoB_dimless*VphiB_dimless;
//    const CFreal VzB_dimless = zB_dimless/rB_dimless*VrB_dimless - rhoB_dimless/rB_dimless*VthetaB_dimless;
//
//    ghostState[4] = BxG_dimless;
//    ghostState[5] = ByG_dimless;
//    ghostState[6] = BzG_dimless;
//
//    // density*r = continuous => densityG = densityI*rI/rG
//    // Continuous BC for r^2*rho: TYPO IN SKRALAN's PAPER
//    //CFreal densityG = rI*rI/(rG*rG)*densityI;
//
//    ghostState[1] = VxG_dimless;
//    ghostState[2] = VyG_dimless;
//    ghostState[3] = VzG_dimless;
//  
//    // Continuous temperature
//    ghostState[7] = PrI_dimless;//*rI_dimless*rI_dimless/(rG_dimless*rG_dimless);
//
//    ghostState[8] = -intState[8];
//
//     // there are two possible boundary conditions for phi
//
//     // 1) a reference value should be imposed for phi
//     //_dataGhostState[MHDProjectionTerm::PHI] = _refPhi;
//
//     // 2)
//     //_dataGhostState[MHDProjectionTerm::PHI] = -_dataInnerState[MHDProjectionTerm::PHI];    
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutletProj3DMHD::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                                     std::vector< std::vector< RealVector* > >& ghostGrads,
                                                     const std::vector< RealVector >& normals,
                                                     const std::vector< RealVector >& coords)
{
  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());
  
  const CFuint nbrGradVars = intGrads[0].size();
  
  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    
    if (m_imposeBDecay)
    {
      // normal
      const RealVector& normal = normals[iState];
      
      const CFreal xI_dimless = coords[iState][XX];
      const CFreal yI_dimless = coords[iState][YY];
      const CFreal zI_dimless = coords[iState][ZZ];
      const CFreal rI_dimless = (coords[iState]).norm2();
    
      // to be checked if imposed minimum is needed (suspected problems at poles otherwise)
      const CFreal rhoI_dimless = std::max(std::sqrt(xI_dimless*xI_dimless + yI_dimless*yI_dimless),1.0e-4);
      
      const CFreal densityI_dimless = m_tempStates[iState][0];
  
      const CFreal VxI_dimless = m_tempStates[iState][1];
      const CFreal VyI_dimless = m_tempStates[iState][2];
      const CFreal VzI_dimless = m_tempStates[iState][3];
  
      const CFreal BxI_dimless = m_tempStates[iState][4];
      const CFreal ByI_dimless = m_tempStates[iState][5];
      const CFreal BzI_dimless = m_tempStates[iState][6];
  
      const CFreal PrI_dimless = m_tempStates[iState][7];

      const CFreal BrI_dimless = xI_dimless/rI_dimless*BxI_dimless + yI_dimless/rI_dimless*ByI_dimless + zI_dimless/rI_dimless*BzI_dimless;
      const CFreal BthetaI_dimless = xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*BxI_dimless + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*ByI_dimless - rhoI_dimless/rI_dimless*BzI_dimless;
      const CFreal BphiI_dimless = -yI_dimless/rhoI_dimless*BxI_dimless + xI_dimless/rhoI_dimless*ByI_dimless;

      const CFreal VrI_dimless = xI_dimless/rI_dimless*VxI_dimless + yI_dimless/rI_dimless*VyI_dimless + zI_dimless/rI_dimless*VzI_dimless;
      const CFreal VthetaI_dimless = xI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*VxI_dimless + yI_dimless*zI_dimless/(rhoI_dimless*rI_dimless)*VyI_dimless - rhoI_dimless/rI_dimless*VzI_dimless;
      const CFreal VphiI_dimless = -yI_dimless/rhoI_dimless*VxI_dimless + xI_dimless/rhoI_dimless*VyI_dimless;
      
      //phi
      *ghostGrads[iState][8] = *intGrads[iState][8]; 
      
      // density
      const RealVector& densGradI =  *intGrads[iState][0];
      RealVector& densGradG =  *ghostGrads[iState][0];
      const CFreal nDensGrad = MathFunctions::innerProd(densGradI, normal);
      
      const CFreal ddensdrG = -2.0*densityI_dimless/rI_dimless;
      
      densGradG = densGradI + 2.0*(ddensdrG - nDensGrad)*normal;
      
      // B
      const CFreal dBrdrG = -2.0*BrI_dimless/rI_dimless;
      //const CFreal dBthetadrG = 0.0;
      //const CFreal dBphidrG = 0.0;
      
      // Bx
      const RealVector& BxGradI =  *intGrads[iState][4];
      RealVector& BxGradG =  *ghostGrads[iState][4];
      const CFreal nBxGrad = MathFunctions::innerProd(BxGradI, normal);
      
      const CFreal dBxdrG = xI_dimless/rI_dimless*dBrdrG;
      
      BxGradG = BxGradI + 2.0*(dBxdrG - nBxGrad)*normal;
      
      // By
      const RealVector& ByGradI =  *intGrads[iState][5];
      RealVector& ByGradG =  *ghostGrads[iState][5];
      const CFreal nByGrad = MathFunctions::innerProd(ByGradI, normal);
      
      const CFreal dBydrG = yI_dimless/rI_dimless*dBrdrG;
      
      ByGradG = ByGradI + 2.0*(dBydrG - nByGrad)*normal;
      
      // Bz
      const RealVector& BzGradI =  *intGrads[iState][6];
      RealVector& BzGradG =  *ghostGrads[iState][6];
      const CFreal nBzGrad = MathFunctions::innerProd(BzGradI, normal);
      
      const CFreal dBzdrG = zI_dimless/rI_dimless*dBrdrG;
      
      BzGradG = BzGradI + 2.0*(dBzdrG - nBzGrad)*normal;
      
      // V
      //const CFreal dVrdrG = 0.0;
      const CFreal dVthetadrG = 2.0*VthetaI_dimless/rI_dimless;
      const CFreal dVphidrG = -VphiI_dimless/rI_dimless;
      
      // Vx
      const RealVector& VxGradI =  *intGrads[iState][1];
      RealVector& VxGradG =  *ghostGrads[iState][1];
      const CFreal nVxGrad = MathFunctions::innerProd(VxGradI, normal);
      
      const CFreal dVxdrG = xI_dimless*zI_dimless/(rI_dimless*rhoI_dimless)*dVthetadrG - yI_dimless/rhoI_dimless*dVphidrG;
      
      VxGradG = VxGradI + 2.0*(dVxdrG - nVxGrad)*normal;
      
      // Vy
      const RealVector& VyGradI =  *intGrads[iState][2];
      RealVector& VyGradG =  *ghostGrads[iState][2];
      const CFreal nVyGrad = MathFunctions::innerProd(VyGradI, normal);
      
      const CFreal dVydrG = yI_dimless*zI_dimless/(rI_dimless*rhoI_dimless)*dVthetadrG + xI_dimless/rhoI_dimless*dVphidrG;
      
      VyGradG = VyGradI + 2.0*(dVydrG - nVyGrad)*normal;
      
      // Vz
      const RealVector& VzGradI =  *intGrads[iState][3];
      RealVector& VzGradG =  *ghostGrads[iState][3];
      const CFreal nVzGrad = MathFunctions::innerProd(VzGradI, normal);
      
      const CFreal dVzdrG = -rhoI_dimless/rI_dimless*dVthetadrG;
      
      VzGradG = VzGradI + 2.0*(dVzdrG - nVzGrad)*normal;
      
      // Pr
      const RealVector& PrGradI =  *intGrads[iState][7];
      RealVector& PrGradG =  *ghostGrads[iState][7];
      const CFreal nPrGrad = MathFunctions::innerProd(PrGradI, normal);
      
      const CFreal dPrdrG = -2.0*PrI_dimless/rI_dimless;
      
      PrGradG = PrGradI + 2.0*(dPrdrG - nPrGrad)*normal;
      
    }
    else
    {
      for (CFuint iGradVar = 0; iGradVar < nbrGradVars; ++iGradVar)
      {
        *ghostGrads[iState][iGradVar] = *intGrads[iState][iGradVar]; 
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCSuperOutletProj3DMHD::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = true;

//  // get MHD 3D varset
//  m_varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
//  if (m_varSet.isNull())
//  {
//    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MHD3DProjectionVarSet in BCSuperOutletProj3DMHD!");
//  }
//
//  // resize the physical data for internal and ghost solution points
//  m_varSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
//  m_varSet->getModel()->resizePhysicalData(m_intSolPhysData  );
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  
  // compute flux point coordinates
  SafePtr< vector<RealVector> > flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  const CFreal nbrFaceFlxPnts = flxLocalCoords->size();
  
  m_tempStates.resize(nbrFaceFlxPnts);
  
  // number of equations
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  for (CFuint iFlx = 0; iFlx < nbrFaceFlxPnts; ++iFlx)
  {
    m_tempStates[iFlx].resize(nbEqs); 
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

