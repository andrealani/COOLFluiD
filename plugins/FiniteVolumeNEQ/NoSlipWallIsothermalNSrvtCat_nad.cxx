#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "FiniteVolumeNEQ/NoSlipWallIsothermalNSrvtCat_nad.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

#include "MathTools/MatrixInverter.hh"

#include "NavierStokes/MultiScalarVarSet.hh"
#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/Euler2DVarSet.hh"

#include <cmath>
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::MathTools;
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NoSlipWallIsothermalNSrvtCat_nad,
                      CellCenterFVMData,
                       FiniteVolumeNEQModule>
NoSlipWallIsothermalNSrvtCat_nadFVMCCProvider("NoSlipWallIsothermalNSrvtCat_nadFVMCC");

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSrvtCat_nad::defineConfigOptions(Config::OptionList& options)
{
 options.addConfigOption< CFreal > ("Tw","Temperature at the wall");
 options.addConfigOption< CFreal > ("Nr","Number of catalised equations at the wall");
 options.addConfigOption< CFuint > ("NewtonLoop","Number of Newton loop that we want to do");


}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallIsothermalNSrvtCat_nad::NoSlipWallIsothermalNSrvtCat_nad(const std::string& name) :
  FVMCC_BC(name),
  socket_walldistance("wallDistance")
{
  addConfigOptionsTo(this);
  
  m_Tw = 0.0;
  setParameter("Tw",&m_Tw);
  
  m_nr = 0.0;
  setParameter("Nr",&m_nr);
  
  m_nl = 2;
  setParameter("NewtonLoop",&m_nl);

   m_coeff = 2.0;
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallIsothermalNSrvtCat_nad::~NoSlipWallIsothermalNSrvtCat_nad()
{
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSrvtCat_nad::setGhostState(GeometricEntity *const face)
{
  // Get the inner and boundary states
  State& innerState = *face->getState(0);
  State& ghostState = *face->getState(1);
  CFreal  RUNIV = 8.314511212;
  SafePtr<PhysicalChemicalLibrary> library =
  PhysicalModelStack::getActive()->getImplementor()->
  getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  assert(library.isNotNull());

  const bool initialization = getMethodData().isInitializationPhase();
    DataHandle<CFreal> WallStorage = socket_walldistance.getDataHandle();
    m_walldistance = WallStorage[innerState.getLocalID()];

  // When we initialize the wall distance is to yet available,
  //so we just extrapolate the ghost states
  if (initialization || m_walldistance == 0.0){
    for (CFuint i = 0; i < m_nbSpecies; ++i) {
      ghostState[i] = innerState[i];
      
    }
    ghostState[m_nbSpecies] = 0.0;
    ghostState[m_nbSpecies+1] = 0.0;
    ghostState[m_nbSpecies+2] = m_Tw;
  }
  else{
    // Get the wall distance
    DataHandle<CFreal> WallStorage = socket_walldistance.getDataHandle();
    m_walldistance = WallStorage[innerState.getLocalID()];

    // Compute the molar fraction: We first compute the massic fraction
    // and then we convert to molar fraction
    compute_mass_fraction(innerState, ghostState, m_yb, m_yi, m_rhob, m_rhoi);
    mass_frac_to_mol_frac(m_yb, m_yi, m_xb, m_xi);

    // This comand is to fill Y and X of mutation with the values of the gostcell
    library->setSpeciesFractions(m_yb);
     
    CFreal m_p = library->pressure(m_rhoi, m_Tw, CFNULL);
  
    // Newton loop (few loops are necessary)
    for (CFuint k = 0; k < m_nl; k++)
    { 
      
      for (CFuint i = 0; i < m_nbSpecies; ++i) 
      {
        m_dx[i] = m_xb[i] - m_xi[i];
      }
      m_dx /= m_walldistance;
      
      library->getDij_fick(m_dx, m_p, m_Tw, m_Dbin, m_Diff_flux);
      CFreal mmt;
      getMolarMass(m_xb,mmt);
      CFreal m_rho=m_p/(RUNIV/mmt*m_Tw);
      getWallRateProd(m_rho, m_yb, m_Tw, m_nu, m_muN, m_muO, m_muN2, m_muNO, m_muO2, m_omega);
      
      // Compute the right hand side that is diff_flux - m_omega
      for (CFuint i = 0; i < m_nbSpecies; ++i) 
      {
        m_b[i] = m_Diff_flux[i] - m_omega[i];
	
      }
      
      CFreal eps=1.0e-7;  
      // We compute the derivative of the wall production and diffusive flux by numerical perturbation
      m_xp = m_xb;
      m_dxp = m_dx;
      for (CFuint i = 0; i < m_nbSpecies; ++i) 
        { 
          for (CFuint j = 0; j < m_nbSpecies; ++j) 
          { 
            m_xp[j] = m_xb[j] + eps;
	    m_dxp[j] = (m_xp[j] - m_xi[i])/ m_walldistance;
            
	    library->getSpeciesMassFractions(m_xp, m_yp);
	    library->setSpeciesFractions(m_yp);
	    getMolarMass(m_xp,mmt);
	    m_rho=m_p/(RUNIV/mmt*m_Tw);
	    
            getWallRateProd(m_rho, m_yp, m_Tw, m_nu, m_muN, m_muO, m_muN2, m_muNO, m_muO2, m_omegap);
            library->getDij_fick(m_dxp, m_p, m_Tw, m_Dbin, m_Diff_fluxp);
	    
	    m_a(i,j)= ((m_Diff_fluxp[i]-m_Diff_flux[i]) -(m_omegap[i] - m_omega[i]))/eps;
            m_xp[j] = m_xb[j];
	    m_dxp[j] = m_dx[j];
            }
        } 
    m_inverter->invert(m_a, m_inva);
    m_b = m_inva*m_b; 
    m_xb -= m_b;
    
    
    library->getSpeciesMassFractions(m_xb, m_yb);
    library->setSpeciesFractions(m_yb);
     
    }
     // We go back to the partial densities
    mol_frac_to_part_dens(m_xb, m_p, m_Tw, m_partrho);

    for (CFuint i = 0; i < m_nbSpecies; ++i) 
      {
        ghostState[i] = max(0.0,(2.0*m_partrho[i] - innerState[i]));
      }
      
  ghostState[m_nbSpecies] = -innerState[m_nbSpecies];
  ghostState[m_nbSpecies+1] = -innerState[m_nbSpecies+1];
  ghostState[m_nbSpecies+2] = max(100.0,2.0*m_Tw-innerState[m_nbSpecies+2]);
  }
}


//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSrvtCat_nad::setup()
 {
   FVMCC_BC::setup();
  _varSet = getMethodData().getUpdateVar().d_castTo<MultiScalarVarSet<Euler2DVarSet> >();
   m_nbSpecies = _varSet->getModel()->getNbScalarVars(0); 
   m_yi.resize(m_nbSpecies);
   m_xi.resize(m_nbSpecies);
   m_yb.resize(m_nbSpecies);
   m_yp.resize(m_nbSpecies);
   m_xb.resize(m_nbSpecies);
   m_xp.resize(m_nbSpecies);
   m_mm.resize(m_nbSpecies);
   m_dx.resize(m_nbSpecies);
   m_dxp.resize(m_nbSpecies);
   m_omega.resize(m_nbSpecies);
   m_omegap.resize(m_nbSpecies);
   m_Dbin.resize(m_nbSpecies,m_nbSpecies);
   m_Diff_flux.resize(m_nbSpecies);
   m_Diff_fluxp.resize(m_nbSpecies);
   m_b.resize(m_nbSpecies);
   m_zero.resize(m_nbSpecies);
   m_zero = 0.0;
   m_mcal.resize(m_nbSpecies);
   m_nu.resize(m_nbSpecies,m_nr);
   m_muN.resize(m_nbSpecies,m_nr);
   m_muO.resize(m_nbSpecies,m_nr);
   m_muN2.resize(m_nbSpecies,m_nr);
   m_muNO.resize(m_nbSpecies,m_nr);
   m_muO2.resize(m_nbSpecies,m_nr);
   m_a.resize(m_nbSpecies,m_nbSpecies);
   m_inva.resize(m_nbSpecies,m_nbSpecies);
   m_partrho.resize(m_nbSpecies);
   m_rhoG.resize(m_nbSpecies);

   m_nu(0,0) = 1.; m_nu(1,0) = 0.; m_nu(2,0) = 0.; m_nu(3,0) = 0.; m_nu(4,0) = 0.;
   m_nu(0,1) = 0.; m_nu(1,1) = 1.; m_nu(2,1) = 0.; m_nu(3,1) = 0.; m_nu(4,1) = 0.;


   m_muN(0,0) = 0.; m_muN(1,0) = 0.; m_muN(2,0) = 1.; m_muN(3,0) = 0.; m_muN(4,0) = 0.;
   m_muN(0,1) = 0.; m_muN(1,1) = 0.; m_muN(2,1) = 0.; m_muN(3,1) = 0.; m_muN(4,1) = 0.; 
                      
                      
   m_muO(0,0) = 0.; m_muO(1,0) = 0.; m_muO(2,0) = 0.; m_muO(3,0) = 0.; m_muO(4,0) = 0.;
   m_muO(0,1) = 0.; m_muO(1,1) = 0.; m_muO(2,1) = 0.; m_muO(3,1) = 0.; m_muO(4,1) = 1.; 


   m_muN2(0,0) = 0.; m_muN2(1,0) = 0.; m_muN2(2,0) = 0.; m_muN2(3,0) = 0.; m_muN2(4,0) = 0.;
   m_muN2(0,1) = 0.; m_muN2(1,1) = 0.; m_muN2(2,1) = 0.; m_muN2(3,1) = 0.; m_muN2(4,1) = 0.; 
                                                                             
                                                                             
   m_muNO(0,0) = 0.; m_muNO(1,0) = 0.; m_muNO(2,0) = 0.; m_muNO(3,0) = 0.; m_muNO(4,0) = 0.;
   m_muNO(0,1) = 0.; m_muNO(1,1) = 0.; m_muNO(2,1) = 0. ; m_muNO(3,1) = 0.; m_muNO(4,1) = 0.; 

   m_muO2(0,0) = 0.; m_muO2(1,0) = 0.; m_muO2(2,0) = 0.; m_muO2(3,0) = 0.; m_muO2(4,0) = 0.;
   m_muO2(0,1) = 0.; m_muO2(1,1) = 0.; m_muO2(2,1) = 0. ; m_muO2(3,1) = 0.; m_muO2(4,1) = 0.; 
        
  m_inverter.reset(MatrixInverter::create(m_nbSpecies, false));

  
  // the temperature ID is equal to the maximum velocity ID + 1
 
  
  m_tempNode.resize(PhysicalModelStack::getActive()->getDim());
  m_midNode.resize(PhysicalModelStack::getActive()->getDim());
  m_tempGhostNode.resize(PhysicalModelStack::getActive()->getDim());
  m_faceNormal.resize(PhysicalModelStack::getActive()->getDim());

                                
}                                     


//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
NoSlipWallIsothermalNSrvtCat_nad::needsSockets()
{
   std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = FVMCC_BC::needsSockets();

   result.push_back(&socket_walldistance);

   return result;
}
//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSrvtCat_nad::compute_mass_fraction(State& innerState, State& ghostState, RealVector& m_yb, RealVector& m_yi, CFreal& m_rhob, CFreal& m_rhoi)
{
  m_rhoi = 0.0;
  m_rhob = 0.0;
  for (CFuint i = 0; i < m_nbSpecies; ++i) {
    m_rhoi += innerState[i];
    m_rhob += ghostState[i];
  }

  for (CFuint i = 0; i < m_nbSpecies; ++i) {
    m_yi[i] = innerState[i]/m_rhoi;
    m_yb[i] = ghostState[i]/m_rhob;
  }

}

//////////////////////////////////////////////////////////////////////////////////

  void NoSlipWallIsothermalNSrvtCat_nad::mass_frac_to_mol_frac(RealVector& m_yb, RealVector& m_yi, RealVector& m_xb, RealVector& m_xi){

  SafePtr<PhysicalChemicalLibrary> library =
  PhysicalModelStack::getActive()->getImplementor()->
  getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  assert(library.isNotNull());

  // First we compute the molar mass of the gas
  library-> getMolarMasses(m_mm);
  CFreal m_gasmmi = 0.;
  CFreal m_gasmmb = 0.;
  CFreal tmpi = 0.0;
  CFreal tmpb = 0.0;

  for (CFuint i = 0; i < m_nbSpecies; ++i) {
    tmpi += m_yi[i]/m_mm[i];
    tmpb += m_yb[i]/m_mm[i];    
  }

  m_gasmmi = 1.0/tmpi;
  m_gasmmb = 1.0/tmpb;

//   CF_DEBUG_OBJ(m_gasmmi);
//   CF_DEBUG_OBJ(m_gasmmb);
  /*  CF_DEBUG_OBJ(m_mm);*/

//   CF_DEBUG_OBJ(m_yb);
  
  for (CFuint i = 0; i < m_nbSpecies; ++i) {
    m_xi[i] = m_yi[i]/m_mm[i];
    m_xb[i] = m_yb[i]/m_mm[i];
  }

  m_xi *= m_gasmmi;
  m_xb *=  m_gasmmb;
// CF_DEBUG_OBJ(m_xb);
}
//////////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSrvtCat_nad::mol_frac_to_part_dens(RealVector& m_xb, CFreal& m_p, CFreal& m_Tw, RealVector& m_partrho){


  SafePtr<PhysicalChemicalLibrary> library =
  PhysicalModelStack::getActive()->getImplementor()->
  getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  assert(library.isNotNull());
  library-> getMolarMasses(m_mm);
  CFreal RUNIV = 8.31451;
  CFreal m_gasmm = 0.;
  CFreal m_rho = 0.;
  
 for (CFuint i = 0; i < m_nbSpecies; ++i) {
    m_gasmm += m_xb[i]*m_mm[i];
 
  }

  CFreal m_Rgb = RUNIV/m_gasmm;
  m_rho = m_p/(m_Rgb*m_Tw);
  for (CFuint i = 0; i < m_nbSpecies; ++i) {
    m_partrho[i] = (m_xb[i]*m_mm[i]/m_gasmm);
  }
  m_partrho *= m_rho;  
}

//////////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSrvtCat_nad::mol_frac_to_mass_frac(RealVector& m_xb, RealVector& m_yb){


  SafePtr<PhysicalChemicalLibrary> library =
  PhysicalModelStack::getActive()->getImplementor()->
  getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  assert(library.isNotNull());
  library-> getMolarMasses(m_mm);
  //  CFreal RUNIV = 8.31451;
  CFreal m_gasmm = 0.;
  // CFreal m_rho = 0.;
  
 for (CFuint i = 0; i < m_nbSpecies; ++i) {
    m_gasmm += m_xb[i]*m_mm[i];
 
  }

  for (CFuint i = 0; i < m_nbSpecies; ++i) {
    m_yb[i] = (m_xb[i]*m_mm[i]/m_gasmm);
  }
}
//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSrvtCat_nad::getWallRateProd(CFreal& m_rhob,
                       RealVector& m_yb,
                       CFreal& m_Tw,
                       RealMatrix& m_nu,
                       RealMatrix& m_muN,
                       RealMatrix& m_muO,
                       RealMatrix& m_muN2,
                       RealMatrix& m_muNO,
                       RealMatrix& m_muO2,
                       RealVector& m_omegawall)
{                      
                       
  const CFreal Runiv = 8.314511212;
  const CFreal Pi = 3.14159265;
  CFreal m_tmp;
  CFreal m_tmp2;

  SafePtr<PhysicalChemicalLibrary> library =
  PhysicalModelStack::getActive()->getImplementor()->
  getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  assert(library.isNotNull());
  library-> getMolarMasses(m_mm);
  library-> getGammaO(m_GammaO);
  library-> getGammaN(m_GammaN);

//  cout << m_GammaO << endl;
//  cout << m_GammaN << endl;
 
  for (CFuint i= 0; i<m_nbSpecies; i++)
  {
    // We first compute the impinging flux times mi (for derivation see with NV)
    m_mcal[i] = m_yb[i]*m_rhob*sqrt(m_Tw*Runiv/(2.0*Pi*m_mm[i]));
   // CF_DEBUG_OBJ( m_yb[i]);
   // CF_DEBUG_OBJ(m_Tw);
   // CF_DEBUG_OBJ(m_rhob);
 //   CF_DEBUG_OBJ(m_mcal[i]);
  }
  
 /* CF_DEBUG_OBJ(m_mcal[0]);
  CF_DEBUG_OBJ(m_mcal[1]);*//*
  CF_DEBUG_OBJ(m_mcal[2]);
  CF_DEBUG_OBJ(m_mcal[3]);
  CF_DEBUG_OBJ(m_mcal[4]);*/
  
  RealVector m_gammaV;
  m_gammaV.resize(m_nbSpecies);
  m_gammaV[0] = m_GammaN;
  m_gammaV[1] = m_GammaO;
  m_gammaV[2] = 0.;
  m_gammaV[3] = 0.;
  m_gammaV[4] = 0.;

  
  for (CFuint i= 0; i<m_nbSpecies; i++)
  {
    // Now we compute the destructing part and the [production part
    m_tmp = 0.0;
    m_tmp2 = 0.0;
    for (CFuint r = 0; r < m_nr; r++)
    {
   
      m_tmp += m_gammaV[i]*m_nu(i,r);
      m_tmp2 += (m_GammaN*m_muN(i,r)*m_mcal[0]  + m_GammaO*m_muO(i,r)*m_mcal[1] );
  
    }
    m_tmp *= m_mcal[i];
    
//     CF_DEBUG_OBJ(m_tmp);
//     CF_DEBUG_OBJ(m_tmp2);
    
    m_omegawall[i] = m_tmp - m_tmp2; 
   
  }
  //  CF_DEBUG_OBJ(m_omegawall[0]);
  //  CF_DEBUG_OBJ(m_omegawall[1]);
  //  CF_DEBUG_OBJ(m_omegawall[2]);
  //  CF_DEBUG_OBJ(m_omegawall[3]);
  //  CF_DEBUG_OBJ(m_omegawall[4]);
}
//////////////////////////////////////////////////////////////////////////////

void NoSlipWallIsothermalNSrvtCat_nad::getMolarMass(RealVector &m_xp,CFreal &m_mmt){
SafePtr<PhysicalChemicalLibrary> library =
  PhysicalModelStack::getActive()->getImplementor()->
  getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  assert(library.isNotNull());
  library-> getMolarMasses(m_mm);
  m_mmt = 0.;
  
 for (CFuint i = 0; i < m_nbSpecies; ++i) {
    m_mmt += m_xp[i]*m_mm[i];
    
 
  }
}
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
