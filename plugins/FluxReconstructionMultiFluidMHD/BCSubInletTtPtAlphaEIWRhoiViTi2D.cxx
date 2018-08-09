#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMultiFluidMHD/BCSubInletTtPtAlphaEIWRhoiViTi2D.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Framework/PhysicalConsts.hh"
#include "Framework/MethodStrategyProvider.hh"
#include <iostream>
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::MultiFluidMHD;
using namespace COOLFluiD::Physics::Maxwell;
using namespace COOLFluiD::Config;


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<BCSubInletTtPtAlphaEIWRhoiViTi2D, FluxReconstructionSolverData, BCStateComputer, FluxReconstructionMultiFluidMHDModule> subInletTtPtAlphaEIWRhoiViTiProvider("BCSubInletTtPtAlphaEIWRhoiViTi2D");

//////////////////////////////////////////////////////////////////////////////

void BCSubInletTtPtAlphaEIWRhoiViTi2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<CFreal> >("Ttot","total temperature");
  options.addConfigOption< std::vector<CFreal> >("Ptot","total pressure");
  options.addConfigOption< std::vector<CFreal> >("alpha","alpha");
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletTtPtAlphaEIWRhoiViTi2D::BCSubInletTtPtAlphaEIWRhoiViTi2D(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData(),
  m_tTotal(),
  m_pTotal(),
  m_alpha()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  m_tTotal = std::vector<CFreal>();
   setParameter("Ttot",&m_tTotal);

  m_pTotal = std::vector<CFreal>();
   setParameter("Ptot",&m_pTotal);

  m_alpha = std::vector<CFreal>();
   setParameter("alpha",&m_alpha);
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletTtPtAlphaEIWRhoiViTi2D::~BCSubInletTtPtAlphaEIWRhoiViTi2D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletTtPtAlphaEIWRhoiViTi2D::computeGhostStates(const vector< State* >& intStates,
                                                    vector< State* >& ghostStates,
                                                    const std::vector< RealVector >& normals,
                                                    const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

  CFLog(VERBOSE, "\n\n\n Error 1 \n");

/*
  cf_assert(m_tTotal.size() == nbrStates);
  cf_assert(m_pTotal.size() == nbrStates);
*/
  CFLog(VERBOSE, "\n\n\n Error 2 \n");

  // get some physical data from the model
  const CFreal gamma = m_varSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);

  CFLog(VERBOSE, "\n\n\n\n\n gamma = " << gamma << "\n\n\n\n\n");


  ///MultiFluidMHD Subsonic Outlet Condition imposing pressure in 2D
  // here a fix is needed in order to have always ghostT > 0
  // if ghostT < 0  then the inner value is set
  const CFuint endEM = 8;

 // **********************TEMPORARY CHANGE********************************* 

  const CFuint firstSpecies = m_varSet->getModel()->getFirstScalarVar(0);
  const CFuint firstVelocity = m_varSet->getModel()->getFirstScalarVar(1);
  const CFuint firstTemperature = m_varSet->getModel()->getFirstScalarVar(2);
  

  CFLog(VERBOSE, "\n\n\n\n\n firstSpecies Size = " << firstSpecies << "\n\n\n\n\n");
  CFLog(VERBOSE, "\n\n\n\n\n firstVelocity Size = " << firstVelocity << "\n\n\n\n\n");
  CFLog(VERBOSE, "\n\n\n\n\n First Temperature Size = " << firstTemperature << "\n\n\n\n\n");
  CFLog(VERBOSE, "\n\n\n\n\n  m_ghostSolPhysData Size = " << m_ghostSolPhysData.size() << "\n\n\n\n\n");
  CFLog(VERBOSE, "\n\n\n\n\n  m_intSolPhysData Size = " << m_intSolPhysData.size() << "\n\n\n\n\n");
  
  //const CFuint tgAlpha = tan(m_alpha);
  

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
      // dereference states
      State& intState   = (*intStates[iState]);
      State& ghostState = (*ghostStates[iState]);

      //cf_assert(intState.size() == 4);
      //cf_assert(ghostState.size() == 4);

      // set the physical data starting from the inner state
      m_varSet->computePhysicalData(intState,m_intSolPhysData);


  CFLog(VERBOSE, "\n\n\n Error 1 \n");

    (m_ghostSolPhysData)[0] = (m_intSolPhysData)[0] /*+ 2*bn*nx*/;	//Bx
    (m_ghostSolPhysData)[1] = (m_intSolPhysData)[1] /*+ 2*bn*ny*/;	//By
    (m_ghostSolPhysData)[2] = (m_intSolPhysData)[2] /*+ 2*bn*nz*/;	//Bz
    (m_ghostSolPhysData)[3] = (m_intSolPhysData)[3] /*- 2*en*nx*/;	//Ex
    (m_ghostSolPhysData)[4] = (m_intSolPhysData)[4] /*- 2*en*ny*/;	//Ey
    (m_ghostSolPhysData)[5] = (m_intSolPhysData)[5] /*- 2*en*nz*/;	//Ez
    (m_ghostSolPhysData)[6] = (m_intSolPhysData)[6];			//Psi
    (m_ghostSolPhysData)[7] = (m_intSolPhysData)[7];			//Phi

    const bool isLeake = true;

  if (isLeake) {

    const CFreal m_electron = m_varSet->getModel()->getMolecularMass1();
    const CFreal m_neutral = m_varSet->getModel()->getMolecularMass2();
    const CFreal m_ion = m_varSet->getModel()->getMolecularMass3();
    const CFreal K_B = PhysicalConsts::Boltzmann();
    const CFreal Rplasma = 2*K_B/m_ion;				// ions gas constant
    const CFreal Rneutral = K_B/m_neutral;

  CFLog(VERBOSE, "\n\n\n Error 2 \n");

    const CFreal uxi = (m_intSolPhysData)[firstVelocity];
    const CFreal uyi = (m_intSolPhysData)[firstVelocity + 1];
    const CFreal uxn = (m_intSolPhysData)[firstVelocity + 2];
    const CFreal uyn = (m_intSolPhysData)[firstVelocity + 3];

  CFLog(VERBOSE, "\n\n\n Error 3 \n");

    const CFreal Ui = std::sqrt(uxi*uxi + uyi*uyi);
    const CFreal Un = std::sqrt(uxn*uxn + uyn*uyn);

  CFLog(VERBOSE, "\n\n\n Error 3 \n");

    const CFreal rhoi = (m_intSolPhysData)[endEM]*(m_intSolPhysData)[firstSpecies];
    const CFreal rhon = (m_intSolPhysData)[endEM]*(m_intSolPhysData)[firstSpecies + 1];

  CFLog(VERBOSE, "\n\n\n Error 3 \n");

    const CFreal Ti = (m_intSolPhysData)[firstTemperature];	 //CFLog(VERBOSE, "\n\n\n\n Temperature i= " << Ti << "\n\n\n);
    const CFreal Tn = (m_intSolPhysData)[firstTemperature + 1];  //(m_intSolPhysData)[firstTemperature + 4];
    const CFreal pi = rhoi*Rplasma*Ti;
    const CFreal pn = rhon*Rneutral*Tn;

  CFLog(VERBOSE, "\n\n\n Error 3 \n");

    const CFreal ai = (m_intSolPhysData)[firstTemperature + 3];  //(m_intSolPhysData)[firstTemperature + 2];
    const CFreal an = (m_intSolPhysData)[firstTemperature + 4];  //(m_intSolPhysData)[firstTemperature + 6];
    const CFreal Mi = Ui/ai;
    const CFreal Mn = Un/an;

  CFLog(VERBOSE, "\n\n\n\n\n Mi = " << Mi << "\n\n\n\n\n");
  CFLog(VERBOSE, "\n\n\n Error 3 \n");

    const CFreal coeffi = 1 + 0.5*(gamma - 1)*Mi*Mi;
    const CFreal coeffn = 1 + 0.5*(gamma - 1)*Mn*Mn;


  CFLog(VERBOSE, "\n\n\n\n\n coeffi = " << coeffi << "\n\n\n\n\n");

  CFLog(VERBOSE, "\n\n\n Error 4 \n");


    const CFreal coeffPow_i = pow(coeffi, gamma/(gamma - 1));
    const CFreal coeffPow_n = pow(coeffn, gamma/(gamma - 1));


  CFLog(VERBOSE, "\n\n\n\n\n coeffPow_i = " << coeffPow_i << "\n\n\n\n\n");



  CFLog(VERBOSE, "\n\n\n Error 6 \n");

    const CFreal Ttoti = Ti*coeffi;
    const CFreal Ttotn = Tn*coeffn;
    //const CFreal coeff2 = (1 + gamma*Mi*Mi);
    const CFreal Ptoti = pi*coeffPow_i;
    const CFreal Ptotn = pn*coeffPow_n;

  CFLog(VERBOSE, "\n\n\n Error 7 \n");

    CFreal tgAlphai = 0.;
    CFreal tgAlphan = 0.;
    if(uxi != 0.){tgAlphai = uyi/uxi;}
    if(uxn != 0.){tgAlphan = uyn/uxn;}


    CFLog(VERBOSE, "\n\n\n Here 2 \n");

    const CFreal Ttotg_i = 2*m_tTotal[0] - Ttoti;

    CFLog(VERBOSE, "\n\n\n Here 3 \n");
    const CFreal Ptotg_i = 2*m_pTotal[0] - Ptoti;



    const CFreal Ttotg_n = 2*m_tTotal[1] - Ttotn;
    const CFreal Ptotg_n = 2*m_pTotal[1] - Ptotn;

    CFLog(VERBOSE, "\n\n\n Here 4 \n");

    const CFreal tgAlphag_i = 2*tan(m_alpha[0]) - tgAlphai;
    const CFreal tgAlphag_n = 2*tan(m_alpha[1]) - tgAlphan;

    CFLog(VERBOSE, "\n\n\n Here 5 \n");

    const CFreal Mg_i = Mi;
    const CFreal Mg_n = Mn;

    CFLog(VERBOSE, "\n\n\n Here 6 \n");

    const CFreal Tg_i = Ttotg_i/coeffi;

    CFLog(VERBOSE, "\n\n\n Here 7 \n");

    const CFreal pg_i = Ptotg_i/coeffPow_i;

    CFLog(VERBOSE, "\n\n\n Here 8 \n");
    const CFreal Tg_n = Ttotg_n/coeffn;
    CFLog(VERBOSE, "\n\n\n Here 9 \n");
    const CFreal pg_n = Ptotg_n/coeffPow_n;
    CFLog(VERBOSE, "\n\n\n Here 10 \n");
    const CFreal rhog_i = pg_i/(Rplasma*Tg_i);

    CFLog(VERBOSE, "\n\n\n Here 11 \n");
    const CFreal rhog_n = pg_n/(Rneutral*Tg_n);
    const CFreal rhog_total = rhog_i + rhog_n;

    CFLog(VERBOSE, "\n\n\n Here 12 \n");

    const CFreal uxg_i = Mg_i*std::sqrt(gamma*Rplasma*Tg_i/(1 + tgAlphag_i*tgAlphag_i));
    const CFreal uyg_i = tgAlphag_i*uxg_i;

    CFLog(VERBOSE, "\n\n\n Here 13 \n");
    const CFreal uxg_n = Mg_n*std::sqrt(gamma*Rneutral*Tg_n/(1 + tgAlphag_n*tgAlphag_n));
    const CFreal uyg_n = tgAlphag_n*uxg_n;

    const CFreal UgUg_i = uxg_i*uxg_i + uyg_i*uyg_i;
    const CFreal UgUg_n = uxg_n*uxg_n + uyg_n*uyg_n;
    const CFreal Hg_i = (gamma/(gamma - 1)*pg_i + 0.5*rhog_i*UgUg_i)/rhog_i;
    const CFreal Hg_n = (gamma/(gamma - 1)*pg_n + 0.5*rhog_n*UgUg_n)/rhog_n;
    const CFreal Ag_i = sqrt(gamma*pg_i/rhog_i);
    const CFreal Ag_n = sqrt(gamma*pg_n/rhog_n);

    (m_ghostSolPhysData)[endEM]            = rhog_total;
    (m_ghostSolPhysData)[firstSpecies]     = rhog_i/rhog_total;
    (m_ghostSolPhysData)[firstSpecies + 1] = rhog_n/rhog_total;

  CFLog(VERBOSE, "\n\n\n Here 1 \n");

    (m_ghostSolPhysData)[firstVelocity]     = uxg_i;
    (m_ghostSolPhysData)[firstVelocity + 1] = uyg_i;
    (m_ghostSolPhysData)[firstVelocity + 2] = uxg_n;
    (m_ghostSolPhysData)[firstVelocity + 3] = uyg_n;
    (m_ghostSolPhysData)[firstTemperature]     = Tg_i;
    (m_ghostSolPhysData)[firstTemperature + 1] = Tg_n;  //pg_i;
    (m_ghostSolPhysData)[firstTemperature + 2] = Ag_i;
    (m_ghostSolPhysData)[firstTemperature + 3] = Ag_n;  //Hg_i;
    (m_ghostSolPhysData)[firstTemperature + 4] = pg_i;
    (m_ghostSolPhysData)[firstTemperature + 5] = pg_n;
    (m_ghostSolPhysData)[firstTemperature + 6] = Hg_i;
    (m_ghostSolPhysData)[firstTemperature + 7] = Hg_n;

  CFLog(VERBOSE, "\n\n\n Here 2 \n");

    m_varSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
 }
  
}
//////////////////////////////////////////////////////////////////////////////

void BCSubInletTtPtAlphaEIWRhoiViTi2D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                            std::vector< std::vector< RealVector* > >& ghostGrads,
                                            const std::vector< RealVector >& normals,
                                            const std::vector< RealVector >& coords)
{

  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // number of gradient variables
  cf_assert(nbrStateGrads > 0);

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    cf_assert(intGrads[iState].size() == 4);

    // normal
    const RealVector& normal = normals[iState];

    // pressure
    RealVector& presGradI = *intGrads  [iState][0];
    RealVector& presGradG = *ghostGrads[iState][0];
    const CFreal nPresGrad = (presGradI[XX]*normal[XX] + presGradI[YY]*normal[YY]);
    presGradG = presGradI - 2.0*nPresGrad*normal;

    // temperature
    RealVector& tempGradI = *intGrads  [iState][4];
    RealVector& tempGradG = *ghostGrads[iState][4];
    const CFreal nTempGrad = (tempGradI[XX]*normal[XX] + tempGradI[YY]*normal[YY]);
    tempGradG = tempGradI - 2.0*nTempGrad*normal;

    // velocity
    RealVector& uGradI = *intGrads  [iState][1];
    RealVector& uGradG = *ghostGrads[iState][1];
    RealVector& vGradI = *intGrads  [iState][2];
    RealVector& vGradG = *ghostGrads[iState][2];
    //RealVector& wGradI = *intGrads  [iState][3];
    //RealVector& wGradG = *ghostGrads[iState][3];

    // internal normal and tangential component
    const RealVector velocNGradI = uGradI*normal[XX] + vGradI*normal[YY];
    const RealVector uTGradI = uGradI - velocNGradI*normal[XX];
    const RealVector vTGradI = vGradI - velocNGradI*normal[YY];
    //const RealVector wTGradI = wGradI - velocNGradI*normal[ZZ];

    // ghost normal and tangential component
    const RealVector velocNGradG = velocNGradI;
    RealVector velocTGradNI(3);
    velocTGradNI[XX] = uTGradI[XX]*normal[XX] + uTGradI[YY]*normal[YY]; // + uTGradI[ZZ]*normal[ZZ];
    velocTGradNI[YY] = vTGradI[XX]*normal[XX] + vTGradI[YY]*normal[YY]; // + vTGradI[ZZ]*normal[ZZ];
    //velocTGradNI[ZZ] = wTGradI[XX]*normal[XX] + wTGradI[YY]*normal[YY] + wTGradI[ZZ]*normal[ZZ];
    const RealVector uTGradG = uTGradI - 2.0*velocTGradNI[XX]*normal;
    const RealVector vTGradG = vTGradI - 2.0*velocTGradNI[YY]*normal;
    //const RealVector wTGradG = wTGradI - 2.0*velocTGradNI[ZZ]*normal;

    // compute ghost velocity gradients
    uGradG = uTGradG + velocNGradG*normal[XX];
    vGradG = vTGradG + velocNGradG*normal[YY];
    //wGradG = wTGradG + velocNGradG*normal[ZZ];
  }

}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletTtPtAlphaEIWRhoiViTi2D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  // get varset
  m_varSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
  if (m_varSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> in BCSubInletEulerTtPtAlpha2D!");
  }

  // resize the physical data for internal and ghost solution points
  m_varSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_varSet->getModel()->resizePhysicalData(m_intSolPhysData);

  // non-dimensionalize pressure and temperature
  //m_tTotal /= m_eulerVarSet->getModel()->getTempRef ();
  //m_pTotal /= m_eulerVarSet->getModel()->getPressRef();
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

