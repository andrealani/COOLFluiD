#include "Framework/MethodStrategyProvider.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/DiffMFMHD3DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Maxwell/Maxwell3DVarSet.hh"
#include "Framework/PhysicalConsts.hh"
#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "BCSubInletTtPtAlphaEIWRhoiViTi3D.hh"
#include "Common/NotImplementedException.hh"

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

Framework::MethodStrategyProvider<BCSubInletTtPtAlphaEIWRhoiViTi3D,FluxReconstructionSolverData,BCStateComputer,FluxReconstructionMultiFluidMHDModule>                                                                                           
  BCSubInletTtPtAlphaEIWRhoiViTi3DProvider("SubInletTtPtAlphaEIWRhoiViTi3D");

//////////////////////////////////////////////////////////////////////////////

void BCSubInletTtPtAlphaEIWRhoiViTi3D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Ttot","total temperature");
  options.addConfigOption< CFreal >("Ptot","total pressure");
  options.addConfigOption< CFreal >("alphaXY","alphaXY");
  options.addConfigOption< CFreal >("alphaXZ","alphaXZ");
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletTtPtAlphaEIWRhoiViTi3D::BCSubInletTtPtAlphaEIWRhoiViTi3D(const std::string& name) :
  BCStateComputer(name),
  m_varSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData(),
  m_tTotal(),
  m_pTotal(),
  m_alphaXY(),
  m_alphaXZ()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  m_tTotal = std::vector<CFreal>();;
   setParameter("Ttot",&m_tTotal);

  m_pTotal = std::vector<CFreal>();;
   setParameter("Ptot",&m_pTotal);

  m_alphaXY = std::vector<CFreal>();;
   setParameter("alphaXY",&m_alphaXY);

  m_alphaXZ = std::vector<CFreal>();;
   setParameter("alphaXZ",&m_alphaXZ);
}

//////////////////////////////////////////////////////////////////////////////

BCSubInletTtPtAlphaEIWRhoiViTi3D::~BCSubInletTtPtAlphaEIWRhoiViTi3D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletTtPtAlphaEIWRhoiViTi3D::computeGhostStates(const vector< State* >& intStates,
                                                    vector< State* >& ghostStates,
                                                    const std::vector< RealVector >& normals,
                                                    const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

  // get some data from the physical model
  const CFreal gamma = m_varSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma -1.0;
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);
  //const CFuint idGassConst = m_varSet->getModel()->getR();

  //const CFreal tgAlphaXY = tan(m_alphaXY);
  //const CFreal tgAlphaXZ = tan(m_alphaXZ);

  ///MultiFluidMHD Subsonic Outlet Condition imposing pressure in 2D
  // here a fix is needed in order to have always ghostT > 0
  // if ghostT < 0  then the inner value is set
  const CFuint endEM = 8;
  const CFuint firstSpecies = m_varSet->getModel()->getFirstScalarVar(0);
  const CFuint firstVelocity = m_varSet->getModel()->getFirstScalarVar(1);
  const CFuint firstTemperature = m_varSet->getModel()->getFirstScalarVar(2);


    //cf_assert(m_tTotal.size() == nbSpecies);
    //cf_assert(m_pTotal.size() == nbSpecies);

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
      // dereference states
      State& intState   = (*intStates[iState]);
      State& ghostState = (*ghostStates[iState]);

      // set the physical data starting from the inner state
      m_varSet->computePhysicalData(intState,m_intSolPhysData);


    (m_ghostSolPhysData)[0] = (m_intSolPhysData)[0] /*+ 2*bn*nx*/;	//Bx
    (m_ghostSolPhysData)[1] = (m_intSolPhysData)[1] /*+ 2*bn*ny*/;	//By
    (m_ghostSolPhysData)[2] = (m_intSolPhysData)[2] /*+ 2*bn*nz*/;	//Bz
    (m_ghostSolPhysData)[3] = (m_intSolPhysData)[3] /*- 2*en*nx*/;	//Ex
    (m_ghostSolPhysData)[4] = (m_intSolPhysData)[4] /*- 2*en*ny*/;	//Ey
    (m_ghostSolPhysData)[5] = (m_intSolPhysData)[5] /*- 2*en*nz*/;	//Ez
    (m_ghostSolPhysData)[6] = (m_intSolPhysData)[6];			//Psi
    (m_ghostSolPhysData)[7] = (m_intSolPhysData)[7];			//Phi

    const bool isLeake = true;

    cf_assert(intState.size() == 5);
    cf_assert(ghostState.size() == 5);

  if (isLeake) {

    const CFreal m_electron = m_varSet->getModel()->getMolecularMass1();
    const CFreal m_neutral = m_varSet->getModel()->getMolecularMass2();
    const CFreal m_ion = m_varSet->getModel()->getMolecularMass3();
    const CFreal K_B = PhysicalConsts::Boltzmann();
    const CFreal Rplasma = 2*K_B/m_ion;				// ions gas constant
    const CFreal Rneutral = K_B/m_neutral;

    const CFreal uxi = (m_intSolPhysData)[firstVelocity];
    const CFreal uyi = (m_intSolPhysData)[firstVelocity + 1];
    const CFreal uzi = (m_intSolPhysData)[firstVelocity + 2];
    const CFreal uxn = (m_intSolPhysData)[firstVelocity + 3];
    const CFreal uyn = (m_intSolPhysData)[firstVelocity + 4];
    const CFreal uzn = (m_intSolPhysData)[firstVelocity + 5];
    const CFreal Ui = std::sqrt(uxi*uxi + uyi*uyi + uzi*uzi);
    const CFreal Un = std::sqrt(uxn*uxn + uyn*uyn + uzn*uzn);
    const CFreal rhoi = (m_intSolPhysData)[endEM]*(m_intSolPhysData)[firstSpecies];
    const CFreal rhon = (m_intSolPhysData)[endEM]*(m_intSolPhysData)[firstSpecies + 2];
    const CFreal Ti = (m_intSolPhysData)[firstTemperature];
    const CFreal Tn = (m_intSolPhysData)[firstTemperature + 4];
    const CFreal pi = rhoi*Rplasma*Ti;
    const CFreal pn = rhon*Rneutral*Tn;
    const CFreal ai = (m_intSolPhysData)[firstTemperature + 2];
    const CFreal an = (m_intSolPhysData)[firstTemperature + 6];
    const CFreal Mi = Ui/ai;
    const CFreal Mn = Un/an;
    const CFreal coeffi = 1 + 0.5*(gamma - 1)*Mi*Mi;
    const CFreal coeffn = 1 + 0.5*(gamma - 1)*Mn*Mn;
    const CFreal coeffPow_i = pow(coeffi, gamma/(gamma - 1));
    const CFreal coeffPow_n = pow(coeffn, gamma/(gamma - 1));
    const CFreal Ttoti = Ti*coeffi;
    const CFreal Ttotn = Tn*coeffn;
    //const CFreal coeff2 = (1 + gamma*Mi*Mi);
    const CFreal Ptoti = pi*coeffPow_i;
    const CFreal Ptotn = pn*coeffPow_n;
    CFreal tgAlphaXYi = 0.;
    CFreal tgAlphaXYn = 0.;
    CFreal tgAlphaXZi = 0.;
    CFreal tgAlphaXZn = 0.;

    if(uxi != 0.){tgAlphaXYi = uyi/uxi;}
    if(uxn != 0.){tgAlphaXYn = uyn/uxn;}

    const CFreal Ttotg_i = 2*m_tTotal[0] - Ttoti;
    const CFreal Ptotg_i = 2*m_pTotal[0] - Ptoti;

    const CFreal Ttotg_n = 2*m_tTotal[1] - Ttotn;
    const CFreal Ptotg_n = 2*m_pTotal[1] - Ptotn;

    const CFreal tgAlphaXYg_i = 2*tan(m_alphaXY[0]) - tgAlphaXYi;
    const CFreal tgAlphaXYg_n = 2*tan(m_alphaXY[1]) - tgAlphaXYn;
    const CFreal tgAlphaXZg_i = 2*tan(m_alphaXZ[0]) - tgAlphaXZi;
    const CFreal tgAlphaXZg_n = 2*tan(m_alphaXZ[1]) - tgAlphaXZn;

    const CFreal Mg_i = Mi;
    const CFreal Mg_n = Mn;
    const CFreal Tg_i = Ttotg_i/coeffi;
    const CFreal pg_i = Ptotg_i/coeffPow_i;
    const CFreal Tg_n = Ttotg_n/coeffn;
    const CFreal pg_n = Ptotg_n/coeffPow_n;
    const CFreal rhog_i = pg_i/(Rplasma*Tg_i);
    const CFreal rhog_n = pg_n/(Rneutral*Tg_n);
    const CFreal rhog_total = rhog_i + rhog_n;

    const CFreal uxg_i = Mg_i*std::sqrt(gamma*Rplasma*Tg_i/(1 + tgAlphaXYg_i*tgAlphaXYg_i));
    const CFreal uyg_i = tgAlphaXYg_i*uxg_i;
    const CFreal uzg_i = tgAlphaXZg_i*uxg_i;
    const CFreal uxg_n = Mg_n*std::sqrt(gamma*Rneutral*Tg_n/(1 + tgAlphaXYg_n*tgAlphaXYg_n));
    const CFreal uyg_n = tgAlphaXYg_n*uxg_n;
    const CFreal uzg_n = tgAlphaXZg_n*uxg_n;

    const CFreal UgUg_i = uxg_i*uxg_i + uyg_i*uyg_i + uzg_i*uzg_i;
    const CFreal UgUg_n = uxg_n*uxg_n + uyg_n*uyg_n + uzg_n*uzg_n;
    const CFreal Hg_i = (gamma/(gamma - 1)*pg_i + 0.5*rhog_i*UgUg_i)/rhog_i;
    const CFreal Hg_n = (gamma/(gamma - 1)*pg_n + 0.5*rhog_n*UgUg_n)/rhog_n;
    const CFreal Ag_i = sqrt(gamma*pg_i/rhog_i);
    const CFreal Ag_n = sqrt(gamma*pg_n/rhog_n);

    (m_ghostSolPhysData)[endEM]            = rhog_total;
    (m_ghostSolPhysData)[firstSpecies]     = rhog_i/rhog_total;
    (m_ghostSolPhysData)[firstSpecies + 1] = rhog_n/rhog_total;

    (m_ghostSolPhysData)[firstVelocity]     = uxg_i;
    (m_ghostSolPhysData)[firstVelocity + 1] = uyg_i;
    (m_ghostSolPhysData)[firstVelocity + 2] = uzg_i;
    (m_ghostSolPhysData)[firstVelocity + 3] = uxg_n;
    (m_ghostSolPhysData)[firstVelocity + 4] = uyg_n;
    (m_ghostSolPhysData)[firstVelocity + 5] = uzg_n;
    (m_ghostSolPhysData)[firstTemperature]     = Tg_i;
    (m_ghostSolPhysData)[firstTemperature + 1] = pg_i;
    (m_ghostSolPhysData)[firstTemperature + 2] = Ag_i;
    (m_ghostSolPhysData)[firstTemperature + 3] = Hg_i;
    (m_ghostSolPhysData)[firstTemperature + 4] = Tg_n;
    (m_ghostSolPhysData)[firstTemperature + 5] = pg_n;
    (m_ghostSolPhysData)[firstTemperature + 6] = Ag_n;
    (m_ghostSolPhysData)[firstTemperature + 7] = Hg_n;

    m_varSet->computeStateFromPhysicalData(m_ghostSolPhysData, ghostState);
  }
 }
}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletTtPtAlphaEIWRhoiViTi3D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
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
    cf_assert(intGrads[iState].size() == 5);

    // normal
    const RealVector& normal = normals[iState];

    // pressure
    RealVector& presGradI = *intGrads  [iState][0];
    RealVector& presGradG = *ghostGrads[iState][0];
    const CFreal nPresGrad = (presGradI[XX]*normal[XX] + presGradI[YY]*normal[YY] + presGradI[ZZ]*normal[ZZ]);
    presGradG = presGradI - 2.0*nPresGrad*normal;

    // temperature
    RealVector& tempGradI = *intGrads  [iState][4];
    RealVector& tempGradG = *ghostGrads[iState][4];
    const CFreal nTempGrad = (tempGradI[XX]*normal[XX] + tempGradI[YY]*normal[YY] + tempGradI[ZZ]*normal[ZZ]);
    tempGradG = tempGradI - 2.0*nTempGrad*normal;

    // velocity
    RealVector& uGradI = *intGrads  [iState][1];
    RealVector& uGradG = *ghostGrads[iState][1];
    RealVector& vGradI = *intGrads  [iState][2];
    RealVector& vGradG = *ghostGrads[iState][2];
    RealVector& wGradI = *intGrads  [iState][3];
    RealVector& wGradG = *ghostGrads[iState][3];

    // internal normal and tangential component
    const RealVector velocNGradI = uGradI*normal[XX] + vGradI*normal[YY] + wGradI*normal[ZZ];
    const RealVector uTGradI = uGradI - velocNGradI*normal[XX];
    const RealVector vTGradI = vGradI - velocNGradI*normal[YY];
    const RealVector wTGradI = wGradI - velocNGradI*normal[ZZ];

    // ghost normal and tangential component
    const RealVector velocNGradG = velocNGradI;
    RealVector velocTGradNI(3);
    velocTGradNI[XX] = uTGradI[XX]*normal[XX] + uTGradI[YY]*normal[YY] + uTGradI[ZZ]*normal[ZZ];
    velocTGradNI[YY] = vTGradI[XX]*normal[XX] + vTGradI[YY]*normal[YY] + vTGradI[ZZ]*normal[ZZ];
    velocTGradNI[ZZ] = wTGradI[XX]*normal[XX] + wTGradI[YY]*normal[YY] + wTGradI[ZZ]*normal[ZZ];
    const RealVector uTGradG = uTGradI - 2.0*velocTGradNI[XX]*normal;
    const RealVector vTGradG = vTGradI - 2.0*velocTGradNI[YY]*normal;
    const RealVector wTGradG = wTGradI - 2.0*velocTGradNI[ZZ]*normal;

    // compute ghost velocity gradients
    uGradG = uTGradG + velocNGradG*normal[XX];
    vGradG = vTGradG + velocNGradG*normal[YY];
    wGradG = wTGradG + velocNGradG*normal[ZZ];
  }

}

//////////////////////////////////////////////////////////////////////////////

void BCSubInletTtPtAlphaEIWRhoiViTi3D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  // get MultiFluidMHD 3D varset
  m_varSet = getMethodData().getUpdateVar().d_castTo< MultiFluidMHDVarSet<Maxwell3DVarSet> >();
  if (m_varSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not MultiFluidMHDVarSet<Maxwell3DProjectionVarSet> in BCSubInletTtPtAlphaEIWRhoiViTi3D!");
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

