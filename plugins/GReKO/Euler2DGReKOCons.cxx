#include "GReKO.hh"
#include "Euler2DGReKOCons.hh"
#include "Common/NotImplementedException.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/NSTurbTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace GReKO {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DGReKOCons, ConvectiveVarSet, GReKOModule, 1>
euler2DGReKOConsProvider("Euler2DGReKOCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DGReKOCons::Euler2DGReKOCons(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<Euler2DVarSet>(term)
{
  const CFuint nbTurbEquations = getModel()->getNbScalarVars(0);

  vector<std::string> names(4 + nbTurbEquations);
  names[0] = "rho";
  names[1] = "rhoU";
  names[2] = "rhoV";
  names[3] = "rhoE";

  names[4] = "rhoK";
  names[5] = "rhoOmega";
 
  names[6] = "rhoGa";
  names[7] = "rhoRe";

  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DGReKOCons::~Euler2DGReKOCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOCons::setup()
{
  MultiScalarVarSet<Euler2DVarSet>::setup();

  // set the equation set data for each of the equation subsets
  // first equation subset
  Euler2DVarSet::getEqSetData().resize(1);
  Euler2DVarSet::getEqSetData()[0].setup(0,0,4);

  // second equation subset
  MultiScalarVarSet<Euler2DVarSet>::getEqSetData().resize(1);
  MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[0].setup
    (1,4,getModel()->getNbScalarVars(0));

  setConstJacob();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOCons::setConstJacob()
{

///@todo fix this,
/*this is only correct for the Euler part, not the 2 turb equations



  vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  (*jacobians)[0](0,1) = 1.0;
  (*jacobians)[0](1,3) = gammaMinus1;

  (*jacobians)[1](0,2) = 1.0;
  (*jacobians)[1](2,3) = gammaMinus1;*/
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOCons::computeJacobians()
{

  ///@todo fix this,
/*this is only correct for the Euler part, not the 2 turb equations




  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal avU = linearData[EulerTerm::VX];
  const CFreal avV = linearData[EulerTerm::VY];
  const CFreal avH = linearData[EulerTerm::H];
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaMinus3 = gamma - 3.;
  const CFreal uu = avU*avU;
  const CFreal uv = avU*avV;
  const CFreal vv = avV*avV;
  const CFreal sumVel2 = 0.5*gammaMinus1*(uu + vv);

  vector<RealMatrix>* const jacobians =
    PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  (*jacobians)[0](1,0) = sumVel2 - uu;
  (*jacobians)[0](1,1) = -gammaMinus3*avU;
  (*jacobians)[0](1,2) = -gammaMinus1*avV;
  (*jacobians)[0](2,0) = -uv;
  (*jacobians)[0](2,1) = avV;
  (*jacobians)[0](2,2) = avU;
  (*jacobians)[0](3,0) = sumVel2*avU - avU*avH;
  (*jacobians)[0](3,1) = -gammaMinus1*uu + avH;
  (*jacobians)[0](3,2) = -gammaMinus1*uv;
  (*jacobians)[0](3,3) = getModel()->getGamma()*avU;

  (*jacobians)[1](1,0) = -uv;
  (*jacobians)[1](1,1) = avV;
  (*jacobians)[1](1,2) = avU;
  (*jacobians)[1](2,0) = sumVel2 - vv;
  (*jacobians)[1](2,1) = -gammaMinus1*avU;
  (*jacobians)[1](2,2) = -gammaMinus3*avV;
  (*jacobians)[1](3,0) = sumVel2*avV - avV*avH;
  (*jacobians)[1](3,1) = -gammaMinus1*uv;
  (*jacobians)[1](3,2) = -gammaMinus1*vv + avH;
  (*jacobians)[1](3,3) = getModel()->getGamma()*avV;*/
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOCons::computeEigenValuesVectors(RealMatrix& rightEv,
                                           RealMatrix& leftEv,
                                           RealVector& eValues,
                                           const RealVector& normal)
{

  //throw Common::NotImplementedException (FromHere(),"Euler2DGReKOCons::computeEigenValuesVectors()");
  ///@todo fix this, this is only correct for the Euler part, not the 2 turb equations
  const RealVector& linearData = getModel()->getPhysicalData();

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal avRho = linearData[EulerTerm::RHO];
  const CFreal avU   = linearData[EulerTerm::VX];
  const CFreal avV   = linearData[EulerTerm::VY];
  const CFreal avH   = linearData[EulerTerm::H];
  const CFreal avA   = linearData[EulerTerm::A];

  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal um = avU*nx + avV*ny;
  const CFreal ra = 0.5*avRho/avA;
  const CFreal avA2 = avA*avA;
  const CFreal coeffM2 = 0.5*gammaMinus1*(avU*avU + avV*avV)/avA2;
  const CFreal invAvRho = 1./avRho;
  const CFreal uDivA = gammaMinus1*avU/avA;
  const CFreal vDivA = gammaMinus1*avV/avA;
  const CFreal rhoA = avRho*avA;

  rightEv(0,0) = 1.;
  rightEv(0,1) = 0.;
  rightEv(0,2) = ra;
  rightEv(0,3) = ra;
  rightEv(1,0) = avU;
  rightEv(1,1) = avRho*ny;
  rightEv(1,2) = ra*(avU + avA*nx);
  rightEv(1,3) = ra*(avU - avA*nx);
  rightEv(2,0) = avV;
  rightEv(2,1) = -avRho*nx;;
  rightEv(2,2) = ra*(avV + avA*ny);
  rightEv(2,3) = ra*(avV - avA*ny);
  rightEv(3,0) = 0.5*(avU*avU +avV*avV);
  rightEv(3,1) = avRho*(avU*ny - avV*nx);
  rightEv(3,2) = ra*(avH + avA*um);
  rightEv(3,3) = ra*(avH - avA*um);

  rightEv(4,4) = 1.;
  rightEv(5,5) = 1.;

  leftEv(0,0) = 1.- coeffM2;
  leftEv(0,1) = uDivA/avA;
  leftEv(0,2) = vDivA/avA;
  leftEv(0,3) = -gammaMinus1/avA2;
  leftEv(1,0) = invAvRho*(avV*nx - avU*ny);
  leftEv(1,1) = invAvRho*ny;
  leftEv(1,2) = -invAvRho*nx;
  leftEv(1,3) = 0.0;
  leftEv(2,0) = avA*invAvRho*(coeffM2 - um/avA);
  leftEv(2,1) = invAvRho*(nx - uDivA);
  leftEv(2,2) = invAvRho*(ny - vDivA);
  leftEv(2,3) = gammaMinus1/rhoA;
  leftEv(3,0) = avA*invAvRho*(coeffM2 + um/avA);
  leftEv(3,1) = -invAvRho*(nx + uDivA);
  leftEv(3,2) = -invAvRho*(ny + vDivA);
  leftEv(3,3) = gammaMinus1/rhoA;

  leftEv(4,4) = 1.;
  leftEv(5,5) = 1.;

  eValues[0] = um;
  eValues[1] = um;
  eValues[2] = um + avA;
  eValues[3] = um - avA;
  eValues[4] = um;
  eValues[5] = um;

  //degugging infos
  CFLogDebugMax( "RightEigenvectors @Euler2DChar::computeEigenValuesVectors" << "\n"
  << rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @Euler2DChar::computeEigenValuesVectors" << "\n"
  << leftEv << "\n");
  CFLogDebugMax( "EigenValues @Euler2DChar::computeEigenValuesVectors" << "\n"
  << eValues << "\n" << "\n");

}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DGReKOCons::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOCons::splitJacobian(RealMatrix& jacobPlus,
                                RealMatrix& jacobMin,
                                RealVector& eValues,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DGReKOCons::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOCons::computePhysicalData(const State& state,
					    RealVector& data)
{
  const CFreal rho  = state[0];
  const CFreal rhoU = state[1];
  const CFreal rhoV = state[2];
  const CFreal rhoE = state[3];
  const CFreal rhoK = state[4];
  const CFreal rhoOmega = state[5];
  const CFreal rhoGa = state[6];
  const CFreal rhoRe = state[7];
  
  const CFreal overRho  = 1./rho;

  const CFreal V2 = (rhoU*rhoU + rhoV*rhoV)*overRho*overRho;
  const CFreal R = getModel()->getR();
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  data[EulerTerm::VX] = rhoU*overRho;
  data[EulerTerm::VY] = rhoV*overRho;
  data[EulerTerm::RHO] = rho;
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::P] = gammaMinus1*(rhoE - 0.5*rho*V2 - rhoK);
  data[EulerTerm::T] = data[EulerTerm::P]/(R*rho);
  data[EulerTerm::H] = (rhoE + data[EulerTerm::P])*overRho;
  data[EulerTerm::E] = rhoE/rho;

//  data[EulerTerm::A] = sqrt(gammaMinus1*(data[EulerTerm::H] - 0.5*V2) + (gammaMinus1*K);
  data[EulerTerm::A] = sqrt(gamma*data[EulerTerm::P]*overRho);
  data[EulerTerm::GAMMA] = gamma;

  const CFuint iK = getModel()->getFirstScalarVar(0);
  data[iK] = rhoK*overRho;
  data[iK+1] = rhoOmega*overRho;
  data[iK+2] = rhoGa*overRho;
  data[iK+3] = rhoRe*overRho;

}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOCons::computeStateFromPhysicalData(const RealVector& data,
						State& state)
{
  state[0] = data[EulerTerm::RHO];
  state[1] = data[EulerTerm::RHO]*data[EulerTerm::VX];
  state[2] = data[EulerTerm::RHO]*data[EulerTerm::VY];
  state[3] = data[EulerTerm::RHO]*data[EulerTerm::H] - data[EulerTerm::P];

  // Set the turbulent variables
  const CFuint firstTurbVar = getModel()->getFirstScalarVar(0);
  const CFuint nbTurbVar = getModel()->getNbScalarVars(0);

  for (CFuint ie = 0; ie < nbTurbVar; ++ie){
    state[4 + ie] = data[EulerTerm::RHO]*data[firstTurbVar + ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler2DGReKOCons::getSpeed(const State& state) const
{
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];
  return sqrt(u*u + v*v);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOCons::setDimensionalValues(const State& state,
                                          RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]*refData[EulerTerm::RHO];
  result[1] = state[1]*refData[EulerTerm::RHO]*refData[EulerTerm::V];
  result[2] = state[2]*refData[EulerTerm::RHO]*refData[EulerTerm::V];
  result[3] = state[3]*refData[EulerTerm::RHO]*refData[EulerTerm::H];

  const CFuint iK = getModel()->getFirstScalarVar(0);
  result[4] = state[4]*refData[EulerTerm::RHO]*refData[iK];
  result[5] = state[5]*refData[EulerTerm::RHO]*refData[iK+1];
  result[6] = state[6]*refData[EulerTerm::RHO]*refData[iK+2];
  result[7] = state[7]*refData[EulerTerm::RHO]*refData[iK+3];
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOCons::setAdimensionalValues(const State& state,
                                          RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  result[0] = state[0]/refData[EulerTerm::RHO];
  result[1] = state[1]/(refData[EulerTerm::RHO]*refData[EulerTerm::V]);
  result[2] = state[2]/(refData[EulerTerm::RHO]*refData[EulerTerm::V]);
  result[3] = state[3]/(refData[EulerTerm::RHO]*refData[EulerTerm::H]);

  const CFuint iK = getModel()->getFirstScalarVar(0);
  result[4] = state[4]/(refData[EulerTerm::RHO]*refData[iK]);
  result[5] = state[5]/(refData[EulerTerm::RHO]*refData[iK+1]);
  result[6] = state[6]/(refData[EulerTerm::RHO]*refData[iK+2]);
  result[7] = state[7]/(refData[EulerTerm::RHO]*refData[iK+3]);
}


//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOCons::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  const CFreal rho  = state[0];
  const CFreal rhoU = state[1];
  const CFreal rhoV = state[2];
  const CFreal rhoE = state[3];
  const CFreal rhoK = state[4];
  const CFreal rhoOmega = state[5];
  const CFreal rhoGa = state[6];
  const CFreal rhoRe = state[7];

  result[0] = rho*refData[EulerTerm::RHO];
  result[1] = rhoU*refData[EulerTerm::RHO]*refData[EulerTerm::V];
  result[2] = rhoV*refData[EulerTerm::RHO]*refData[EulerTerm::V];
  result[3] = rhoE*refData[EulerTerm::RHO]*refData[EulerTerm::H];

  const CFuint iK = getModel()->getFirstScalarVar(0);
  result[4] = rhoK*refData[EulerTerm::RHO]*refData[iK];
  result[5] = rhoOmega*refData[EulerTerm::RHO]*refData[iK+1];
  result[6] = rhoGa*refData[EulerTerm::RHO]*refData[iK+2];
  result[7] = rhoRe*refData[EulerTerm::RHO]*refData[iK+3];


  extra.resize(5);

  const CFreal V2 = (rhoU*rhoU + rhoV*rhoV)/(rho*rho);
  cout << " V2       "   << endl;
  const CFreal R = getModel()->getR();
  const CFreal gamma = getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  const CFreal p = gammaMinus1*(rhoE - 0.5*rho*V2);
  const CFreal T = p/(R*rho);

  const CFreal pdim = p * getModel()->getPressRef();
  const CFreal Tdim = T  * getModel()->getTempRef();
  const CFreal rhoDim   = result[0];
  const CFreal Kdim     = result[4];
  const CFreal Omegadim = result[5];
  const CFreal Gadim    = result[6];
  const CFreal Redim    = result[7];


   //Why the Visosities are calculated here
   //TO be Modified mut is wrong
    SafePtr<NSTurbTerm> nsTurbTerm  =
    PhysicalModelStack::getActive()->getImplementor()->
    getDiffusiveTerm().d_castTo<NSTurbTerm>();
  
  const CFreal mu_dim = nsTurbTerm->getDynViscosityDim(pdim,Tdim) * refData[NSTurbTerm::MU];
  const CFreal mut_dim = rhoDim * Kdim / Omegadim;
  
  extra[0] = mu_dim;
  extra[1] = mut_dim;
  extra[2] = mut_dim/mu_dim;

  const CFreal Tu       = 100 * (std::sqrt(2*Kdim/3))/(std::sqrt(V2));
  //const CFreal Tu        = std::max(Tu1,0.027497);
  //const CFreal Tu        = std::max(Tu1,0.00000027497);
  const CFreal overTu    = 1/Tu;
   
           if (Tu<=1.3) {
              CFreal  Rethetat = (1173.51-589.428*Tu + 0.2196*overTu*overTu);
  extra[3] = Rethetat;
                 }
          else {
                const CFreal lamco5   = Tu - 0.5658;
                const CFreal pwtu   = -0.671;
               CFreal Rethetat = 331.5*std::pow(lamco5,pwtu);
  extra[3] = Rethetat;
               }
  
  extra[4] = Tu;

}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> Euler2DGReKOCons::getExtraVarNames() const
{

  vector<std::string> names(5);
  names[0] = "mu";
  names[1] = "muT";
  names[2] = "muT/mu";
  names[3] = "Rethetat";
  names[4] = "Tu";

  return names;

}

//////////////////////////////////////////////////////////////////////////////

void Euler2DGReKOCons::computePerturbedPhysicalData
(const Framework::State& currState,
 const RealVector& bData,
 RealVector& data,
 CFuint iVar)
{
  cf_assert(iVar < PhysicalModelStack::getActive()->getNbEq());
  
  const CFuint nbEqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor().getTotalNbEqSS();
  
  if(nbEqSS == 1){
    computePhysicalData(currState, data);
  }
  
  if(nbEqSS == 2){
    if (iVar < 4) {  // u + du // v + dv // T + dT // p + dp
      computePhysicalData(currState, data);
      const CFuint iK = getModel()->getFirstScalarVar(0);
      data[iK] = bData[iK];
      data[iK+1] = bData[iK+1];
    }
  }
  else { // k + dk   omega + dOmega
    data = bData; // all the components are set == default except for K
    
    // Terms for turbulent intensity
    const CFuint iK = getModel()->getFirstScalarVar(0);
    data[iK]= currState[4]/currState[0];
    data[iK+1]= currState[5]/currState[0];
  }
  
  if(nbEqSS == 3){
    if (iVar < 4) {  // u + du // v + dv // T + dT // p + dp
      computePhysicalData(currState, data);
      const CFuint iK = getModel()->getFirstScalarVar(0);
      data[iK] = bData[iK];
      data[iK+1] = bData[iK+1];
    }
    else { // k + dk  ||  omega + dOmega
      if (iVar == 4) { // k + dk
	data = bData; // all the components are set == default except for K
	
	// Terms for turbulent intensity
	const CFuint iK = getModel()->getFirstScalarVar(0);
	data[iK]= currState[4]/currState[0];
      }
      
      if (iVar == 5) { // omega + domega
	data = bData; // all the components are set == default except for Omega
	
	// Terms for turbulent intensity
	const CFuint iK = getModel()->getFirstScalarVar(0);
	data[iK+1]= currState[5]/currState[0];
      }
    }
  }  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace GReKO

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
