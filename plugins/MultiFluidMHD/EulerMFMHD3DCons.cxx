#include "Environment/ObjectProvider.hh"
#include "Common/StringOps.hh"
#include "MultiFluidMHD/MultiFluidMHD.hh"
#include "MultiFluidMHD/EulerMFMHD3DCons.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerMFMHD3DCons, ConvectiveVarSet, MultiFluidMHDModule, 1>
mf3DConsProvider("EulerMFMHD3DCons");

//////////////////////////////////////////////////////////////////////////////

EulerMFMHD3DCons::EulerMFMHD3DCons(Common::SafePtr<BaseTerm> term) :
  MultiFluidMHDVarSet<Maxwell3DProjectionVarSet>(term),
  _rightEv(),
  _leftEv(),
  _m_i()
{    
  const CFuint endEM = 8;
  const CFuint nbSpecies    = getModel()->getNbScalarVars(0);
  const CFuint nbMomentum   = getModel()->getNbScalarVars(1);
  const CFuint nbEnergyEqs  = getModel()->getNbScalarVars(2);
  const CFuint totalNbEqs = nbSpecies + nbMomentum + nbEnergyEqs + 8;  
  const CFuint dim = 3;
  
  vector<std::string> names(totalNbEqs);
  
  names[0] = "Bx";
  names[1] = "By";
  names[2] = "Bz";
  names[3] = "Ex";
  names[4] = "Ey";
  names[5] = "Ez";
  names[6] = "Psi";
  names[7] = "Phi";
  
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[endEM + ie] = "rho" + StringOps::to_str(ie);
  }

  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[endEM + nbSpecies + dim*ie]     = "rhoU" + StringOps::to_str(ie);
    names[endEM + nbSpecies + dim*ie + 1] = "rhoV" + StringOps::to_str(ie);
    names[endEM + nbSpecies + dim*ie + 2] = "rhoW" + StringOps::to_str(ie);
  } 
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[endEM + nbSpecies + nbMomentum + ie] = "rhoE" + StringOps::to_str(ie);
  } 
  
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

EulerMFMHD3DCons::~EulerMFMHD3DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD3DCons::setup()
{
  MultiFluidMHDVarSet<Maxwell3DProjectionVarSet>::setup();
  
  setConstJacob();
  _m_i.resize(3);
}
      
//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD3DCons::setConstJacob()
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD3DCons::computeProjectedJacobian(const RealVector& normal,
						   RealMatrix& jacob)
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD3DCons::computeJacobians()
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD3DCons::computeEigenValuesVectors(RealMatrix& rightEv,
					    RealMatrix& leftEv,
					    RealVector& eValues,
					    const RealVector& normal)
{
}


//////////////////////////////////////////////////////////////////////////////

CFuint EulerMFMHD3DCons::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD3DCons::splitJacobian(RealMatrix& jacobPlus,
					RealMatrix& jacobMin,
					RealVector& eValues,
					const RealVector& normal)
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD3DCons::computePhysicalData(const State& state, RealVector& data)
{
  
  const CFuint nbSpecies    = getModel()->getNbScalarVars(0);
  const CFuint nbMomentum   = getModel()->getNbScalarVars(1);
  const CFuint nbEnergyEqs  = getModel()->getNbScalarVars(2);
  const CFuint endEM = 8;
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);  
  const CFuint firstVelocity = getModel()->getFirstScalarVar(1);   
  const CFuint firstTemperature = getModel()->getFirstScalarVar(2); 

  const bool isLeake = getModel()->isLeake();

  // plasma + neutrals model
  if(isLeake){ std::cout<<"EulerMFMHD3DCons::computePhysicalData NOT IMPLEMENTED \n";}
 
///Used for debugging     
  
  data[PTERM::BX] = state[0];
  data[PTERM::BY] = state[1]; 
  data[PTERM::BZ] = state[2];  
  data[PTERM::EX] = state[3];  
  data[PTERM::EY] = state[4];  
  data[PTERM::EZ] = state[5];  
  data[PTERM::PSI] = state[6];  
  data[PTERM::PHI] = state[7];  
  
  
  //set the total density
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    rho += state[endEM + ie];
  }

  data[PTERM::RHO] = rho;
  const CFreal ovRho = 1./rho;
  
  //set the energy parameters
  const CFreal gamma = getModel()->getGamma();
  const CFreal K_gas = getModel()->getK();
  const CFreal m_e = getModel()->getMolecularMass1();
  const CFreal m_n = getModel()->getMolecularMass2();
  const CFreal m_p = getModel()->getMolecularMass3(); 
  
  //set the molar masses of the species (should be changed in the future)
  _m_i[0] = m_e;
  _m_i[1] = m_n;
  _m_i[2] = m_p;
  
  
  //set the species mass fraction
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    const CFreal rhoi = state[endEM + ie];
    data[firstSpecies + ie] = rhoi*ovRho;
  //set the species velocities in 3D 
//   const CFuint firstVelocity = getModel()->getFirstScalarVar(1); 
    const CFuint dim = 3;
    const CFreal ui = state[endEM + nbSpecies + dim*ie]/rhoi;
    const CFreal vi = state[endEM + nbSpecies + dim*ie + 1]/rhoi;
    const CFreal wi = state[endEM + nbSpecies + dim*ie + 2]/rhoi;
    data[firstVelocity + dim*ie] = ui;
    data[firstVelocity + dim*ie + 1] = vi;
    data[firstVelocity + dim*ie + 2] = wi;
 
  //set the energy physical data  
  const CFuint firstTemperature = getModel()->getFirstScalarVar(2);
    
    
    const CFreal mi = _m_i[ie];    
    const CFreal V2 = ui*ui + vi*vi + wi*wi;
    const CFreal c_p = (gamma/(gamma-1))*(K_gas/mi);
    const CFreal R_gas = K_gas/mi;
    const CFreal c_v = c_p - R_gas;
    const CFreal Ti = (state[endEM + nbSpecies + dim*nbSpecies + ie] - rhoi*V2)/(rhoi*c_v);
    
    data[firstTemperature + 4*ie] = Ti;//Temperature
    data[firstTemperature + 4*ie + 1] = Ti*R_gas*rhoi;//pressure
    data[firstTemperature + 4*ie + 2] = sqrt(gamma*R_gas*Ti);//sound speed
    data[firstTemperature + 4*ie + 3] = 0.5*V2 + c_p*Ti;//total enthaply of species i   
  }    
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD3DCons::computeStateFromPhysicalData(const RealVector& data,
					       State& state)
{

  const bool isLeake = getModel()->isLeake();

  // plasma + neutrals model
  if(isLeake){ std::cout<<"EulerMFMHD3DCons::computePhysicalData NOT IMPLEMENTED \n";}

  const CFuint nbSpecies    = getModel()->getNbScalarVars(0);
  const CFuint endEM = 8;  
  
  state[0] = data[PTERM::BX];
  state[1] = data[PTERM::BY]; 
  state[2] = data[PTERM::BZ];  
  state[3] = data[PTERM::EX];  
  state[4] = data[PTERM::EY];  
  state[5] = data[PTERM::EZ];  
  state[6] = data[PTERM::PSI];  
  state[7] = data[PTERM::PHI];  
  
  const CFreal rho = data[PTERM::RHO];
  
  //set the species densities
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    state[endEM + ie] = data[firstSpecies + ie]*rho; //rhoi = Rho*yi
  }  
  
  //set the species velocities in 3D 
  const CFuint firstVelocity = getModel()->getFirstScalarVar(1);
  const CFuint dim = 3;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    state[endEM + nbSpecies + dim*ie]     = data[firstVelocity + dim*ie]*data[firstSpecies + ie]*rho;
    state[endEM + nbSpecies + dim*ie + 1] = data[firstVelocity + dim*ie + 1]*data[firstSpecies + ie]*rho;    
    state[endEM + nbSpecies + dim*ie + 2] = data[firstVelocity + dim*ie + 2]*data[firstSpecies + ie]*rho;    
  }
  //set the energy parameters
  const CFreal gamma = getModel()->getGamma();
  const CFreal K_gas = getModel()->getK();
  const CFreal m_e = getModel()->getMolecularMass1();
  const CFreal m_n = getModel()->getMolecularMass2();
  const CFreal m_p = getModel()->getMolecularMass3(); 
  
  //set the molar masses of the species (should be changed in the future)
  _m_i[0] = m_e;
  _m_i[1] = m_n;
  _m_i[2] = m_p;
 
 //set the Energies
  const CFuint firstTemperature = getModel()->getFirstScalarVar(2);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    const CFreal V2 = data[firstVelocity + dim*ie]*data[firstVelocity + dim*ie] + data[firstVelocity + dim*ie + 1]*data[firstVelocity + dim*ie + 1] + 
		      data[firstVelocity + dim*ie + 2]*data[firstVelocity + dim*ie + 2];
    const CFreal c_p = (gamma/(gamma-1))*(K_gas/_m_i[ie]);
    const CFreal R_gas = K_gas/_m_i[ie];
    const CFreal c_v = c_p - R_gas;
    
    state[endEM + nbSpecies + dim*nbSpecies + ie] = data[firstSpecies + ie]*rho*(c_v*data[firstTemperature + 4*ie] + 0.5*V2);
  }       
}

//////////////////////////////////////////////////////////////////////////////

CFreal EulerMFMHD3DCons::getSpeed(const State& state) const
{
 
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD3DCons::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies    = getModel()->getNbScalarVars(0);
  const CFuint nbMomentum   = getModel()->getNbScalarVars(1);
  const CFuint nbEnergyEqs  = getModel()->getNbScalarVars(2);
  const CFuint totalNbEquations = nbSpecies + nbMomentum + nbEnergyEqs + 8; 
  
  for (CFuint i = 0; i < totalNbEquations; ++i) {
    result[i] = state[i]*refData[i];   
  }
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD3DCons::setAdimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFuint nbSpecies    = getModel()->getNbScalarVars(0);
  const CFuint nbMomentum   = getModel()->getNbScalarVars(1);
  const CFuint nbEnergyEqs  = getModel()->getNbScalarVars(2);
  const CFuint totalNbEquations = nbSpecies + nbMomentum + nbEnergyEqs + 8; 
  
  for (CFuint i = 0; i < totalNbEquations; ++i) {
    result[i] = state[i]/refData[i];   
  }  

}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD3DCons::computePerturbedPhysicalData(const Framework::State& state,
					       const RealVector& pdataBkp,
					       RealVector& pdata,
					       CFuint iVar)
{
  throw Common::NotImplementedException
      (FromHere(), "EulerMFMHD3DCons::computePerturbedPhysicalData() not implemented");
}

//////////////////////////////////////////////////////////////////////////////

bool EulerMFMHD3DCons::isValid(const RealVector& data)
{
 //  bool correct = true;
//   enum index {RHO, RHOU,RHOV,RHOE};

//   const CFreal rho = data[RHO];
//   const CFreal ovRho = 1./rho;
//   const CFreal u = data[RHOU]*ovRho;
//   const CFreal v = data[RHOV]*ovRho;
//   const CFreal V2 = u*u + v*v;

//   const CFreal gamma = getModel()->getGamma();
//   const CFreal gammaMinus1 = gamma - 1.;
//   /// Next call returns R/M, in dimensional units.
//   const CFreal R = getModel()->getR();

//   const CFreal p = gammaMinus1*(data[RHOE] - 0.5*rho*V2);

//   const CFreal T = p*ovRho/R;
  
//   const CFreal a2 = gamma*p*ovRho;

//   if( ( p < 0.) || ( T < 0.) || ( a2 < 0.) ){
//   return correct = false;
//   }

//   return correct;
  return true;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
