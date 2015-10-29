#include "Environment/ObjectProvider.hh"
#include "Common/StringOps.hh"
#include "MultiFluidMHD/MultiFluidMHD.hh"
#include "MultiFluidMHD/EulerMFMHD2DCons.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerMFMHD2DCons, ConvectiveVarSet, MultiFluidMHDModule, 1>
mf2DConsProvider("EulerMFMHD2DCons");

//////////////////////////////////////////////////////////////////////////////

EulerMFMHD2DCons::EulerMFMHD2DCons(Common::SafePtr<BaseTerm> term) :
  MultiFluidMHDVarSet<Maxwell2DProjectionVarSet>(term),
  _rightEv(),
  _leftEv(),
  _m_i()
{    
  const CFuint endEM = 8;
  const CFuint nbSpecies    = getModel()->getNbScalarVars(0);
  const CFuint nbMomentum   = getModel()->getNbScalarVars(1);
  const CFuint nbEnergyEqs  = getModel()->getNbScalarVars(2);
  const CFuint totalNbEqs = nbSpecies + nbMomentum + nbEnergyEqs + 8;   
  
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
    names[endEM + nbSpecies + 2*ie]     = "rhoU" + StringOps::to_str(ie);
    names[endEM + nbSpecies + 2*ie + 1] = "rhoV" + StringOps::to_str(ie);
  } 
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[endEM + nbSpecies + nbMomentum + ie] = "rhoE" + StringOps::to_str(ie);
  } 
  
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

EulerMFMHD2DCons::~EulerMFMHD2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD2DCons::setup()
{
  MultiFluidMHDVarSet<Maxwell2DProjectionVarSet>::setup();
  
  setConstJacob();
  _m_i.resize(3);
}
      
//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD2DCons::setConstJacob()
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD2DCons::computeProjectedJacobian(const RealVector& normal,
						   RealMatrix& jacob)
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD2DCons::computeJacobians()
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD2DCons::computeEigenValuesVectors(RealMatrix& rightEv,
					    RealMatrix& leftEv,
					    RealVector& eValues,
					    const RealVector& normal)
{
}


//////////////////////////////////////////////////////////////////////////////

CFuint EulerMFMHD2DCons::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD2DCons::splitJacobian(RealMatrix& jacobPlus,
					RealMatrix& jacobMin,
					RealVector& eValues,
					const RealVector& normal)
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD2DCons::computePhysicalData(const State& state, RealVector& data)
{
  
  const CFuint nbSpecies    = getModel()->getNbScalarVars(0);
  const CFuint nbMomentum   = getModel()->getNbScalarVars(1);
  const CFuint nbEnergyEqs  = getModel()->getNbScalarVars(2);
  const CFuint endEM = 8;
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);  
  const CFuint firstVelocity = getModel()->getFirstScalarVar(1);   
  const CFuint firstTemperature = getModel()->getFirstScalarVar(2); 
  
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
    //set the species velocities in 2D
    const CFreal ui = state[endEM + nbSpecies + 2*ie]/rhoi;
    const CFreal vi = state[endEM + nbSpecies + 2*ie + 1]/rhoi;
    data[firstVelocity + 2*ie] = ui;
    data[firstVelocity + 2*ie + 1] = vi;
 
    //set the energy physical data
    const CFuint firstTemperature = getModel()->getFirstScalarVar(2);
    
    const CFreal mi = _m_i[ie];    
    const CFreal V2 = ui*ui + vi*vi;
    const CFreal c_p = (gamma/(gamma-1))*(K_gas/mi);
    const CFreal R_gas = K_gas/mi;
    const CFreal c_v = c_p - R_gas;
    const CFreal Ti = (state[endEM + nbSpecies + 2*nbSpecies + ie] - rhoi*V2)/(rhoi*c_v);
    
    data[firstTemperature + 4*ie] = Ti;//Temperature
    data[firstTemperature + 4*ie + 1] = Ti*R_gas*rhoi;//pressure
    data[firstTemperature + 4*ie + 2] = sqrt(gamma*R_gas*Ti);//sound speed
    data[firstTemperature + 4*ie + 3] = 0.5*V2 + c_p*Ti;//total enthaply of species i   
  }    
  //cout << "EulerMFMHD2DCons::computePhysicalData" << endl;
  //for (CFuint ie = 0; ie < firstTemperature + 4*nbSpecies; ++ie) {
    //cout << "data["<< ie <<"] = "<< data[ie] << endl;
  //}
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD2DCons::computeStateFromPhysicalData(const RealVector& data,
					       State& state)
{
   
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
  
  //set the species velocities in 2D 
  const CFuint firstVelocity = getModel()->getFirstScalarVar(1);  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    state[endEM + nbSpecies + 2*ie] = data[firstVelocity + 2*ie]*data[firstSpecies + ie]*rho;
    state[endEM + nbSpecies + 2*ie + 1] = data[firstVelocity + 2*ie + 1]*data[firstSpecies + ie]*rho;    
    
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
    const CFreal V2 = data[firstVelocity + 2*ie]*data[firstVelocity + 2*ie] + data[firstVelocity + 2*ie + 1]*data[firstVelocity + 2*ie + 1];
    const CFreal c_p = (gamma/(gamma-1))*(K_gas/_m_i[ie]);
    const CFreal R_gas = K_gas/_m_i[ie];
    const CFreal c_v = c_p - R_gas;
    
    state[endEM + nbSpecies + 2*nbSpecies + ie] = data[firstSpecies + ie]*rho*(c_v*data[firstTemperature + 4*ie] + 0.5*V2);
  } 
  //cout << "EulerMFMHD2DCons::computeStateFromPhysicalData" << endl;
  //for (CFuint ie = 0; ie < endEM + 3*nbSpecies; ++ie) {
    //cout << "state["<< ie <<"] = "<< state[ie] << endl;
  //}
}

//////////////////////////////////////////////////////////////////////////////

CFreal EulerMFMHD2DCons::getSpeed(const State& state) const
{
 
  return 0.;
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD2DCons::setDimensionalValues(const State& state,
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

void EulerMFMHD2DCons::setAdimensionalValues(const State& state,
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

void EulerMFMHD2DCons::computePerturbedPhysicalData(const Framework::State& state,
					       const RealVector& pdataBkp,
					       RealVector& pdata,
					       CFuint iVar)
{
  throw Common::NotImplementedException
      (FromHere(), "EulerMFMHD2DCons::computePerturbedPhysicalData() not implemented");
}

//////////////////////////////////////////////////////////////////////////////

bool EulerMFMHD2DCons::isValid(const RealVector& data)
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
