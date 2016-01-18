#include "MultiFluidMHD/MultiFluidMHD.hh"
#include "EulerMFMHD3DRhoiViTiToCons.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/PhysicalConsts.hh"
#include "Environment/ObjectProvider.hh"
#include "EulerMFMHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerMFMHD3DRhoiViTiToCons, VarSetTransformer, MultiFluidMHDModule, 1> eulerMFMHD3DRhoiViTiToConsProvider("EulerMFMHD3DRhoiViTiToCons");

//////////////////////////////////////////////////////////////////////////////

EulerMFMHD3DRhoiViTiToCons::EulerMFMHD3DRhoiViTiToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerMFMHDTerm>()),
  _m_i()
{  
  _m_i.resize(3);
}

//////////////////////////////////////////////////////////////////////////////

EulerMFMHD3DRhoiViTiToCons::~EulerMFMHD3DRhoiViTiToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD3DRhoiViTiToCons::transform(const State& state, State& result)
{
  const CFuint nbSpecies    = _model->getNbScalarVars(0);
  const CFuint nbMomentum   = _model->getNbScalarVars(1);
  const CFuint nbEnergyEqs  = _model->getNbScalarVars(2);
  const CFuint endEM = 8;
  
  //Electro Magnetic Field Needs no tranformation   
  result[0] = state[0];
  result[1] = state[1];  
  result[2] = state[2];
  result[3] = state[3];
  result[4] = state[4];
  result[5] = state[5]; 
  result[6] = state[6];
  result[7] = state[7];  
  
  //The densities continue without transformation
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    result[endEM + ie] = state[endEM + ie];
  }
  
  //Set the momentum RhoiUi, RhoiVi and RhoiWi
  const CFuint dim = 3;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) { 
    const CFreal rho_i = state[endEM + ie];
    const CFreal uID   = endEM + nbSpecies + dim*ie;
    const CFreal vID   = uID + 1; 
    const CFreal wID   = uID + 2; 
    result[uID] = rho_i*state[uID];
    result[vID] = rho_i*state[vID];
    result[wID] = rho_i*state[wID];
  }

  const bool isLeake = _model->isLeake();
  
  if (isLeake){
    //ions
    const CFreal m_p = _model->getMolecularMass1();    
    const CFreal gamma = _model->getGamma();
    const CFreal K_B = PhysicalConsts::Boltzmann(); 
    
    const CFreal R_i = K_B/m_p;				// ions gas constant
    const CFreal R_p = 2.*R_i;				// Plasma gas constant (ions + electrons)
    const CFreal Cv_p = R_p/(gamma-1.);
    const CFreal Cp_p = gamma*Cv_p;
    
    const CFreal u_i = state[endEM + nbSpecies];
    const CFreal v_i = state[endEM + nbSpecies + 1];
    const CFreal w_i = state[endEM + nbSpecies + 2];
    const CFreal V2_i = u_i*u_i + v_i*v_i + w_i*w_i;
    const CFuint TiID = endEM + (1 + dim)*nbSpecies;
    
    result[TiID] = state[endEM]*(Cv_p*state[TiID] + 0.5*V2_i);    
    
    //neutrals
    const CFreal m_n = _model->getMolecularMass2();
    
    const CFreal R_n = K_B/m_n;				// neutrals gas constant
    const CFreal Cv_n = R_n/(gamma-1.);      // Cv for neutrals 
    const CFreal Cp_n = gamma*Cv_n;          // Cp for neutrals
    
    const CFreal u_n = state[endEM + nbSpecies + 3];
    const CFreal v_n = state[endEM + nbSpecies + 4];
    const CFreal w_n = state[endEM + nbSpecies + 5];
    const CFreal V2_n = u_n*u_n + v_n*v_n + w_n*w_n;
    const CFuint TnID = TiID + 1;
    
    result[TnID] = state[endEM + 1]*(Cv_n*state[TnID] + 0.5*V2_n);
//     std::cout << "EulerMFMHD3DRhoiViTiToCons::transform() => result = "<< result <<"\n";
    
  }
  else{
    //set the energy parameters
    const CFreal gamma = _model->getGamma();
    const CFreal K_gas = _model->getK();
    const CFreal m_e = _model->getMolecularMass3();
    const CFreal m_n = _model->getMolecularMass2();
    const CFreal m_p = _model->getMolecularMass1(); 
    const CFuint dim = 3;
    
  //   set the molar masses of the species (should be changed in the future)

    _m_i[0] = m_e;
    _m_i[1] = m_n;
    _m_i[2] = m_p;  
    
  //   Set RhoiEi = Rhoi*(Cv*Ti + 0.5*V²)
    
    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    
      const CFreal V2 = state[nbSpecies + dim*ie]*state[nbSpecies + dim*ie] + state[nbSpecies + dim*ie + 1]*state[nbSpecies + dim*ie + 1] + 
		        state[nbSpecies + dim*ie + 2]*state[nbSpecies + dim*ie + 2];
      const CFreal c_p = (gamma/(gamma-1))*(K_gas/_m_i[ie]);
      const CFreal R_gas = K_gas/_m_i[ie];
      const CFreal c_v = c_p - R_gas;  

      result[endEM + nbSpecies + nbMomentum + ie] = state[endEM + ie]*(c_v*state[endEM + (1 + dim)*nbSpecies + ie] + 0.5*V2); 
      
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD3DRhoiViTiToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFuint nbSpecies    = _model->getNbScalarVars(0);
  const CFuint nbMomentum   = _model->getNbScalarVars(1);
  const CFuint nbEnergyEqs  = _model->getNbScalarVars(2);
  const CFuint endEM = 8;

  //Electro Magnetic Field Needs no tranformation   
  result[0] = data[0];
  result[1] = data[1];  
  result[2] = data[2];
  result[3] = data[3];
  result[4] = data[4];
  result[5] = data[5]; 
  result[6] = data[6];
  result[7] = data[7];  

  //Density of the mixture
  const CFreal rho = data[PTERM::RHO];
  //The densities Rhoi = Rho*yi

  const CFuint firstSpecies = _model->getFirstScalarVar(0);

  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    result[endEM + ie] = data[firstSpecies + ie]*rho;
  }

  //Set the momentum RhoiUi and RhoiVi  
  const CFuint firstVelocity = _model->getFirstScalarVar(1);
  const CFuint dim = 3;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) { 
    result[endEM + nbSpecies + dim*ie]     = data[firstVelocity + dim*ie]*rho*data[firstSpecies + ie];
    result[endEM + nbSpecies + dim*ie + 1] = data[firstVelocity + dim*ie + 1]*rho*data[firstSpecies + ie];
    result[endEM + nbSpecies + dim*ie + 2] = data[firstVelocity + dim*ie + 2]*rho*data[firstSpecies + ie];
  }

  const bool isLeake = _model->isLeake();

  if (isLeake){
    const CFuint firstTemperature = _model->getFirstScalarVar(2);
    //ions
    const CFreal m_p = _model->getMolecularMass1();    
    const CFreal gamma = _model->getGamma();
    const CFreal K_B = _model->getK(); 
    
    const CFreal R_i = K_B/m_p;				// ions gas constant
    const CFreal R_p = 2*R_i;				// Plasma gas constant (ions + electrons)
    const CFreal Cp_p = (gamma/(gamma-1))*R_p;	
    const CFreal Cv_p = (1/(gamma-1))*R_p;  
    
    const CFreal u_i = data[firstVelocity];
    const CFreal v_i = data[firstVelocity + 1];
    const CFreal w_i = data[firstVelocity + 2];
    const CFreal V2_i = u_i*u_i + v_i*v_i + w_i*w_i;
    
    result[endEM + (1 + dim)*nbSpecies] = data[firstSpecies]*rho*(Cv_p*data[firstTemperature] + 0.5*V2_i);    
    //neutrals
    const CFreal m_n = _model->getMolecularMass2();
    
    const CFreal R_n = K_B/m_n;				// neutrals gas constant
    const CFreal Cp_n = (gamma/(gamma-1))*R_n;	
    const CFreal Cv_n = (1/(gamma-1))*R_n;
    const CFreal u_n = data[firstVelocity + 3];
    const CFreal v_n = data[firstVelocity + 4];
    const CFreal w_n = data[firstVelocity + 5]; 
    const CFreal V2_n = u_n*u_n + v_n*v_n + w_n*w_n;
    
    result[endEM + (1 + dim)*nbSpecies + 1] = data[firstSpecies + 1]*rho*(Cv_n*data[firstTemperature + 4] + 0.5*V2_n);
//     std::cout << "EulerMFMHD3DRhoiViTiToCons::transformFromRef() => result = "<< result <<"\n";
  }
 
  else {
    //set the energy parameters
    const CFreal gamma = _model->getGamma();
    const CFreal K_gas = _model->getK();
    const CFreal m_e = _model->getMolecularMass3();
    const CFreal m_n = _model->getMolecularMass2();
    const CFreal m_p = _model->getMolecularMass1(); 
    const CFuint firstTemperature = _model->getFirstScalarVar(2);
    
    //set the molar masses of the species (should be changed in the future)

    _m_i[0] = m_e;
    _m_i[1] = m_n;
    _m_i[2] = m_p;  
    
    //Set RhoiEi = Rhoi*(Cv*Ti + 0.5*V²)
    
    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    
      const CFreal V2 = data[firstVelocity + dim*ie]*data[firstVelocity + dim*ie] + data[firstVelocity + dim*ie + 1]*data[firstVelocity + dim*ie + 1] +
			data[firstVelocity + dim*ie + 2]*data[firstVelocity + dim*ie + 2];
      const CFreal c_p = (gamma/(gamma-1))*(K_gas/_m_i[ie]);
      const CFreal R_gas = K_gas/_m_i[ie];
      const CFreal c_v = c_p - R_gas;  
      
      result[endEM + nbSpecies + nbMomentum + ie] = rho*data[firstSpecies + ie]*(c_v* data[firstTemperature + 4*ie] + 0.5*V2); 
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
