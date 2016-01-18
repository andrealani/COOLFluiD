#include "MultiFluidMHD/MultiFluidMHD.hh"
#include "EulerMFMHD2DHalfConsToRhoiViTi.hh"
#include "Framework/PhysicalModel.hh"
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

Environment::ObjectProvider<EulerMFMHD2DHalfConsToRhoiViTi, VarSetTransformer, MultiFluidMHDModule, 1> eulerMFMHD2DHalfConsToRhoiViTiProvider("EulerMFMHD2DHalfConsToRhoiViTi");

//////////////////////////////////////////////////////////////////////////////

EulerMFMHD2DHalfConsToRhoiViTi::EulerMFMHD2DHalfConsToRhoiViTi(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerMFMHDTerm>())
{  
  _m_i.resize(3);
}

//////////////////////////////////////////////////////////////////////////////

EulerMFMHD2DHalfConsToRhoiViTi::~EulerMFMHD2DHalfConsToRhoiViTi()
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD2DHalfConsToRhoiViTi::transform(const State& state, State& result)
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
  
  
  //The densities continue wthout transformation
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    result[endEM + ie] = state[endEM + ie];
  }
  
  //Set the velocities Ui = RhoiUi/Rhoi, Vi = RhoiVi/Rhoi  and Wi = RhoiWi/Rhoi
  const CFuint dim = 3;
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) { 
    result[endEM + nbSpecies + dim*ie]     = state[endEM + nbSpecies + dim*ie]/state[endEM + ie];
    result[endEM + nbSpecies + dim*ie + 1] = state[endEM + nbSpecies + dim*ie + 1]/state[endEM + ie];
    result[endEM + nbSpecies + dim*ie + 2] = state[endEM + nbSpecies + dim*ie + 2]/state[endEM + ie];
  }
  
  const bool isLeake = _model->isLeake();

  // plasma + neutrals model
  if(isLeake){

    //ions
    //set the energy parameters
    const CFreal gamma = _model->getGamma();	// gamma = 5/3
    const CFreal K_B = _model->getK();		// Boltzmann constant
    
    const CFreal m_p = _model->getMolecularMass1();
    
    const CFreal R_i = K_B/m_p;				// ions gas constant
    const CFreal R_p = 2*R_i;				// Plasma gas constant (ions + electrons)
    const CFreal Cp_p = (gamma/(gamma-1))*R_p;	
    const CFreal Cv_p = (1/(gamma-1))*R_p;
    const CFreal u_i = result[endEM + nbSpecies];
    const CFreal v_i = result[endEM + nbSpecies + 1];
    const CFreal w_i = result[endEM + nbSpecies + 2];
    const CFreal V2_i = u_i*u_i + v_i*v_i + w_i*w_i;
    
    result[endEM + nbSpecies + dim*nbSpecies] = (state[endEM + (1 + dim)*nbSpecies]/state[endEM] - 0.5*V2_i)/Cv_p; 
    
    //neutrals
    //set the energy parameters
    const CFreal m_n = _model->getMolecularMass2();
    
    const CFreal R_n = K_B/m_n;				// neutrals gas constant
    const CFreal Cp_n = (gamma/(gamma-1))*R_n;	
    const CFreal Cv_n = (1/(gamma-1))*R_n;
    const CFreal u_n = result[endEM + nbSpecies + 3];
    const CFreal v_n = result[endEM + nbSpecies + 4];
    const CFreal w_n = result[endEM + nbSpecies + 5];
    const CFreal V2_n = u_n*u_n + v_n*v_n + w_n*w_n;
    
    result[endEM +  nbSpecies + dim*nbSpecies + 1] = (state[endEM + (1 + dim)*nbSpecies + 1]/state[endEM + 1] - 0.5*V2_n)/Cv_n;    
    
  }
  
  else{   
    //set the energy parameters
    const CFreal gamma = _model->getGamma();
    const CFreal K_gas = _model->getK();
    const CFreal m_e = _model->getMolecularMass3();
    const CFreal m_n = _model->getMolecularMass2();
    const CFreal m_p = _model->getMolecularMass1(); 
    
    //set the molar masses of the species (should be changed in the future)

    _m_i[0] = m_e;
    _m_i[1] = m_n;
    _m_i[2] = m_p;  
    
    //Set Ti = (RhoiEi/Rhoi - 0.5*V²)/Cv
    
    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    
      const CFreal u = state[endEM + nbSpecies + dim*ie]/state[endEM + ie];
      const CFreal v = state[endEM + nbSpecies + dim*ie + 1]/state[endEM + ie];
      const CFreal w = state[endEM + nbSpecies + dim*ie + 2]/state[endEM + ie];
      const CFreal V2 = u*u + v*v + w*w;
      const CFreal c_p = (gamma/(gamma-1))*(K_gas/_m_i[ie]);
      const CFreal R_gas = K_gas/_m_i[ie];
      const CFreal c_v = c_p - R_gas;  
      
      result[endEM + (1 + dim)*nbSpecies + ie] = (state[endEM + (1 + dim)*nbSpecies + ie]/state[endEM + ie] - 0.5*V2)/c_v; 
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void EulerMFMHD2DHalfConsToRhoiViTi::transformFromRef(const RealVector& data, State& result)
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
  
  //Set the momentum Ui, Vi and Wi
  const CFuint dim = 3;  
  const CFuint firstVelocity = _model->getFirstScalarVar(1);  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) { 
    result[endEM + nbSpecies + dim*ie]     = data[firstVelocity + dim*ie];
    result[endEM + nbSpecies + dim*ie + 1] = data[firstVelocity + dim*ie + 1];
    result[endEM + nbSpecies + dim*ie + 2] = data[firstVelocity + dim*ie + 2];
  }
    
  
  //Set Ti 
  const CFuint firstTemperature = _model->getFirstScalarVar(2);  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
  
    result[endEM + (1 + dim)*nbSpecies + ie] = data[firstTemperature + 4*ie]; 
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
