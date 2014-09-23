#include "MultiFluidMHD/MultiFluidMHD.hh"
#include "MultiFluidMHD/DiffMFMHD3DRhoiViTi.hh"
#include "Common/StringOps.hh"
#include "Environment/ObjectProvider.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MultiFluidMHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<DiffMFMHD3DRhoiViTi, DiffusiveVarSet, MultiFluidMHDModule, 2>
nsMHD3DRhoiViTiProvider("MultiFluidMHD3DRhoiViTi");

//////////////////////////////////////////////////////////////////////////////

DiffMFMHD3DRhoiViTi::DiffMFMHD3DRhoiViTi(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  DiffMFMHD3DVarSet(name, model),
  m_eulerModel(model->getConvectiveTerm().d_castTo<PTERM>()),
  m_density()
{
  const CFuint endEM = 8;
  const CFuint nbSpecies = m_eulerModel->getNbScalarVars(0);
  const CFuint nbMomentum   = m_eulerModel->getNbScalarVars(1);
  const CFuint nbEnergyEqs  = m_eulerModel->getNbScalarVars(2);
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
    names[endEM + nbSpecies + dim*ie]     = "U" + StringOps::to_str(ie);
    names[endEM + nbSpecies + dim*ie + 1] = "V" + StringOps::to_str(ie);
    names[endEM + nbSpecies + dim*ie + 2] = "W" + StringOps::to_str(ie);
  } 
  
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[endEM + nbSpecies + nbMomentum + ie] = "T" + StringOps::to_str(ie);
  } 
  
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

DiffMFMHD3DRhoiViTi::~DiffMFMHD3DRhoiViTi()
{
}

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHD3DRhoiViTi::setGradientVars(const vector<RealVector*>& states,
					  RealMatrix& values,
					  const CFuint stateSize)

{
  using namespace std;
  using namespace COOLFluiD::Framework;
  
  const CFuint nbValues = values.nbRows();
  for (CFuint i = 0; i < nbValues; ++i) {
    for (CFuint j = 0; j < stateSize; ++j) {
      values(i,j) = (*states[j])[i];
    }
  }
}

// {
//   const CFuint nbStates = stateSize;
//   const CFuint nbSpecies = m_eulerModel->getNbScalarVars(0);
//   const CFuint endEM = 8;
//   const CFreal K_gas = m_eulerModel->getK();
//   const CFreal gamma = m_eulerModel->getGamma();
//   m_m_i[0] = m_eulerModel->getMolecularMass1();
//   m_m_i[1] = m_eulerModel->getMolecularMass2();
//   m_m_i[2] = m_eulerModel->getMolecularMass3();
//   
//   cf_assert(values.nbRows() == PhysicalModelStack::getActive()->getNbEq());
//   
//   for (CFuint i = 0; i < nbStates; ++i) {
//     
//     const RealVector& state = *states[i];
//     
//     for (CFuint ie = 0; ie < nbSpecies; ++ie) {
//       
//       const CFreal mass = m_m_i[ie];
//       const CFreal R = K_gas/mass; 
//       const CFreal c_p = (gamma/(gamma-1))*R;
//       const CFreal c_v = c_p - R;
//     
//       const CFreal rhoi = state[endEM + ie];
//       const CFreal ui = state[endEM + nbSpecies + 2*ie];
//       const CFreal vi = state[endEM + nbSpecies + 2*ie + 1];
//       const CFreal V2 = ui*ui + vi*vi;
//       const CFreal Ti = state[endEM + 3*nbSpecies + ie];
//     
//        values(endEM + nbSpecies + 2*ie,i) = ui;					//U-component of ie species
//        values(endEM + nbSpecies + 2*ie + 1,i) = vi;				//V-component of ie species
//        values(endEM + nbSpecies + 2*nbSpecies + ie,i) = Ti;		//Temperature of ie species
//        values(endEM + ie,i) = R*rhoi*Ti;			//pressure = R*rhoi*T
//     }
//   }
// }

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHD3DRhoiViTi::setGradientVarGradients(const vector<RealVector*>& states,
                                                  const vector< vector<RealVector*> >& stateGradients,
                                                  vector< vector<RealVector*> >& gradVarGradients,
                                                  const CFuint stateSize)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  
  cf_assert(stateGradients.size() > 0);
  const CFuint nbValues = stateGradients[0].size();
  for (CFuint i = 0; i < nbValues; ++i)
  {
    for (CFuint j = 0; j < stateSize; ++j)
    {
      *gradVarGradients[j][i] = *stateGradients[j][i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHD3DRhoiViTi::setStateGradients(const vector<RealVector*>& states,
                                            const vector< vector<RealVector*> >& gradVarGradients,
                                            vector< vector<RealVector*> >& stateGradients,
                                            const CFuint stateSize)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  
  cf_assert(stateGradients.size() > 0);
  const CFuint nbValues = stateGradients[0].size();
  for (CFuint i = 0; i < nbValues; ++i)
    {
      for (CFuint j = 0; j < stateSize; ++j)
	{
	  *stateGradients[j][i] = *gradVarGradients[j][i];
	}
  }
}

//////////////////////////////////////////////////////////////////////////////

RealVector& DiffMFMHD3DRhoiViTi::getDynViscosityVec(const RealVector& state,
						    const vector<RealVector*>& gradients)
{
  return getModel().getDynViscosityDim();
}

//////////////////////////////////////////////////////////////////////////////

CFreal DiffMFMHD3DRhoiViTi::getDynViscosity(const RealVector& state,
                                            const vector<RealVector*>& gradients)
{
  return getModel().getDynViscosityDim().max();
}

//////////////////////////////////////////////////////////////////////////////

CFreal DiffMFMHD3DRhoiViTi::getThermConductivity(const RealVector& state,
                                      const CFreal& dynViscosity)
{
  return getModel().getThermConductivityDim().max();
}

//////////////////////////////////////////////////////////////////////////////

CFreal DiffMFMHD3DRhoiViTi::getDensity(const RealVector& state)
{
  const CFuint nbSpecies = m_eulerModel->getNbScalarVars(0);
  const CFuint endEM = 8;
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    rho += state[endEM + ie];
  }
  return rho;
}

//////////////////////////////////////////////////////////////////////////////

void DiffMFMHD3DRhoiViTi::setGradientState(const RealVector& state)
{
  cf_assert(_gradState.size() == state.size());
//   const CFuint endEM = 8;
//   const CFuint nbSpecies = m_eulerModel->getNbScalarVars(0);
//   const CFuint nbMomentum   = m_eulerModel->getNbScalarVars(1);
//   const CFuint nbEnergyEqs  = m_eulerModel->getNbScalarVars(2);
//   const CFuint totalNbEqs = nbSpecies + nbMomentum + nbEnergyEqs + 8;
  
//   const CFreal K_gas = m_eulerModel->getK();
//   const CFreal gamma = m_eulerModel->getGamma();
//   m_m_i[0] = m_eulerModel->getMolecularMass1();
//   m_m_i[1] = m_eulerModel->getMolecularMass2();
//   m_m_i[2] = m_eulerModel->getMolecularMass3();
  
//   for (CFuint i = 0; i< totalNbEqs; ++i) {
    
//     const CFreal mass = m_m_i[i];
//     const CFreal R = K_gas/mass; 
//     const CFreal c_p = (gamma/(gamma-1))*R;
//     const CFreal c_v = c_p - R;
//     
//     const CFreal rhoi = state[i];
//     const CFreal ui = state[nbSpecies + 2*i]/rhoi;
//     const CFreal vi = state[nbSpecies + 2*i + 1]/rhoi;
//     const CFreal V2 = ui*ui + vi*vi;
//     const CFreal Ti = (1/(state[i]*c_v))*(state[nbSpecies + 2*nbSpecies + i] - 0.5*state[i]*V2);
//     
//     _gradState[i] = rhoi*R*Ti;						//Pressure
//     _gradState[nbSpecies + i] = ui;					//U-component
//     _gradState[nbSpecies + 2*i] = vi;					//V-component
//     _gradState[3*nbSpecies + i] = Ti;					//Temperature
  _gradState = state;
//   }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MultiFluidMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
