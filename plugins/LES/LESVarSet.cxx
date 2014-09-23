#include "LES/LESVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

/// global variable containing a pointer to the instance of the class above
/// @note this pointer will be set in the setup function of the instance
/// @see LES2DVarSet and LES3DVarSet
/// @warning there should be only one instance of the LESVarSet class
/// This is because there is only one pointer to this instance.
COOLFluiD::Common::SafePtr< COOLFluiD::LES::LESVarSet > lesVarSetInstance = CFNULL;

//////////////////////////////////////////////////////////////////////////////

extern "C" {

#include "les_interface.h"

/// definition of some functions declared in les_interface.h
/// provided by the solver

//////////////////////////////////////////////////////////////////////////////

void les_set_stencil_size ( unsigned int stencil_size )
{
}

//////////////////////////////////////////////////////////////////////////////

void les_put_volume(double* volume)
{
  lesVarSetInstance->putVolume(volume);
}

//////////////////////////////////////////////////////////////////////////////

void les_put_velocity(double* vel)
{
  lesVarSetInstance->putVelocity(vel);
}

//////////////////////////////////////////////////////////////////////////////

void les_put_grad_velocity(double* grad_vel)
{
  lesVarSetInstance->putVelGrads(grad_vel);
}

//////////////////////////////////////////////////////////////////////////////

void les_put_temperature(double* temperature)
{
  lesVarSetInstance->putTemperature(temperature);
}

//////////////////////////////////////////////////////////////////////////////

void les_put_grad_temperature(double* grad_t)
{
  lesVarSetInstance->putTempGrad(grad_t);
}

//////////////////////////////////////////////////////////////////////////////

void les_put_density(double* density)
{
  lesVarSetInstance->putDensity(density);
}

//////////////////////////////////////////////////////////////////////////////

void les_put_cp(double* cp)
{
  lesVarSetInstance->putCp(cp);
}

//////////////////////////////////////////////////////////////////////////////

void les_put_wall_distance(double* distance)
{
  *distance = 0.;
/// @todo implement this one
}

} // extern "C"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace LES {

//////////////////////////////////////////////////////////////////////////////

LESVarSet::LESVarSet(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  Physics::NavierStokes::NavierStokesVarSet(name, model),
  m_dimensionality(),
  m_nbrEqs(),
  m_nbrEqsMin1(),
  m_eulerModel(model->getConvectiveTerm().d_castTo<Physics::NavierStokes::EulerTerm>()),
  m_currState(),
  m_currGrad(),
  m_stateVolume()
{
}

//////////////////////////////////////////////////////////////////////////////

LESVarSet::~LESVarSet() {}

//////////////////////////////////////////////////////////////////////////////

void LESVarSet::setup()
{
  Physics::NavierStokes::NavierStokesVarSet::setup();

  // get the dimensionality and the number of equations
  m_dimensionality = PhysicalModelStack::getActive()->getDim ();
  m_nbrEqs         = PhysicalModelStack::getActive()->getNbEq();
  m_nbrEqsMin1 = m_nbrEqs - 1;

  // set pointer to current object global variable
  lesVarSetInstance = this;

  // initialize the les library
  les_initialize (m_dimensionality);
}

//////////////////////////////////////////////////////////////////////////////

void LESVarSet::putVelocity(double* velocity)
{
  for (CFuint iDir = 0; iDir < m_dimensionality; ++iDir)
  {
    velocity[iDir] = _gradState[iDir+1];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LESVarSet::putVelGrads(double* velGrads)
{
  CFuint idx = 0;
  for (CFuint iVel = 0; iVel < m_dimensionality; ++iVel)
  {
    for (CFuint iDir = 0; iDir < m_dimensionality; ++iDir, ++idx)
    {
      velGrads[idx] = (*(*m_currGrad)[iVel+1])[iDir];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void LESVarSet::putTemperature(double* temperature)
{
  *temperature = _gradState[m_nbrEqsMin1];
}

//////////////////////////////////////////////////////////////////////////////

void LESVarSet::putTempGrad(double* tempGrad)
{
  for (CFuint iDir = 0; iDir < m_dimensionality; ++iDir)
  {
    tempGrad[iDir] = (*(*m_currGrad)[m_nbrEqsMin1])[iDir];
  }
}

//////////////////////////////////////////////////////////////////////////////

void LESVarSet::putDensity(double* density)
{
  *density = getDensity(*m_currState);
}

//////////////////////////////////////////////////////////////////////////////

void LESVarSet::putCp(double* cp)
{
  const CFreal gamma = m_eulerModel->getGamma();
  *cp = gamma*m_eulerModel->getR()/(gamma - 1.);
}

//////////////////////////////////////////////////////////////////////////////

void LESVarSet::putVolume(double* vol)
{
  *vol = getVolume();
}

//////////////////////////////////////////////////////////////////////////////

CFreal LESVarSet::getDynSGSViscosity(const RealVector& state,
                                     const vector<RealVector*>& gradients,
                                     const CFreal& volume)
{
  // set the gradient variables in the gradient state (not the gradients themselves!)
  // in other words, transforms state to [p u v w T] and stores in _gradState
  // this is not necessary since the eddy_viscosity model only needs gradients
  // setGradientState(state);
  
  // If state, gradients, and volume are dimensional,
  // it will return dimensional dynamic SGS viscosity.
  // Otherwise the adimensional dynamic SGS viscosity.
  
  // set pointer to current state
  m_currState = &state;

  // set pointer to current gradients (to enable passing the data to the library
  m_currGrad = &gradients;
  
  // set the volume of the cell
  setVolume(volume);
    
  // get eddy viscosity and thermal conductivity
  CFreal eddyDynVisc   = 0.0;
  les_compute_eddy_dynvisc(&eddyDynVisc  );
  
  return eddyDynVisc;
}

//////////////////////////////////////////////////////////////////////////////

CFreal LESVarSet::getSGSKineticEnergy(const RealVector& state,
                                      const vector<RealVector*>& gradients,
                                      const CFreal& volume)
{
  // set the gradient variables in the gradient state (not the gradients themselves!)
  // in other words, transforms state to [p u v w T] and stores in _gradState
  // this is not necessary since the eddy_viscosity model only needs gradients
  // setGradientState(state);
  
  // If state, gradients, and volume are dimensional,
  // it will return dimensional SGS kinetic energy
  // Otherwise the adimensional SGS kinetic energy
  CFreal Cv = 0.094;
  CFreal rho = state[0]/(m_eulerModel->getRdim()*state[3]);
  CFreal nuT = getDynSGSViscosity(state,gradients,volume)/rho;
                                       
  // set pointer to current state
  m_currState = &state;
  CFreal dim = static_cast<CFreal>(PhysicalModelStack::getActive()->getDim());
  return std::pow(nuT/(Cv*std::pow(volume,1./dim)),2);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LES

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
