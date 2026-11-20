#include <algorithm>

#include "Common/CFLog.hh"
#include "MHD/MHD.hh"
#include "MHD/MHDProjectionDiffPrimRes.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MHDProjectionDiffPrimRes, DiffusiveVarSet, MHDModule, 2>
mhd2DProjectionDiffPrimResProvider("MHD2DProjectionDiffPrimRes");

Environment::ObjectProvider<MHDProjectionDiffPrimRes, DiffusiveVarSet, MHDModule, 2>
mhd3DProjectionDiffPrimResProvider("MHD3DProjectionDiffPrimRes");

//////////////////////////////////////////////////////////////////////////////

void MHDProjectionDiffPrimRes::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal, Config::DynamicOption<> >
    ("eta","Magnetic resistivity. Default 0.");
  options.addConfigOption< CFreal, Config::DynamicOption<> >
    ("nu","Kinematic viscosity. Default 0.");
}

//////////////////////////////////////////////////////////////////////////////

MHDProjectionDiffPrimRes::MHDProjectionDiffPrimRes(const std::string& name,
                                                   Common::SafePtr<Framework::PhysicalModelImpl> model) :
  MHDProjectionDiffVarSet(name, model),
  _dim(0)
{
  addConfigOptionsTo(this);

  _eta = 0.0;
  setParameter("eta", &_eta);

  _nu = 0.0;
  setParameter("nu", &_nu);

  vector<std::string> names(9);
  names[0] = "rho";
  names[1] = "u";
  names[2] = "v";
  names[3] = "w";
  names[4] = "Bx";
  names[5] = "By";
  names[6] = "Bz";
  names[7] = "p";
  names[8] = "phi";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

MHDProjectionDiffPrimRes::~MHDProjectionDiffPrimRes()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHDProjectionDiffPrimRes::setup()
{
  MHDProjectionDiffVarSet::setup();

  _dim = PhysicalModelStack::getActive()->getDim();

  CFLog(INFO, "MHDProjectionDiffPrimRes: " << _dim << "D, eta = " << _eta << ", nu = " << _nu << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void MHDProjectionDiffPrimRes::setGradientVars(const vector<RealVector*>& states,
                                                RealMatrix& values,
                                                const CFuint stateSize)
{
  // Copy all state variables; for index 7 convert p -> T_adim = p/rho
  for (CFuint i = 0; i < 9; ++i) {
    for (CFuint j = 0; j < stateSize; ++j) {
      if (i == 7) {
        values(i,j) = (*states[j])[i] / (*states[j])[0];
      } else {
        values(i,j) = (*states[j])[i];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void MHDProjectionDiffPrimRes::setGradientVarGradients(const vector<RealVector*>& states,
                                                        const vector< vector<RealVector*> >& stateGradients,
                                                        vector< vector<RealVector*> >& gradVarGradients,
                                                        const CFuint stateSize)
{
  throw Common::NotImplementedException
    (FromHere(), "MHDProjectionDiffPrimRes::setGradientVarGradients()");
}

//////////////////////////////////////////////////////////////////////////////

void MHDProjectionDiffPrimRes::setStateGradients(const vector<RealVector*>& states,
                                                  const vector< vector<RealVector*> >& gradVarGradients,
                                                  vector< vector<RealVector*> >& stateGradients,
                                                  const CFuint stateSize)
{
  throw Common::NotImplementedException
    (FromHere(), "MHDProjectionDiffPrimRes::setStateGradients()");
}

//////////////////////////////////////////////////////////////////////////////

CFreal MHDProjectionDiffPrimRes::getDynViscosity(const RealVector& state,
                                                  const vector<RealVector*>& gradients)
{
  // effective dynamic diffusivity, used only for the diffusive time-step estimate
  return state[0]*std::max(_eta, _nu);
}

//////////////////////////////////////////////////////////////////////////////

CFreal MHDProjectionDiffPrimRes::getDensity(const RealVector& state)
{
  return state[0];
}

//////////////////////////////////////////////////////////////////////////////

void MHDProjectionDiffPrimRes::setGradientState(const RealVector& state)
{
  cf_assert(_gradState.size() == state.size());
  _gradState = state;
  // Convert p -> T_adim = p/rho (consistent with setGradientVars)
  _gradState[7] = state[7] / state[0];
}

//////////////////////////////////////////////////////////////////////////////

void MHDProjectionDiffPrimRes::computeFluxJacobian(const RealVector& state,
                                                    const RealVector& gradientJacob,
                                                    const RealVector& normal,
                                                    const CFreal& radius,
                                                    RealMatrix& fluxJacob)
{
  throw Common::NotImplementedException
    (FromHere(), "MHDProjectionDiffPrimRes::computeFluxJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void MHDProjectionDiffPrimRes::setComposition(const RealVector& state,
                                               const bool isPerturb,
                                               const CFuint iVar)
{
}

//////////////////////////////////////////////////////////////////////////////

RealVector& MHDProjectionDiffPrimRes::getFlux(const RealVector& state,
                                               const vector<RealVector*>& gradients,
                                               const RealVector& normal,
                                               const CFreal& radius)
{
  setGradientState(state);

  const CFreal rho = state[0];

  // dynamic diffusivity used for the diffusive time-step estimate
  (getModel().getPhysicalData())[MHDProjectionDiffTerm::MU] = rho*std::max(_eta, _nu);

  _flux = 0.;

  // gradients[1..3] = grad(u, v, w), gradients[4..6] = grad(Bx, By, Bz)
  // viscosity: F(rho v_i) = nu rho grad(v_i) . n
  CFreal viscFlux[3] = {0., 0., 0.};
  if (_nu > 0.) {
    for (CFuint i = 0; i < 3; ++i) {
      const RealVector& gradV = *gradients[1+i];
      CFreal gradVn = 0.;
      for (CFuint j = 0; j < _dim; ++j) {
        gradVn += gradV[j]*normal[j];
      }
      viscFlux[i] = _nu*rho*gradVn;
      _flux[1+i] = viscFlux[i];
    }
  }

  // resistivity: F(B_i) = eta (grad(B_i) . n - sum_j d_i B_j n_j), with d_z = 0 in 2D
  CFreal resFlux[3] = {0., 0., 0.};
  if (_eta > 0.) {
    for (CFuint i = 0; i < 3; ++i) {
      const RealVector& gradB = *gradients[4+i];
      CFreal f = 0.;
      for (CFuint j = 0; j < _dim; ++j) {
        f += gradB[j]*normal[j];
      }
      if (i < _dim) {
        CFreal dBn = 0.;
        for (CFuint j = 0; j < _dim; ++j) {
          dBn += (*gradients[4+j])[i]*normal[j];
        }
        f -= dBn;
      }
      resFlux[i] = _eta*f;
      _flux[4+i] = resFlux[i];
    }
  }

  // energy: F(E) = B . F(B) + v . F(rho v)
  _flux[7] = state[4]*resFlux[0] + state[5]*resFlux[1] + state[6]*resFlux[2] +
             state[1]*viscFlux[0] + state[2]*viscFlux[1] + state[3]*viscFlux[2];

  return _flux;
}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& MHDProjectionDiffPrimRes::getFlux(const RealVector& state,
                                               const vector<RealVector*>& gradients,
                                               const CFreal& radius)
{
  throw Common::NotImplementedException
    (FromHere(), "MHDProjectionDiffPrimRes::getFlux() (flux vector)");
  return _fluxVec;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
