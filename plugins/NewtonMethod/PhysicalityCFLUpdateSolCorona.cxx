// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonMethodMHD.hh"
#include "NewtonMethod/PhysicalityCFLUpdateSolCorona.hh"

#include "Common/BadValueException.hh"
#include "Common/CFLog.hh"
#include "Common/PE.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "MathTools/MathConsts.hh"

#include <algorithm>
#include <cmath>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PhysicalityCFLUpdateSolCorona, NewtonIteratorData, NewtonMethodMHDModule>
physicalityCFLUpdateSolCoronaProvider("PhysicalityCFLUpdateSolCorona");

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSolCorona::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("DensityBC","Reference density of the corona inner boundary, sets the Alfven speed cap (densityBoundaryValue of StdUpdateSolPP2).");
  options.addConfigOption< CFreal >("AlfvenFactor","Multiplier of the reference Alfven speed in the density cap.");
  options.addConfigOption< CFreal >("MinTemperature","Temperature floor in K.");
  options.addConfigOption< CFreal >("DampingRadius","Radius below which the velocity damping applies.");
  options.addConfigOption< CFreal >("MeanMolecularWeight","Mean molecular weight used in the temperature.");
  options.addConfigOption< std::vector<CFreal> >("RefValues","Corona normalization [Bref, rhoref, pref, vref] in SI units.");
}

//////////////////////////////////////////////////////////////////////////////

PhysicalityCFLUpdateSolCorona::PhysicalityCFLUpdateSolCorona(const std::string& name) :
  PhysicalityCFLUpdateSol(name),
  m_refValues(),
  m_bMagMax(1.),
  m_layoutChecked(false)
{
  addConfigOptionsTo(this);

  m_rhoBC = 0.;
  setParameter("DensityBC",&m_rhoBC);

  m_alfvenFactor = 2.;
  setParameter("AlfvenFactor",&m_alfvenFactor);

  m_minT = 1.0e4;
  setParameter("MinTemperature",&m_minT);

  m_dampingRadius = 1.102;
  setParameter("DampingRadius",&m_dampingRadius);

  m_mu = 1.27;
  setParameter("MeanMolecularWeight",&m_mu);

  setParameter("RefValues",&m_refValues);
}

//////////////////////////////////////////////////////////////////////////////

PhysicalityCFLUpdateSolCorona::~PhysicalityCFLUpdateSolCorona()
{
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSolCorona::setup()
{
  PhysicalityCFLUpdateSol::setup();

  if (m_refValues.size() == 0) {
    // COCONUT corona normalization: Bref [T], rhoref [kg/m3], pref [Pa], vref [m/s]
    m_refValues.resize(4);
    m_refValues[0] = 2.2e-4;
    m_refValues[1] = 1.67e-13;
    m_refValues[2] = 0.03851;
    m_refValues[3] = 4.80e5;
  }
  if (m_refValues.size() != 4) {
    throw BadValueException(FromHere(),
      "PhysicalityCFLUpdateSolCorona: RefValues needs 4 entries [Bref, rhoref, pref, vref]");
  }
  for (CFuint i = 0; i < m_refValues.size(); ++i) {
    if (!(m_refValues[i] > 0.)) {
      throw BadValueException(FromHere(), "PhysicalityCFLUpdateSolCorona: RefValues must be > 0");
    }
  }
  if (!(m_rhoBC > 0.)) {
    throw BadValueException(FromHere(),
      "PhysicalityCFLUpdateSolCorona: DensityBC must be > 0 (the densityBoundaryValue of the StdUpdateSolPP2 case)");
  }
  if (!(m_alfvenFactor > 0.)) {
    throw BadValueException(FromHere(), "PhysicalityCFLUpdateSolCorona: AlfvenFactor must be > 0");
  }
  if (!(m_minT >= 0.)) {
    throw BadValueException(FromHere(), "PhysicalityCFLUpdateSolCorona: MinTemperature must be >= 0");
  }
  if (!(m_mu > 0.)) {
    throw BadValueException(FromHere(), "PhysicalityCFLUpdateSolCorona: MeanMolecularWeight must be > 0");
  }
  // the limiters read [rho, u, v, w, Bx, By, Bz, p, ...] out of a 3D state
  if (PhysicalModelStack::getActive()->getDim() != DIM_3D ||
      PhysicalModelStack::getActive()->getNbEq() < 8) {
    throw BadValueException(FromHere(),
      "PhysicalityCFLUpdateSolCorona: a 3D model with at least 8 equations is needed");
  }
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSolCorona::checkStateLayout()
{
  const std::vector<std::string>& names = m_varSet->getVarNames();
  if (names.size() < 8 || names[0] != "rho" || names[7] != "p") {
    std::string found;
    for (CFuint i = 0; i < names.size(); ++i) { found += names[i] + " "; }
    throw BadValueException(FromHere(),
      "PhysicalityCFLUpdateSolCorona: the update variables have to be ordered as "
      "[rho, u, v, w, Bx, By, Bz, p, ...] (MHD3DProjectionPrim or PrimE), found [" +
      found + "]");
  }
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSolCorona::beforeUpdate()
{
  if (!m_layoutChecked) {
    checkStateLayout();
    m_layoutChecked = true;
  }

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  const CFuint nbStates = states.size();

  // floor of 1 in normalized field units, as in StdUpdateSolPP2
  m_bMagMax = 1.;
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    const State& state = *states[iState];
    if (state.isParUpdatable()) {
      const CFreal bMag = std::sqrt(state[4]*state[4] + state[5]*state[5] +
                                    state[6]*state[6]);
      m_bMagMax = std::max(m_bMagMax, bMag);
    }
  }

#ifdef CF_HAVE_MPI
  if (PE::GetPE().IsParallel()) {
    // StdUpdateSolPP2 keeps this maximum rank local, which makes its density
    // cap depend on the partitioning
    const std::string nsp = getMethodData().getNamespace();
    const CFreal bMagMaxLocal = m_bMagMax;
    MPI_Allreduce(&bMagMaxLocal, &m_bMagMax, 1, MPI_DOUBLE, MPI_MAX,
                  PE::GetPE().GetCommunicator(nsp));
  }
#endif
}

//////////////////////////////////////////////////////////////////////////////

void PhysicalityCFLUpdateSolCorona::afterUpdate()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  SafePtr<FilterState> filterState = getMethodData().getFilterState();

  const CFreal bRef   = m_refValues[0];
  const CFreal rhoRef = m_refValues[1];
  const CFreal pRef   = m_refValues[2];
  const CFreal vRef   = m_refValues[3];
  const CFreal mu0 = 1.2566e-6;  // vacuum permeability [H/m]
  const CFreal kB  = 1.38e-23;   // Boltzmann constant [J/K]
  const CFreal mH  = 1.67e-27;   // hydrogen mass [kg]
  // widths of the two tanh ramps, in m/s and in K: each cap is approached
  // smoothly over its width, so a state sitting just under one is left alone
  const CFreal dVA = 2.0e3;
  const CFreal dT  = 2.0e1;
  // damping ramp of StdUpdateSolPP2: a floor of 0.3 at the inner boundary,
  // reaching 1 about 0.1 Rs above it
  const CFreal rampFloor = 0.3;
  const CFreal rampRate  = 8.68;
  // the damping sound speed uses this gamma in StdUpdateSolPP2, whatever
  // gamma the physical model itself runs with
  const CFreal gammaDamping = 1.667;
  const CFreal pi = MathConsts::CFrealPi();

  const CFreal maxVA = m_alfvenFactor*m_bMagMax*bRef/std::sqrt(m_rhoBC*rhoRef*mu0);

  const CFuint nbStates = states.size();
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    State& state = *states[iState];
    if (!state.isParUpdatable()) { continue; }

    // 1. Alfven speed cap: raise the density wherever the local Alfven speed
    //    runs past maxVA. The 1e-14 keeps vA finite in a field free state.
    const CFreal bMag = std::sqrt(state[4]*state[4] + state[5]*state[5] +
                                  state[6]*state[6]) + 1.e-14;
    const CFreal rhoDim = state[0]*rhoRef;
    const CFreal vA = bMag*bRef/std::sqrt(rhoDim*mu0);
    CFreal s = 0.5 + 0.5*std::tanh(pi*(vA - maxVA)/dVA);
    state[0] = (s*(bMag*bRef)*(bMag*bRef)/(maxVA*maxVA)/mu0 + rhoDim*(1. - s))/rhoRef;

    // 2. temperature floor at fixed density. n kB = 2 rho kB/(mu mH) counts
    //    one electron per ion, as the COCONUT temperature does.
    const CFreal nkB = 2.*state[0]*rhoRef*kB/(m_mu*mH);
    const CFreal pDim = state[7]*pRef;
    const CFreal T = pDim/nkB;
    s = 0.5 + 0.5*std::tanh(pi*(m_minT - T)/dT);
    state[7] = (s*m_minT*nkB + pDim*(1. - s))/pRef;

    // 3. velocity damping below the damping radius: nothing faster than the
    //    local sound speed at the inner boundary, released with the height
    const Node& coord = state.getCoordinates();
    const CFreal r = std::sqrt(coord[XX]*coord[XX] + coord[YY]*coord[YY] +
                               coord[ZZ]*coord[ZZ]);
    if (r < m_dampingRadius) {
      const CFreal vDim = std::sqrt(state[1]*state[1] + state[2]*state[2] +
                                    state[3]*state[3])*vRef;
      if (vDim > 0.) {
        const CFreal cSound = std::sqrt(gammaDamping*state[7]*pRef/(state[0]*rhoRef));
        const CFreal ramp = std::min(1., rampFloor + std::tanh((r - 1.)*rampRate));
        const CFreal damping = std::min(1., (cSound/vDim)*ramp);
        state[1] *= damping;
        state[2] *= damping;
        state[3] *= damping;
      }
    }

    // the base class has already filtered the updated state, so the clipping
    // above has to pass the same filter
    filterState->filter(state);
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD
