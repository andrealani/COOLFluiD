// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionNavierStokes_NSBoundaryState_hh
#define COOLFluiD_FluxReconstructionNavierStokes_NSBoundaryState_hh

//////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <vector>

#include "Common/BadValueException.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/RealVector.hh"
#include "Framework/State.hh"
#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Helpers of the Navier-Stokes boundary conditions: the boundary state built
 * from prescribed primitive values, and the wall heat flux prescribed on the
 * normal temperature gradient.
 *
 * @author Rayan Dhib
 */

//////////////////////////////////////////////////////////////////////////////

/**
 * Computes a Navier-Stokes boundary state from the primitive variables
 * (p, u, v[, w], T):
 *
 *   rho = p/(R T),   H = gamma/(gamma-1) p/rho + 0.5 |u|^2,   E = H - p/rho,
 *   a = sqrt(gamma p/rho)
 *
 * The other physical data are those of the interior state.
 *
 * @param varSet     Euler variable set of the update variables
 * @param intState   interior state
 * @param primState  boundary values of (p, u, v[, w], T)
 * @param bndState   the boundary state in update variables
 */
inline void computeNSBoundaryState(Physics::NavierStokes::EulerVarSet& varSet,
                                   const Framework::State& intState,
                                   const RealVector& primState,
                                   RealVector& bndState)
{
  using Physics::NavierStokes::EulerTerm;

  RealVector physData;
  varSet.getModel()->resizePhysicalData(physData);
  varSet.computePhysicalData(intState,physData);

  const CFuint dim = primState.size() - 2;
  const CFreal p = primState[0];
  const CFreal T = primState[dim+1];
  const CFreal rho = p/(varSet.getModel()->getR()*T);
  const CFreal gamma = varSet.getModel()->getGamma();

  CFreal speed2 = 0.;
  for (CFuint iDim = 0; iDim < dim; ++iDim)
  {
    physData[EulerTerm::VX+iDim] = primState[1+iDim];
    speed2 += primState[1+iDim]*primState[1+iDim];
  }

  physData[EulerTerm::RHO] = rho;
  physData[EulerTerm::P] = p;
  physData[EulerTerm::T] = T;
  physData[EulerTerm::V] = std::sqrt(speed2);
  physData[EulerTerm::H] = gamma/(gamma-1.)*p/rho + 0.5*speed2;
  physData[EulerTerm::E] = physData[EulerTerm::H] - p/rho;
  physData[EulerTerm::A] = std::sqrt(gamma*p/rho);

  Framework::State state;
  varSet.computeStateFromPhysicalData(physData,state);
  bndState = *state.getData();
}

//////////////////////////////////////////////////////////////////////////////

/**
 * Prescribes the conductive heat flux at a wall on the temperature gradient:
 *
 *   n.grad(T) = -q/(C_Q lambda_eff)
 *
 * with -C_Q lambda_eff the heat flux of the variable set for a unit normal
 * temperature gradient. q is positive into the wall and n points out of the
 * fluid. For q = 0 the normal component of grad(T) is removed.
 *
 * @param varSet    diffusive variable set
 * @param state     boundary state
 * @param grads     boundary gradients of the gradient variables
 * @param normal    unit normal
 * @param tempID    index of T in the gradient variables
 * @param heatFlux  prescribed heat flux q
 */
inline void prescribeNSWallHeatFlux(Physics::NavierStokes::NavierStokesVarSet& varSet,
                                    const RealVector& state,
                                    std::vector< RealVector* >& grads,
                                    const RealVector& normal,
                                    const CFuint tempID,
                                    const CFreal heatFlux)
{
  if (heatFlux == 0.)
  {
    const CFreal normalGradT = MathTools::MathFunctions::innerProd(*grads[tempID],normal);
    *grads[tempID] -= normalGradT*normal;
    return;
  }

  // heat flux of a unit normal temperature gradient, -C_Q lambda_eff
  const RealVector gradT = *grads[tempID];
  *grads[tempID] = normal;
  const CFreal unitGradHeatFlux = varSet.getHeatFlux(state,grads,normal);
  *grads[tempID] = gradT;

  if (!(std::isfinite(unitGradHeatFlux) && unitGradHeatFlux < 0.))
  {
    throw Common::BadValueException(FromHere(),"Prescribed wall heat flux needs a finite positive effective conductivity.");
  }

  // replace the normal component of grad(T) by -q/(C_Q lambda_eff)
  const CFreal bndNormalGradT = heatFlux/unitGradHeatFlux;
  const CFreal normalGradT = MathTools::MathFunctions::innerProd(gradT,normal);
  *grads[tempID] += (bndNormalGradT-normalGradT)*normal;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionNavierStokes_NSBoundaryState_hh
