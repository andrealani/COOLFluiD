// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"

#include "Common/BadValueException.hh"

#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/PositivityIDPNS.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PositivityIDPNS, FluxReconstructionSolverData, FluxReconstructionNavierStokesModule>
    PositivityIDPNSProvider("PositivityIDPNS");

//////////////////////////////////////////////////////////////////////////////

void PositivityIDPNS::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

PositivityIDPNS::PositivityIDPNS(const std::string& name) :
  BasePositivityIDP(name),
  m_gammaMinusOne(0.),
  m_allIndices(),
  m_iE(0),
  m_consUpdateVars(true)
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

PositivityIDPNS::~PositivityIDPNS()
{
}

//////////////////////////////////////////////////////////////////////////////

void PositivityIDPNS::constraintsAtPoint(const RealVector& cons,
                                         CFreal& rho, CFreal& p, CFreal& B2) const
{
  rho = cons[0];

  // no magnetic field in this physics
  B2 = 0.;

  // momentum carries one component per spatial dimension
  CFreal m2 = 0.;
  for (CFuint i = 1; i <= m_dim; ++i)
  {
    m2 += cons[i]*cons[i];
  }

  // p = (gamma-1) (E - |m|^2/(2 rho))
  //
  // With rho <= 0 the -|m|^2/(2 rho) term flips sign and p comes out
  // spuriously large and positive, which would silently disarm the pressure
  // constraint. Return a strictly negative value instead so the caller sees
  // the point as inadmissible; the sequential path in the base class handles
  // it by limiting density first.
  if (rho <= 0.)
  {
    p = -1.;
    return;
  }

  p = m_gammaMinusOne*(cons[m_iE] - 0.5*m2/rho);
}

//////////////////////////////////////////////////////////////////////////////

void PositivityIDPNS::consToUpdate(const RealVector& cons, RealVector& update) const
{
  if (m_consUpdateVars)
  {
    for (CFuint i = 0; i < m_nbrEqs; ++i)
    {
      update[i] = cons[i];
    }
    return;
  }

  // Prim: (rho, u, v, [w], p)
  const CFreal rho = cons[0];

  update[0] = rho;

  CFreal m2 = 0.;
  for (CFuint i = 1; i <= m_dim; ++i)
  {
    update[i] = cons[i]/rho;
    m2 += cons[i]*cons[i];
  }

  update[m_iE] = m_gammaMinusOne*(cons[m_iE] - 0.5*m2/rho);
}

//////////////////////////////////////////////////////////////////////////////

const std::vector< CFuint >& PositivityIDPNS::scaledIndices(const ScaleMode mode) const
{
  // every component is hydrodynamic here, so both modes scale all of them
  return m_allIndices;
}

//////////////////////////////////////////////////////////////////////////////

void PositivityIDPNS::setupPhysics()
{
  CFAUTOTRACE;

  // The state layout below is the single species conservative one. Temperature
  // based update sets put temperature where this expects energy or pressure,
  // so refuse rather than constrain the wrong quantity.
  const std::string updateVarStr = getMethodData().getUpdateVarStr();
  if (updateVarStr == "Cons")
  {
    m_consUpdateVars = true;
  }
  else if (updateVarStr == "Prim")
  {
    m_consUpdateVars = false;
  }
  else
  {
    throw BadValueException (FromHere(),
      "PositivityIDPNS requires UpdateVar = Cons or Prim, got " + updateVarStr);
  }

  SafePtr< EulerTerm > eulerTerm = PhysicalModelStack::getActive()->getImplementor()
                                     ->getConvectiveTerm().d_castTo< EulerTerm >();
  m_gammaMinusOne = eulerTerm->getGamma() - 1.0;

  // Layout: (rho, rhoU, rhoV, [rhoW], rhoE), so energy sits right after the
  // momentum components. Extra transported scalars, from turbulence or from a
  // species model, would shift that slot and are not handled here.
  m_iE = m_dim + 1;

  if (m_nbrEqs != m_dim + 2)
  {
    throw BadValueException (FromHere(),
      "PositivityIDPNS expects the single species Euler layout with dim+2 equations. "
      "Models carrying extra transported scalars need their own subclass.");
  }

  m_allIndices.clear();
  for (CFuint i = 0; i < m_nbrEqs; ++i)
  {
    m_allIndices.push_back(i);
  }

  CFLog(NOTICE, "PositivityIDPNS: gamma = " << eulerTerm->getGamma()
        << ", dim = " << m_dim
        << ", nbEqs = " << m_nbrEqs
        << " (m at 1.." << m_dim
        << ", E at " << m_iE
        << ", update vars " << updateVarStr << ")\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
