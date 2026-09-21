// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MultiScalarTerm.hh"

#include "Common/BadValueException.hh"
#include "Common/NotImplementedException.hh"

#include "MathTools/MathConsts.hh"

#include "NavierStokes/EulerTerm.hh"

#include "FluxReconstructionNEQ/FluxReconstructionNEQ.hh"
#include "FluxReconstructionNEQ/PositivityIDPNEQ.hh"

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

MethodCommandProvider<PositivityIDPNEQ, FluxReconstructionSolverData, FluxReconstructionNEQModule>
    PositivityIDPNEQProvider("PositivityIDPNEQ");

//////////////////////////////////////////////////////////////////////////////

void PositivityIDPNEQ::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

PositivityIDPNEQ::PositivityIDPNEQ(const std::string& name) :
  BasePositivityIDP(name),
  m_nbSpecies(0),
  m_speciesIndices()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

PositivityIDPNEQ::~PositivityIDPNEQ()
{
}

//////////////////////////////////////////////////////////////////////////////

void PositivityIDPNEQ::constraintsAtPoint(const RealVector& cons,
                                          CFreal& rho, CFreal& p, CFreal& B2) const
{
  // the binding density constraint is the smallest partial density: because
  // the limiter is a convex blend towards the cell mean, one theta driven by
  // the minimum makes every species positive at once.
  rho = cons[0];
  for (CFuint is = 1; is < m_nbSpecies; ++is)
  {
    rho = std::min(rho, cons[is]);
  }

  // no thermal constraint, see the class comment. Report a value that can
  // never bind, and the pressure constraint is switched off in setupPhysics()
  // as well, so this is belt and braces.
  p  = MathTools::MathConsts::CFrealMax();
  B2 = 0.;
}

//////////////////////////////////////////////////////////////////////////////

void PositivityIDPNEQ::consToUpdate(const RealVector& cons, RealVector& update) const
{
  // Unreachable: writeBackStates() is overridden so that only the partial
  // densities, which are shared between the two variable sets, are written
  // back. A general conservative -> RhoivtTv mapping needs a Newton inversion
  // for T and Tv against the chemistry library and does not exist in this
  // codebase. Fail loudly rather than return a wrong state if the base class
  // ever starts calling this.
  throw NotImplementedException (FromHere(),
    "PositivityIDPNEQ::consToUpdate(): conservative to RhoivtTv is not implemented. "
    "This limiter scales only the partial densities and writes them back directly.");
}

//////////////////////////////////////////////////////////////////////////////

const std::vector< CFuint >& PositivityIDPNEQ::scaledIndices(const ScaleMode mode) const
{
  // only the partial densities are scaled, in either mode. There is no
  // magnetic field here, so the Hydro/Full distinction that exists to protect
  // div B has nothing to separate.
  return m_speciesIndices;
}

//////////////////////////////////////////////////////////////////////////////

void PositivityIDPNEQ::writeBackStates()
{
  // the first m_nbSpecies entries of RhoivtTv and of the conservative set are
  // the same partial densities, so the scaled values go straight back with no
  // transformation. u, v, T and Tv were never scaled and are left untouched.
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    for (CFuint is = 0; is < m_nbSpecies; ++is)
    {
      (*((*m_cellStates)[iSol]))[is] = m_consSol[iSol][is];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void PositivityIDPNEQ::setupPhysics()
{
  CFAUTOTRACE;

  // The write back above copies the scaled conservative densities straight
  // into the stored state, which is only correct when the update set starts
  // with the same partial densities.
  const std::string updateVarStr = getMethodData().getUpdateVarStr();
  if (updateVarStr != "RhoivtTv" && updateVarStr != "Rhoivt")
  {
    throw BadValueException (FromHere(),
      "PositivityIDPNEQ requires UpdateVar = RhoivtTv or Rhoivt, got " + updateVarStr);
  }

  SafePtr< MultiScalarTerm< EulerTerm > > term =
    PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm().d_castTo< MultiScalarTerm< EulerTerm > >();

  m_nbSpecies = term->getNbScalarVars(0);

  if (m_nbSpecies == 0 || m_nbSpecies > m_nbrEqs)
  {
    throw BadValueException (FromHere(),
      "PositivityIDPNEQ: inconsistent number of species.");
  }

  m_speciesIndices.clear();
  for (CFuint is = 0; is < m_nbSpecies; ++is)
  {
    m_speciesIndices.push_back(is);
  }

  // the thermal constraint would need sensible internal energy, which needs
  // species formation energies that PLATO does not provide. Leaving it on
  // would silently do nothing, so make the limitation explicit.
  if (m_enablePressure)
  {
    CFLog(NOTICE, "PositivityIDPNEQ: switching the pressure constraint off, it is not supported for NEQ.\n");
    m_enablePressure = false;
  }

  CFLog(NOTICE, "PositivityIDPNEQ: nbSpecies = " << m_nbSpecies
        << ", nbEqs = " << m_nbrEqs
        << ", dim = " << m_dim
        << ", scaling the partial densities only\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
