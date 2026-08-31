#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"

#include "Common/BadValueException.hh"

#include "MHD/MHD3DProjectionVarSet.hh"
#include "MHD/MHD2DProjectionVarSet.hh"
#include "MHD/MHDTerm.hh"

#include "FluxReconstructionMHD/FluxReconstructionMHD.hh"
#include "FluxReconstructionMHD/PositivityIDPMHD.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<PositivityIDPMHD, FluxReconstructionSolverData, FluxReconstructionMHDModule>
    PositivityIDPMHDProvider("PositivityIDPMHD");

//////////////////////////////////////////////////////////////////////////////

void PositivityIDPMHD::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

PositivityIDPMHD::PositivityIDPMHD(const std::string& name) :
  BasePositivityIDP(name),
  m_gammaMinusOne(0.),
  m_hydroIndices(),
  m_fullIndices(),
  m_iE(0),
  m_iB(0)
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

PositivityIDPMHD::~PositivityIDPMHD()
{
}

//////////////////////////////////////////////////////////////////////////////

void PositivityIDPMHD::constraintsAtPoint(const RealVector& cons,
                                          CFreal& rho, CFreal& p, CFreal& B2) const
{
  rho = cons[0];

  // momentum and B always carry 3 components, also in 2D
  CFreal m2 = 0.;
  for (CFuint i = 1; i < 4; ++i)
  {
    m2 += cons[i]*cons[i];
  }

  B2 = 0.;
  for (CFuint i = m_iB; i < m_iB + 3; ++i)
  {
    B2 += cons[i]*cons[i];
  }

  // p = (gamma-1) (E - |m|^2/(2 rho) - |B|^2/2)
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

  p = m_gammaMinusOne*(cons[m_iE] - 0.5*m2/rho - 0.5*B2);
}

//////////////////////////////////////////////////////////////////////////////

void PositivityIDPMHD::consToUpdate(const RealVector& cons, RealVector& update) const
{
  const CFreal rho = cons[0];

  update[0] = rho;

  // momentum and B always carry 3 components, also in 2D
  CFreal m2 = 0.;
  for (CFuint i = 1; i < 4; ++i)
  {
    update[i] = cons[i]/rho;
    m2 += cons[i]*cons[i];
  }

  CFreal B2 = 0.;
  for (CFuint i = m_iB; i < m_iB + 3; ++i)
  {
    update[i] = cons[i];
    B2 += cons[i]*cons[i];
  }

  update[m_iE] = m_gammaMinusOne*(cons[m_iE] - 0.5*m2/rho - 0.5*B2);

  // phi, the GLM cleaning variable, is the same in both variable sets
  update[m_nbrEqs-1] = cons[m_nbrEqs-1];
}

//////////////////////////////////////////////////////////////////////////////

const std::vector< CFuint >& PositivityIDPMHD::scaledIndices(const ScaleMode mode) const
{
  return (mode == HYDRO) ? m_hydroIndices : m_fullIndices;
}

//////////////////////////////////////////////////////////////////////////////

void PositivityIDPMHD::setupPhysics()
{
  CFAUTOTRACE;

  // The state layout below assumes primitive update variables. With Cons the
  // index of the pressure slot would mean total energy instead, so refuse
  // rather than silently constrain the wrong quantity.
  const std::string updateVarStr = getMethodData().getUpdateVarStr();
  if (updateVarStr != "Prim")
  {
    throw BadValueException (FromHere(),
      "PositivityIDPMHD requires UpdateVar = Prim, got " + updateVarStr);
  }

  // d_castTo throws on failure rather than returning null, so branch on the
  // dimension first and only attempt the cast that can succeed.
  CFreal gamma = 0.;
  std::string potentialBType;

  if (m_dim == 2)
  {
    SafePtr< MHD2DProjectionVarSet > varSet2D =
      getMethodData().getUpdateVar().d_castTo< MHD2DProjectionVarSet >();
    gamma          = varSet2D->getModel()->getGamma();
    potentialBType = varSet2D->getModel()->getPotentialBType();
  }
  else
  {
    SafePtr< MHD3DProjectionVarSet > varSet3D =
      getMethodData().getUpdateVar().d_castTo< MHD3DProjectionVarSet >();
    gamma          = varSet3D->getModel()->getGamma();
    potentialBType = varSet3D->getModel()->getPotentialBType();
  }

  m_gammaMinusOne = gamma - 1.0;

  // The B0/B1 split is not used in this workflow
  // B in the state carries the whole field. If it were ever switched on, the
  // wave speeds would use B_total while the limiter constrains B1 alone, so
  // the magnetic terms here would be wrong. Refuse instead of being subtly off.
  if (potentialBType != "")
  {
    throw BadValueException (FromHere(),
      "PositivityIDPMHD does not support the B0/B1 potential field split (PotentialBType = "
      + potentialBType + "). The limiter assumes B in the state carries the full field.");
  }

  // Layout is dimension independent: (rho, rhoU, rhoV, rhoW, Bx, By, Bz, rhoE, phi).
  // The 2D projection model keeps the z components of velocity and B even
  // though they are typically zero, so it also carries 9 equations. Indexing
  // off m_dim would put B at 3 in 2D and silently treat rhoW as Bx.
  m_iB = 4;
  m_iE = 7;

  if (m_nbrEqs != 9)
  {
    throw BadValueException (FromHere(),
      "PositivityIDPMHD expects the 9 equation MHD projection layout (rho, rhoU, rhoV, rhoW, Bx, By, Bz, rhoE, phi), in 2D as well as 3D.");
  }

  // Hydro mode scales the hydrodynamic components and freezes B and phi, so
  // the discrete div B is untouched. Full mode scales everything, which always
  // has an admissible base state but perturbs div B.
  m_hydroIndices.clear();
  m_hydroIndices.push_back(0);
  for (CFuint i = 1; i < 4; ++i)
  {
    m_hydroIndices.push_back(i);
  }
  m_hydroIndices.push_back(m_iE);

  m_fullIndices.clear();
  for (CFuint i = 0; i < m_nbrEqs; ++i)
  {
    m_fullIndices.push_back(i);
  }

  CFLog(NOTICE, "PositivityIDPMHD: gamma = " << gamma
        << ", dim = " << m_dim
        << ", nbEqs = " << m_nbrEqs
        << " (m at 1..3, B at " << m_iB << ".." << m_iB+2
        << ", E at " << m_iE
        << ", phi at " << m_nbrEqs-1 << ")\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
