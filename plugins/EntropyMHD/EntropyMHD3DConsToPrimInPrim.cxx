#include "EntropyMHD/EntropyMHD.hh"
#include "EntropyMHD3DConsToPrimInPrim.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace EntropyMHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EntropyMHD3DConsToPrimInPrim, VarSetMatrixTransformer, EntropyMHDModule, 1>
entropyMHD3DConsToPrimInPrimProvider("EntropyMHD3DConsToPrimInPrim");

//////////////////////////////////////////////////////////////////////////////

EntropyMHD3DConsToPrimInPrim::EntropyMHD3DConsToPrimInPrim(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model)
{
}

//////////////////////////////////////////////////////////////////////////////

EntropyMHD3DConsToPrimInPrim::~EntropyMHD3DConsToPrimInPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DConsToPrimInPrim::setMatrix(const RealVector& state)
{
  // state is in Prim variables: (rho, u, v, w, Bx, By, Bz, sigma, phi)
  //
  // dW/dU where W=Prim, U=Cons:
  //   W[0] = rho            -> row 0: (1, 0, 0, 0, ...)
  //   W[1] = rho*u/rho = u  -> row 1: (-u/rho, 1/rho, 0, 0, ...)
  //   W[2] = rho*v/rho = v  -> row 2: (-v/rho, 0, 1/rho, 0, ...)
  //   W[3] = rho*w/rho = w  -> row 3: (-w/rho, 0, 0, 1/rho, ...)
  //   W[4..6] = B           -> identity
  //   W[7] = sigma          -> identity (sigma same in both)
  //   W[8] = phi            -> identity

  const CFreal ovrbar = 1.0 / state[0];
  const CFreal ubar = state[1];
  const CFreal vbar = state[2];
  const CFreal wbar = state[3];

  _transMatrix(0,0) = 1.0;
  _transMatrix(1,0) = -ubar * ovrbar;
  _transMatrix(1,1) = ovrbar;
  _transMatrix(2,0) = -vbar * ovrbar;
  _transMatrix(2,2) = ovrbar;
  _transMatrix(3,0) = -wbar * ovrbar;
  _transMatrix(3,3) = ovrbar;
  _transMatrix(4,4) = 1.0;
  _transMatrix(5,5) = 1.0;
  _transMatrix(6,6) = 1.0;
  _transMatrix(7,7) = 1.0;
  _transMatrix(8,8) = 1.0;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace EntropyMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
