#include "EntropyMHD/EntropyMHD.hh"
#include "EntropyMHD3DConsToPrim.hh"
#include "EntropyMHDTerm.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace EntropyMHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EntropyMHD3DConsToPrim, VarSetTransformer, EntropyMHDModule, 1>
entropyMHD3DConsToPrimProvider("EntropyMHD3DConsToPrim");

//////////////////////////////////////////////////////////////////////////////

EntropyMHD3DConsToPrim::EntropyMHD3DConsToPrim(Common::SafePtr<PhysicalModelImpl> model)
  : VarSetTransformer(model)
{
}

//////////////////////////////////////////////////////////////////////////////

EntropyMHD3DConsToPrim::~EntropyMHD3DConsToPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DConsToPrim::transform(const State& state, State& result)
{
  // Cons = (rho, rhoU, rhoV, rhoW, Bx, By, Bz, sigma, phi)
  // Prim = (rho, u, v, w, Bx, By, Bz, sigma, phi)
  // Only velocities change (momentum -> velocity). sigma is the same.
  const CFreal rho = state[0];

  result[0] = rho;
  result[1] = state[1] / rho;  // rhoU -> u
  result[2] = state[2] / rho;  // rhoV -> v
  result[3] = state[3] / rho;  // rhoW -> w
  result[4] = state[4];        // Bx
  result[5] = state[5];        // By
  result[6] = state[6];        // Bz
  result[7] = state[7];        // sigma (unchanged)
  result[8] = state[8];        // phi
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DConsToPrim::transformFromRef(const RealVector& data, State& result)
{
  // From pdata to primitive state: sigma = p * rho^(1-gamma)
  result[0] = data[EntropyMHDTerm::RHO];
  result[1] = data[EntropyMHDTerm::VX];
  result[2] = data[EntropyMHDTerm::VY];
  result[3] = data[EntropyMHDTerm::VZ];
  result[4] = data[EntropyMHDTerm::BX];
  result[5] = data[EntropyMHDTerm::BY];
  result[6] = data[EntropyMHDTerm::BZ];
  result[7] = data[EntropyMHDTerm::T];  // sigma stored in T slot
  result[8] = data[EntropyMHDTerm::PHI];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace EntropyMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
