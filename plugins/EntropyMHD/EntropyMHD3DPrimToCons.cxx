#include "EntropyMHD/EntropyMHD.hh"
#include "EntropyMHD3DPrimToCons.hh"
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

Environment::ObjectProvider<EntropyMHD3DPrimToCons, VarSetTransformer, EntropyMHDModule, 1>
entropyMHD3DPrimToConsProvider("EntropyMHD3DPrimToCons");

//////////////////////////////////////////////////////////////////////////////

EntropyMHD3DPrimToCons::EntropyMHD3DPrimToCons(Common::SafePtr<PhysicalModelImpl> model)
  : VarSetTransformer(model)
{
}

//////////////////////////////////////////////////////////////////////////////

EntropyMHD3DPrimToCons::~EntropyMHD3DPrimToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DPrimToCons::transform(const State& state, State& result)
{
  // Prim = (rho, u, v, w, Bx, By, Bz, sigma, phi)
  // Cons = (rho, rhoU, rhoV, rhoW, Bx, By, Bz, sigma, phi)
  // Only velocities change (velocity -> momentum). sigma is the same.
  const CFreal rho = state[0];

  result[0] = rho;
  result[1] = rho * state[1];  // u -> rhoU
  result[2] = rho * state[2];  // v -> rhoV
  result[3] = rho * state[3];  // w -> rhoW
  result[4] = state[4];        // Bx
  result[5] = state[5];        // By
  result[6] = state[6];        // Bz
  result[7] = state[7];        // sigma (unchanged)
  result[8] = state[8];        // phi
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DPrimToCons::transformFromRef(const RealVector& data, State& result)
{
  // From pdata to conservative state: sigma stored in T slot
  const CFreal rho = data[EntropyMHDTerm::RHO];

  result[0] = rho;
  result[1] = rho * data[EntropyMHDTerm::VX];
  result[2] = rho * data[EntropyMHDTerm::VY];
  result[3] = rho * data[EntropyMHDTerm::VZ];
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
