#include "EntropyMHD/EntropyMHD.hh"
#include "EntropyMHD3DCons.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace EntropyMHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EntropyMHD3DCons, ConvectiveVarSet, EntropyMHDModule, 1>
entropyMHD3DProjectionConsProvider("EntropyMHD3DCons");

//////////////////////////////////////////////////////////////////////////////

EntropyMHD3DCons::EntropyMHD3DCons(Common::SafePtr<BaseTerm> term)
  : EntropyMHD3DVarSet(term)
{
  vector<std::string> names(9);
  names[0] = "rho";
  names[1] = "rhoU";
  names[2] = "rhoV";
  names[3] = "rhoW";
  names[4] = "Bx";
  names[5] = "By";
  names[6] = "Bz";
  names[7] = "sigma";
  names[8] = "phi";
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

EntropyMHD3DCons::~EntropyMHD3DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DCons::setup()
{
  EntropyMHD3DVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DCons::computePhysicalData(const State& state,
                                               RealVector& data)
{
  // Conservative state: (rho, rhoU, rhoV, rhoW, Bx, By, Bz, sigma, phi)
  // where sigma = p * rho^(1-gamma) is the entropy variable
  const CFreal rho   = state[0];
  const CFreal rhoU  = state[1];
  const CFreal rhoV  = state[2];
  const CFreal rhoW  = state[3];
  const CFreal Bx    = state[4];
  const CFreal By    = state[5];
  const CFreal Bz    = state[6];
  const CFreal sigma = state[7];
  const CFreal phi   = state[8];

  const CFreal rhoInv = 1.0 / rho;
  const CFreal u = rhoU * rhoInv;
  const CFreal v = rhoV * rhoInv;
  const CFreal w = rhoW * rhoInv;

  const CFreal gamma = getModel()->getGamma();

  // Recover pressure: p = sigma * rho^(gamma-1)
  const CFreal p = sigma * std::pow(rho, gamma - 1.0);

  data[EntropyMHDTerm::RHO]   = rho;
  data[EntropyMHDTerm::P]     = p;
  data[EntropyMHDTerm::VX]    = u;
  data[EntropyMHDTerm::VY]    = v;
  data[EntropyMHDTerm::VZ]    = w;
  data[EntropyMHDTerm::BX]    = Bx;
  data[EntropyMHDTerm::BY]    = By;
  data[EntropyMHDTerm::BZ]    = Bz;
  data[EntropyMHDTerm::GAMMA] = gamma;
  data[EntropyMHDTerm::PHI]   = phi;
  data[EntropyMHDTerm::T]     = sigma;  // sigma slot (named T for historical reasons)

  // Coordinates (read by EntropyMHDSourceTerm for gravity)
  const RealVector& coord = state.getCoordinates();
  data[EntropyMHDTerm::XP] = coord[XX];
  data[EntropyMHDTerm::YP] = coord[YY];
  data[EntropyMHDTerm::ZP] = coord[ZZ];
}

//////////////////////////////////////////////////////////////////////////////

void EntropyMHD3DCons::computeStateFromPhysicalData(const RealVector& data,
                                                        State& state)
{
  const CFreal rho = data[EntropyMHDTerm::RHO];

  state[0] = rho;
  state[1] = rho * data[EntropyMHDTerm::VX];
  state[2] = rho * data[EntropyMHDTerm::VY];
  state[3] = rho * data[EntropyMHDTerm::VZ];
  state[4] = data[EntropyMHDTerm::BX];
  state[5] = data[EntropyMHDTerm::BY];
  state[6] = data[EntropyMHDTerm::BZ];
  // sigma = p * rho^(1-gamma)
  state[7] = data[EntropyMHDTerm::P] * std::pow(rho, 1.0 - data[EntropyMHDTerm::GAMMA]);
  state[8] = data[EntropyMHDTerm::PHI];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace EntropyMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
