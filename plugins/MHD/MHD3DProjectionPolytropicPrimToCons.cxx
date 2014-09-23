#include "MHD/MHD.hh"
#include "MHD3DProjectionPolytropicPrimToCons.hh"
#include "MHDProjectionPolytropicTerm.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MHD3DProjectionPolytropicPrimToCons, VarSetTransformer, MHDModule, 1> mhd3DProjectionPolytropicPrimToConsProvider("MHD3DProjectionPolytropicPrimToCons");

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionPolytropicPrimToCons::MHD3DProjectionPolytropicPrimToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<MHDProjectionPolytropicTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionPolytropicPrimToCons::~MHD3DProjectionPolytropicPrimToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPolytropicPrimToCons::transform(const State& state, State& result)
{
  const CFreal rho = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal Bx = state[4];
  const CFreal By = state[5];
  const CFreal Bz = state[6];
  const CFreal phi = state[7];

  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = rho*w;
  result[4] = Bx;
  result[5] = By;
  result[6] = Bz;
  result[7] = phi;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPolytropicPrimToCons::transformFromRef(const RealVector& data, 
						 State& result) 
{
  result[0] = data[MHDProjectionPolytropicTerm::RHO];
  result[1] = data[MHDProjectionPolytropicTerm::RHO]*data[MHDProjectionPolytropicTerm::VX];
  result[2] = data[MHDProjectionPolytropicTerm::RHO]*data[MHDProjectionPolytropicTerm::VY];
  result[3] = data[MHDProjectionPolytropicTerm::RHO]*data[MHDProjectionPolytropicTerm::VZ];
  result[4] = data[MHDProjectionPolytropicTerm::BX];
  result[5] = data[MHDProjectionPolytropicTerm::BY];
  result[6] = data[MHDProjectionPolytropicTerm::BZ];
  result[7] = data[MHDProjectionPolytropicTerm::PHI];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
