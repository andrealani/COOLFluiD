#include "MHD/MHD.hh"
#include "MHD3DProjectionPolytropicConsToPrim.hh"
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

Environment::ObjectProvider<MHD3DProjectionPolytropicConsToPrim, VarSetTransformer, MHDModule, 1> mhd3DProjectionPolytropicConsToPrimProvider("MHD3DProjectionPolytropicConsToPrim");

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionPolytropicConsToPrim::MHD3DProjectionPolytropicConsToPrim(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<MHDProjectionPolytropicTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionPolytropicConsToPrim::~MHD3DProjectionPolytropicConsToPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPolytropicConsToPrim::transform(const State& state, State& result)
{
  const CFreal rho = state[0];
  const CFreal u = state[1]/rho;
  const CFreal v = state[2]/rho;
  const CFreal w = state[3]/rho;
  const CFreal Bx = state[4];
  const CFreal By = state[5];
  const CFreal Bz = state[6];
  const CFreal phi = state[7];

  result[0] = rho;
  result[1] = u;
  result[2] = v;
  result[3] = w;
  result[4] = Bx;
  result[5] = By;
  result[6] = Bz;
  result[7] = phi;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPolytropicConsToPrim::transformFromRef(const RealVector& data, 
						 State& result) 
{
  result[0] = data[MHDProjectionPolytropicTerm::RHO];
  result[1] = data[MHDProjectionPolytropicTerm::VX];
  result[2] = data[MHDProjectionPolytropicTerm::VY];
  result[3] = data[MHDProjectionPolytropicTerm::VZ];
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
