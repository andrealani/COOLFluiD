#include "MHD/MHD.hh"
#include "MHD/MHD3DProjectionEpsConsToPrim.hh"
#include "MHD/MHDProjectionEpsTerm.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MHD3DProjectionEpsConsToPrim, VarSetTransformer, MHDModule, 1> mhd3DProjectionEpsConsToPrimProvider("MHD3DProjectionEpsConsToPrim");

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionEpsConsToPrim::MHD3DProjectionEpsConsToPrim(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<MHDProjectionEpsTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionEpsConsToPrim::~MHD3DProjectionEpsConsToPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionEpsConsToPrim::transform(const State& state, State& result)
{
  const CFreal gammaMinus1 = _model->getGamma() - 1.;
  const CFreal rho = state[0];
  const CFreal u = state[1]/rho;
  const CFreal v = state[2]/rho;
  const CFreal w = state[3]/rho;
  const CFreal Bx = state[4];
  const CFreal By = state[5];
  const CFreal Bz = state[6];
  const CFreal p = gammaMinus1 *
    (state[7] - 0.5*(rho*(u*u + v*v + w*w)
                       + Bx*Bx + By*By + Bz*Bz));
  const CFreal phi = state[8];

  result[0] = rho;
  result[1] = u;
  result[2] = v;
  result[3] = w;
  result[4] = Bx;
  result[5] = By;
  result[6] = Bz;
  result[7] = p;
  result[8] = phi;
  result[9] = state[9];
  result[10] = state[10];
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionEpsConsToPrim::transformFromRef(const RealVector& data, 
						 State& result) 
{
  result[0] = data[MHDProjectionTerm::RHO];
  result[1] = data[MHDProjectionTerm::VX];
  result[2] = data[MHDProjectionTerm::VY];
  result[3] = data[MHDProjectionTerm::VZ];
  result[4] = data[MHDProjectionTerm::BX];
  result[5] = data[MHDProjectionTerm::BY];
  result[6] = data[MHDProjectionTerm::BZ];
  result[7] = data[MHDProjectionTerm::P];
  result[8] = data[MHDProjectionTerm::PHI];
  result[9] = data[MHDProjectionEpsTerm::EPSP];
  result[10] = data[MHDProjectionEpsTerm::EPSM]; 
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
