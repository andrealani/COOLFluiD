#include "MHD/MHD.hh"
#include "MHD2DProjectionConsToPrim.hh"
#include "MHDProjectionTerm.hh"
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

Environment::ObjectProvider<MHD2DProjectionConsToPrim, VarSetTransformer, MHDModule, 1> mhd2DProjectionConsToPrimProvider("MHD2DProjectionConsToPrim");

//////////////////////////////////////////////////////////////////////////////

MHD2DProjectionConsToPrim::MHD2DProjectionConsToPrim(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<MHDProjectionTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD2DProjectionConsToPrim::~MHD2DProjectionConsToPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DProjectionConsToPrim::transform(const State& state, 
					  State& result)
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
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DProjectionConsToPrim::transformFromRef(const RealVector& data, State& result)
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
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
