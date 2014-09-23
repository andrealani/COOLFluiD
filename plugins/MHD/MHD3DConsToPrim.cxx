#include "MHD/MHD.hh"
#include "MHD3DConsToPrim.hh"
#include "MHDPhysicalModel.hh"
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

Environment::ObjectProvider<MHD3DConsToPrim, VarSetTransformer, MHDModule, 1> mhd3DConsToPrimToPrimCProvider("MHD3DConsToPrim");

//////////////////////////////////////////////////////////////////////////////

MHD3DConsToPrim::MHD3DConsToPrim(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<MHDTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DConsToPrim::~MHD3DConsToPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DConsToPrim::transform(const State& state, State& result)
{
  const CFreal rho = state[0];
  const CFreal u   = state[1]/rho;
  const CFreal v   = state[2]/rho;
  const CFreal w   = state[3]/rho;
  const CFreal Bx  = state[4];
  const CFreal By  = state[5];
  const CFreal Bz  = state[6];
  const CFreal gammaMinus1 = _model->getGamma() - 1.0;
  const CFreal p   = gammaMinus1*
    (state[7] - 0.5*(rho*(u*u + v*v + w*w)
				    + Bx*Bx + By*By + Bz*Bz))  ;

  result[0] = rho ;
  result[1] = u ;
  result[2] = v ;
  result[3] = w ;
  result[4] = Bx ;
  result[5] = By ;
  result[6] = Bz ;
  result[7] = p ;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DConsToPrim::transformFromRef(const RealVector& data, State& result)
{
  result[0] = data[MHDTerm::RHO];
  result[1] = data[MHDTerm::VX];
  result[2] = data[MHDTerm::VY];
  result[3] = data[MHDTerm::VZ];
  result[4] = data[MHDTerm::BX];
  result[5] = data[MHDTerm::BY];
  result[6] = data[MHDTerm::BZ];
  result[7] = data[MHDTerm::P];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
