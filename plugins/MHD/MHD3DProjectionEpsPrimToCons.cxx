#include "MHD/MHD.hh"
#include "MHD/MHD3DProjectionEpsPrimToCons.hh"
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

Environment::ObjectProvider<MHD3DProjectionEpsPrimToCons, VarSetTransformer, MHDModule, 1> mhd3DProjectionEpsPrimToConsProvider("MHD3DProjectionEpsPrimToCons");

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionEpsPrimToCons::MHD3DProjectionEpsPrimToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<MHDProjectionEpsTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionEpsPrimToCons::~MHD3DProjectionEpsPrimToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionEpsPrimToCons::transform(const State& state, State& result)
{
  const CFreal gammaMinus1 = _model->getGamma() - 1.;
  const CFreal rho = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal Bx = state[4];
  const CFreal By = state[5];
  const CFreal Bz = state[6];
  const CFreal p = state[7];
  const CFreal phi = state[8];

  const CFreal rhoV2 = rho*(u*u + v*v + w*w);
  const CFreal B2 = (Bx*Bx + By*By + Bz*Bz);

  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = rho*w;
  result[4] = Bx;
  result[5] = By;
  result[6] = Bz;
  result[7] = (p/gammaMinus1) + 0.5*(rhoV2+B2);
  result[8] = phi;
  result[9] = state[9];
  result[10] = state[10];
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionEpsPrimToCons::transformFromRef(const RealVector& data, 
						 State& result) 
{
  const CFreal gammaMinus1 = _model->getGamma() - 1.;

  const CFreal rhoV2 = data[MHDProjectionTerm::RHO]*(data[MHDProjectionTerm::VX]*data[MHDProjectionTerm::VX] + 
	data[MHDProjectionTerm::VY]*data[MHDProjectionTerm::VY] + 
        data[MHDProjectionTerm::VZ]*data[MHDProjectionTerm::VZ]);

  const CFreal B2 = data[MHDProjectionTerm::BX]*data[MHDProjectionTerm::BX] + 
        data[MHDProjectionTerm::BY]*data[MHDProjectionTerm::BY] + 
        data[MHDProjectionTerm::BZ]*data[MHDProjectionTerm::BZ];

  result[0] = data[MHDProjectionTerm::RHO];
  result[1] = data[MHDProjectionTerm::RHO]*data[MHDProjectionTerm::VX];
  result[2] = data[MHDProjectionTerm::RHO]*data[MHDProjectionTerm::VY];
  result[3] = data[MHDProjectionTerm::RHO]*data[MHDProjectionTerm::VZ];
  result[4] = data[MHDProjectionTerm::BX];
  result[5] = data[MHDProjectionTerm::BY];
  result[6] = data[MHDProjectionTerm::BZ];
  result[7] = (data[MHDProjectionTerm::P]/gammaMinus1) + 0.5*(rhoV2+B2);
  result[8] = data[MHDProjectionTerm::PHI];
  result[9] = data[MHDProjectionEpsTerm::EPSP];
  result[10] = data[MHDProjectionEpsTerm::EPSM];
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
