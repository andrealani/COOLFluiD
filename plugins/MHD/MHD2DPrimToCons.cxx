#include "MHD/MHD.hh"
#include "MHD2DPrimToCons.hh"
#include "MHDTerm.hh"
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

Environment::ObjectProvider<MHD2DPrimToCons, VarSetTransformer, MHDModule, 1> 
mhd2DPrimToConsProvider("MHD2DPrimToCons");

//////////////////////////////////////////////////////////////////////////////

MHD2DPrimToCons::MHD2DPrimToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<MHDTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD2DPrimToCons::~MHD2DPrimToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DPrimToCons::transform(const State& state, State& result)
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
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DPrimToCons::transformFromRef(const RealVector& data, 
						 State& result) 
{
  const CFreal gammaMinus1 = _model->getGamma() - 1.;

  const CFreal rhoV2 = data[MHDTerm::RHO]*(data[MHDTerm::VX]*data[MHDTerm::VX] + 
	data[MHDTerm::VY]*data[MHDTerm::VY] + 
        data[MHDTerm::VZ]*data[MHDTerm::VZ]);

  const CFreal B2 = data[MHDTerm::BX]*data[MHDTerm::BX] + 
        data[MHDTerm::BY]*data[MHDTerm::BY] + 
        data[MHDTerm::BZ]*data[MHDTerm::BZ];

  result[0] = data[MHDTerm::RHO];
  result[1] = data[MHDTerm::RHO]*data[MHDTerm::VX];
  result[2] = data[MHDTerm::RHO]*data[MHDTerm::VY];
  result[3] = data[MHDTerm::RHO]*data[MHDTerm::VZ];
  result[4] = data[MHDTerm::BX];
  result[5] = data[MHDTerm::BY];
  result[6] = data[MHDTerm::BZ];
  result[7] = (data[MHDTerm::P]/gammaMinus1) + 0.5*(rhoV2+B2);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
