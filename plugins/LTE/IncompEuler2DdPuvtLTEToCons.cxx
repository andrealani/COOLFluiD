#include "LTE.hh"
#include "IncompEuler2DdPuvtLTEToCons.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/IncompEulerTerm.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<IncompEuler2DdPuvtLTEToCons, VarSetTransformer, LTEModule, 1> incompEuler2DdPuvtLTEToConsProvider("IncompEuler2DdPuvtLTEToCons");

//////////////////////////////////////////////////////////////////////////////

IncompEuler2DdPuvtLTEToCons::IncompEuler2DdPuvtLTEToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<IncompEulerTerm>()),
  _dhe(3)
{
  cf_assert(_model.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

IncompEuler2DdPuvtLTEToCons::~IncompEuler2DdPuvtLTEToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DdPuvtLTEToCons::transform(const State& state, 
					    State& result)
{
  const CFreal dp = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal T = state[3];

  const RealVector& refData = _model->getReferencePhysicalData();

  const CFreal p = _model->getThermodynamPressInf() + dp; // *getModel()->getPressRef();
  CFreal pdim = p; // *refData[IncompEulerTerm::P];
  CFreal Tdim = T*refData[IncompEulerTerm::T];

  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  // set the composition
  library->setComposition(Tdim,pdim);
  library->setDensityEnthalpyEnergy(Tdim,pdim,_dhe);
  const CFreal rho = _dhe[0]/refData[IncompEulerTerm::RHO];
  const CFreal V2 = u*u + v*v;

  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = rho*(_dhe[2]/refData[IncompEulerTerm::E] + 0.5*V2);
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DdPuvtLTEToCons::transformFromRef(const RealVector& data, 
						   State& result)
{
  const CFreal rho = data[IncompEulerTerm::RHO];

  result[0] = rho;
  result[1] = rho*data[IncompEulerTerm::VX];
  result[2] = rho*data[IncompEulerTerm::VY];
  result[3] = rho*data[IncompEulerTerm::E];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
