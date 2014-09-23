#include "LTE.hh"
#include "Euler2DPuvtLTEToCons.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "NavierStokes/EulerTerm.hh"

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

Environment::ObjectProvider<Euler2DPuvtLTEToCons, VarSetTransformer, LTEModule, 1> euler2DPuvtLTEToConsProvider("Euler2DPuvtLTEToCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtLTEToCons::Euler2DPuvtLTEToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>()),
  _dhe(3)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtLTEToCons::~Euler2DPuvtLTEToCons()
{
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtLTEToCons::transform(const State& state, State& result)
{
  const RealVector& refData = _model->getReferencePhysicalData();

  const CFreal u = state[1];
  const CFreal v = state[2];
  CFreal p = _model->getPressureFromState(state[0]);
  CFreal T = state[3];
  CFreal pdim = p*refData[EulerTerm::P];
  CFreal Tdim = T*_model->getTempRef();

  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  // set the composition
  library->setComposition(Tdim,pdim);
  library->setDensityEnthalpyEnergy(Tdim,pdim,_dhe);
  const CFreal rho = _dhe[0]/refData[EulerTerm::RHO];
  const CFreal V2 = u*u + v*v;

  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = rho*(_dhe[2]/refData[EulerTerm::H] + 0.5*V2);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtLTEToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFreal rho = data[EulerTerm::RHO];

  result[0] = rho;
  result[1] = rho*data[EulerTerm::VX];
  result[2] = rho*data[EulerTerm::VY];
  result[3] = rho*data[EulerTerm::E];
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LTE

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
