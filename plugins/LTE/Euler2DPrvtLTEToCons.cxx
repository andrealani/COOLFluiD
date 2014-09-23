#include "LTE.hh"
#include "Euler2DPrvtLTEToCons.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/EulerTerm.hh"
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

Environment::ObjectProvider<Euler2DPrvtLTEToCons, VarSetTransformer, LTEModule, 1> euler2DPrvtLTEToConsProvider("Euler2DPrvtLTEToCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DPrvtLTEToCons::Euler2DPrvtLTEToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>()),
  _dhe(3)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPrvtLTEToCons::~Euler2DPrvtLTEToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrvtLTEToCons::transform(const State& state, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const RealVector& refData = _model->getReferencePhysicalData();
  
  CFreal p = _model->getPressureFromState(state[0]);
  CFreal T = state[3];
  CFreal pdim = p*refData[EulerTerm::P];
  CFreal Tdim = T*_model->getTempRef();
  
  // set the composition
  library->setComposition(Tdim,pdim);
  library->setDensityEnthalpyEnergy(Tdim,pdim,_dhe);
  
  const CFreal rho = _dhe[0]/refData[EulerTerm::RHO];
  const CFreal rhou = state[1];
  const CFreal rhov = state[2];
  const CFreal V2 = (rhou*rhou + rhov*rhov)/(rho*rho);
  
  result[0] = rho;
  result[1] = rhou;
  result[2] = rhov;
  result[3] = rho*(_dhe[2]/refData[EulerTerm::H] + 0.5*V2);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPrvtLTEToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFreal rho = data[EulerTerm::RHO];

  result[0] = rho;
  result[1] = rho*data[EulerTerm::VX];
  result[2] = rho*data[EulerTerm::VY];
  result[3] = rho*data[EulerTerm::E];
  //rho*data.avH - data.avP; gives slightly different result
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LTE

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
