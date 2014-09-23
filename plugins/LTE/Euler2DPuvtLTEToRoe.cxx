#include "LTE.hh"
#include "Euler2DPuvtLTEToRoe.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Environment/ObjectProvider.hh"
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

Environment::ObjectProvider<Euler2DPuvtLTEToRoe, VarSetTransformer, LTEModule, 1> euler2DPuvtLTEToRoeProvider("Euler2DPuvtLTEToRoe");

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtLTEToRoe::Euler2DPuvtLTEToRoe(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>()),
  _dhe(3)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtLTEToRoe::~Euler2DPuvtLTEToRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtLTEToRoe::transform(const State& state, State& result)
{
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
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal h = _dhe[1]/refData[EulerTerm::H];
  const CFreal sqRho = sqrt(rho);
  
  result[0] = sqRho;
  result[1] = sqRho*u;
  result[2] = sqRho*v;
  result[3] = sqRho*(h + 0.5 * (u*u + v*v));
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtLTEToRoe::transformFromRef(const RealVector& data, State& result)
{
  const CFreal sqRho = sqrt(data[EulerTerm::RHO]);
  
  result[0] = sqRho;
  result[1] = sqRho*data[EulerTerm::VX];
  result[2] = sqRho*data[EulerTerm::VY];
  result[3] = sqRho*data[EulerTerm::H];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
