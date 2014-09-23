#include "LTE.hh"
#include "Euler3DLTEDemixPvtToCons.hh"
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

Environment::ObjectProvider<Euler3DLTEDemixPvtToCons, VarSetTransformer, LTEModule,1> euler3DLTEDemixPvtToConsProvider("Euler3DLTEDemixPvtToCons");

//////////////////////////////////////////////////////////////////////////////

Euler3DLTEDemixPvtToCons::Euler3DLTEDemixPvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerLTEDemixTerm>()),
  _dhe(3),
  _ye()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler3DLTEDemixPvtToCons::~Euler3DLTEDemixPvtToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixPvtToCons::transform(const State& state, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  const CFuint nbElements = _model->getNbScalarVars(0);
  const RealVector& refData = _model->getReferencePhysicalData();
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  CFreal p = _model->getPressureFromState(state[0]);
  CFreal T = state[4];
  CFreal pdim = p*refData[EulerTerm::P];
  CFreal Tdim = T*_model->getTempRef();
  
  // set the composition
  _ye.resize(nbElements);
  for(CFuint ie = 0; ie < nbElements; ++ie){
    _ye[ie] = state[5 + ie];
  }

  library->setElemFractions(_ye);
  library->setComposition(Tdim,pdim);
  library->setDensityEnthalpyEnergy(Tdim,pdim,_dhe);
  const CFreal rho = _dhe[0]/refData[EulerTerm::RHO];
  const CFreal V2 = u*u + v*v + w*w;
  
  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = rho*w;
  result[4] = rho*(_dhe[2]/refData[EulerTerm::H] + 0.5*V2);
  // States for the elemental fractinos
  for(CFuint ie = 0; ie < nbElements; ++ie){
    result[5 + ie] = rho*_ye[ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixPvtToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFuint nbElements = _model->getNbScalarVars(0);
  const CFuint firstElement = _model->getFirstScalarVar(0);
  const CFreal rho = data[EulerTerm::RHO];

  result[0] = rho;
  result[1] = rho*data[EulerTerm::VX];
  result[2] = rho*data[EulerTerm::VY];
  result[3] = rho*data[EulerTerm::VZ];
  result[4] = rho*data[EulerTerm::E];
  //rho*data.avH - data.avP; gives slightly different result
  // States for the elemental fractions
  for(CFuint ie = 0; ie < nbElements; ++ie){
    result[5 + ie] = rho*data[firstElement + ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
