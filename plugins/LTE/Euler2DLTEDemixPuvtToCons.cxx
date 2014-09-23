#include "LTE.hh"
#include "Euler2DLTEDemixPuvtToCons.hh"
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

Environment::ObjectProvider<Euler2DLTEDemixPuvtToCons, VarSetTransformer, LTEModule, 1> 
euler2DLTEDemixPuvtToConsProvider("Euler2DLTEDemixPuvtToCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DLTEDemixPuvtToCons::Euler2DLTEDemixPuvtToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerLTEDemixTerm>()),
  _dhe(3),
  _ye()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DLTEDemixPuvtToCons::~Euler2DLTEDemixPuvtToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixPuvtToCons::transform(const State& state, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  const CFuint nbElements = _model->getNbScalarVars(0);
  const CFuint nbEqs =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor().getNbEqsSS();

  const RealVector& refData = _model->getReferencePhysicalData();
  CFreal p = _model->getPressureFromState(state[0]);
  CFreal T = state[3];
  CFreal pdim = p*refData[EulerTerm::P];
  CFreal Tdim = T*_model->getTempRef();
  CFreal rho = 0.0;
  
  // set the composition
  _ye.resize(nbElements);
  for(CFuint ie = 0; ie < nbElements; ++ie){
    _ye[ie] = state[4 + ie];
  }
  
  library->setElemFractions(_ye);
  library->setComposition(Tdim,pdim);
  
  if (nbEqs >= 4) {
    const CFreal u = state[1];
    const CFreal v = state[2];
    library->setDensityEnthalpyEnergy(Tdim,pdim,_dhe);
    rho = _dhe[0]/refData[EulerTerm::RHO];
    const CFreal V2 = u*u + v*v;
    
    result[0] = rho;
    result[1] = rho*u;
    result[2] = rho*v;
    result[3] = rho*(_dhe[2]/refData[EulerTerm::H] + 0.5*V2);
  }
  
  if (nbEqs < 4) {
    rho = library->density(Tdim,pdim)/refData[EulerTerm::RHO];
  }
  
  if (nbEqs != 4) {
    // States for the elemental fractions
    for(CFuint ie = 0; ie < nbElements; ++ie){
      result[4 + ie] = rho*_ye[ie];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixPuvtToCons::transformFromRef(const RealVector& data, State& result)
{
  const CFuint nbElements = _model->getNbScalarVars(0);
  const CFuint firstElement = _model->getFirstScalarVar(0);
  const CFreal rho = data[EulerTerm::RHO];
  const CFuint nbEqs =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor().getNbEqsSS();
  
  if (nbEqs >= 4) {
    result[0] = rho;
    result[1] = rho*data[EulerTerm::VX];
    result[2] = rho*data[EulerTerm::VY];
    result[3] = rho*data[EulerTerm::E];
  }
  
  if (nbEqs != 4) {
    // States for the elemental fractinos
    for(CFuint ie = 0; ie < nbElements; ++ie){
      result[4 + ie] = rho*data[firstElement + ie];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
