#include "NEQ.hh"
#include "Euler2DNEQPvtyToCons.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/EulerPhysicalModel.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQPvtyToCons, VarSetTransformer, NEQModule, 1>
euler2DNEQPvtyToConsProvider("Euler2DNEQPvtyToCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQPvtyToCons::Euler2DNEQPvtyToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerNEQTerm>()),
  _dhe(3),
  _ye()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQPvtyToCons::~Euler2DNEQPvtyToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPvtyToCons::transform(const State& state, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  const CFuint nbSpecies = _model->getNbScalarVars(0);
  const EquationSubSysDescriptor& eqSS =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint iEqSS = eqSS.getEqSS();

  const RealVector& refData = _model->getReferencePhysicalData();
  CFreal p = state[0];
  CFreal T = state[3];
  CFreal pdim = p*refData[EulerTerm::P];
  CFreal Tdim = T*_model->getTempRef();
  CFreal rho = 0.0;

  _ye.resize(nbSpecies);
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    _ye[ie] = state[4 + ie];
  }

  library->setSpeciesFractions(_ye);

 if (iEqSS == 0) {
   const CFreal u = state[1];
   const CFreal v = state[2];
   library->setDensityEnthalpyEnergy(Tdim, pdim,_dhe);
   rho = _dhe[0]/refData[EulerTerm::RHO];
   const CFreal V2 = u*u + v*v;
   
   result[0] = rho;
   result[1] = rho*u;
   result[2] = rho*v;
   result[3] = rho*(_dhe[2]/refData[EulerTerm::H] + 0.5*V2);
 }
 
 if (iEqSS == 1) {
   rho = library->density(Tdim,pdim)/refData[EulerTerm::RHO];
 }

 const CFuint nbEqs = eqSS.getNbEqsSS();
 if ((nbEqs == 4 + nbSpecies) || (iEqSS == 1)) {
   // States for the species mass fractions
   for(CFuint ie = 0; ie < nbSpecies; ++ie){
     result[4 + ie] = rho*_ye[ie];
   }
 }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQPvtyToCons::transformFromRef(const RealVector& data, State& result) 
{
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  const CFuint firstSpecies = _model->getFirstScalarVar(0);
  const CFreal rho = data[EulerTerm::RHO];
  const EquationSubSysDescriptor& eqSS =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint iEqSS = eqSS.getEqSS();

  if (iEqSS == 0) {
    result[0] = rho;
    result[1] = rho*data[EulerTerm::VX];
    result[2] = rho*data[EulerTerm::VY];
    result[3] = rho*data[EulerTerm::E];
  }

  const CFuint nbEqs = eqSS.getNbEqsSS();
  if ((nbEqs == 4 + nbSpecies) || (iEqSS == 1)) {
    // States for the species mass fractions
    for(CFuint ie = 0; ie < nbSpecies; ++ie){
      result[4 + ie] = rho*data[firstSpecies + ie];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
