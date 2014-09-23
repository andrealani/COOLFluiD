#include "NEQ.hh"
#include "Euler1DNEQRhoivtTvToCons.hh"
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

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler1DNEQRhoivtTvToCons, VarSetTransformer, NEQModule, 1>
euler1DNEQRhoivtTvToConsProvider("Euler1DNEQRhoivtTvToCons");

//////////////////////////////////////////////////////////////////////////////

Euler1DNEQRhoivtTvToCons::Euler1DNEQRhoivtTvToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _ye(),
  _dhe(),
  _tvDim(),
  _evDim()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler1DNEQRhoivtTvToCons::~Euler1DNEQRhoivtTvToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivtTvToCons::transform(const State& state, State& result)
{

  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  const CFuint nbSpecies = _model->getNbScalarVars(0);
  // partial densities
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    result[ie] = state[ie];
  }

  // Set the mixture density (sum of the partial densities)
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    rho += state[ie];
  }

  // Set the species
  const CFreal ovRho = 1./rho;
  _ye.resize(nbSpecies);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ye[ie] = state[ie]*ovRho;
  }

  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!!
  library->setSpeciesFractions(_ye);

  const RealVector& refData =  _model->getReferencePhysicalData();
  CFreal rhoDim = rho*refData[EulerTerm::RHO];
  CFreal T = state[nbSpecies + 1];
  CFreal Tdim = T*refData[EulerTerm::T];

  const CFreal u = state[nbSpecies];
  const CFreal V2 = u*u;
  // set the vibrational temperature
  const CFuint nbTv = _model->getNbScalarVars(1);
  _dhe.resize(3 + nbTv);
  _tvDim.resize(nbTv);
  _evDim.resize(nbTv);

  const CFuint startTv = nbSpecies + 2;
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    _tvDim[ie] = state[startTv + ie]*refData[EulerTerm::T];
  }

  CFreal pdim = library->pressure(rhoDim, Tdim, &_tvDim[0]);
  const CFreal ovHref = 1./refData[EulerTerm::H];
  // vibrational temperatures
  library->setDensityEnthalpyEnergy(Tdim, _tvDim, pdim,_dhe);

  const CFuint nbTe = library->getNbTe();
  const CFuint nbTvH = nbTv - nbTe; 

  // data stores the moleculare vibrational energy multiplied 
  // by the molecules mass fractions
  if (nbTvH != 0) {
      for(CFuint i = 0; i < nbTvH; ++i) {
          result[startTv + i] = rho*_dhe[3 + i]*ovHref; 
      }
  }

  if (nbTe == 1) {
    result[startTv + nbTvH] = rho*_dhe[3 + nbTvH]*ovHref;
  }

  result[nbSpecies] = rho*u;
  result[nbSpecies + 1] = rho*(_dhe[2]*ovHref + 0.5*V2);
}

//////////////////////////////////////////////////////////////////////////////

void Euler1DNEQRhoivtTvToCons::transformFromRef(const RealVector& data, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  const CFuint firstSpecies = _model->getFirstScalarVar(0);
  const CFreal rho = data[EulerTerm::RHO];
//   const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
// unused //  const CFuint iEqSS = eqSS.getEqSS();
// unused //  const CFuint nbEqSS = eqSS.getTotalNbEqSS();
//   const RealVector& refData =  _model->getReferencePhysicalData();

  // partial densities
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    result[ie] = rho*data[firstSpecies + ie];
  }

  result[nbSpecies] = rho*data[EulerTerm::VX];
  result[nbSpecies + 1] = rho*data[EulerTerm::E];

  // set the vibrational temperature
  const CFuint nbTv = _model->getNbScalarVars(1);
  const CFuint startTv = nbSpecies + 2;
  const CFuint firstTv = _model->getFirstScalarVar(1);

  for(CFuint ie = 0; ie < nbTv; ++ie) {
    // vibrational temperatures (in data there is y_i*E_i)
    result[startTv + ie] = rho*data[firstTv + ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
