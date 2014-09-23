#include "NEQ.hh"
#include "Euler2DNEQRhoivtTvToRoe.hh"
#include "Environment/ObjectProvider.hh"
#include "NavierStokes/EulerPhysicalModel.hh"
#include "Framework/PhysicalChemicalLibrary.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQRhoivtTvToRoe, VarSetTransformer, NEQModule, 1>
euler2DNEQRhoivtTvToRoeProvider("Euler2DNEQRhoivtTvToRoe");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRhoivtTvToRoe::Euler2DNEQRhoivtTvToRoe
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

Euler2DNEQRhoivtTvToRoe::~Euler2DNEQRhoivtTvToRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhoivtTvToRoe::transform(const State& state, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  const CFuint nbSpecies = _model->getNbScalarVars(0);
// unused //  const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
// unused //  const CFuint iEqSS = eqSS.getEqSS();
// unused //  const CFuint nbEqSS = eqSS.getTotalNbEqSS();
  
  // Set the mixture density (sum of the partial densities)
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    rho += state[ie];
  }
  
  const CFreal sqRho = sqrt(rho);
  const CFreal ovSqRho = 1./sqRho;
  // partial densities
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    result[ie] = state[ie]*ovSqRho;
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
  CFreal T = state[nbSpecies + 2];
  CFreal Tdim = T*refData[EulerTerm::T];
 
// unused //  CFreal p = pdim/refData[EulerTerm::P];
  const CFreal u = state[nbSpecies];
  const CFreal v = state[nbSpecies + 1];
  const CFreal V2 = u*u + v*v;
  // set the vibrational temperature
  const CFuint nbTv = _model->getNbScalarVars(1);
  cf_assert(nbTv == 1);
  _dhe.resize(3 + nbTv);
  _tvDim.resize(nbTv);
  _evDim.resize(nbTv);
  
  const CFuint startTv = nbSpecies + 3;
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    _tvDim[ie] = state[startTv + ie]*refData[EulerTerm::T];
  }
  
  CFreal pdim = library->pressure(rhoDim, Tdim, &_tvDim[0]);
  
  const CFreal ovHref = 1./refData[EulerTerm::H];
  // vibrational temperatures
  library->setDensityEnthalpyEnergy(Tdim, _tvDim, pdim,_dhe, true);
  
  if (nbTv > 1) {
    cf_assert(false);
    // static vector<CFuint> moleculesIDs;
    //     library->setMoleculesIDs(moleculesIDs);
    //     cf_assert(moleculesIDs.size() == nbTv);
    
    //     for(CFuint ie = 0; ie < nbTv; ++ie) {
    //       result[startTv + ie] = _ye[moleculesIDs[ie]]*rho*_dhe[3 + ie]*ovHref;
    //     }
  }
  else if (nbTv == 1) {
    result[startTv] = sqRho*_dhe[3]*ovHref;
  }
  
  result[nbSpecies] = sqRho*u;
  result[nbSpecies + 1] = sqRho*v;
  result[nbSpecies + 2] = sqRho*(_dhe[1]*ovHref + 0.5*V2);
  
  // set extra data needed for the linearization in Roe variables
  // cf_assert(_extraValues != CFNULL);
//   SafePtr<PhysicalChemicalLibrary::ExtraData> eData = library->getExtraData();
//   RealVector& extra = (*_extraValues)[_iState];
//   extra.resize(nbSpecies + 2);
//   for (CFuint ie = 0; ie < nbSpecies; ++ie) {
//     extra[ie] = state[ie];
//   }
//   extra[nbSpecies] = eData->energyTr.sum();
//   extra[nbSpecies + 1] = p;
}

//////////////////////////////////////////////////////////////////////////////
      
void Euler2DNEQRhoivtTvToRoe::transformFromRef(const RealVector& data, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  const CFuint firstSpecies = _model->getFirstScalarVar(0);
  const CFreal sqRho = sqrt(data[EulerTerm::RHO]);
// unused //  const EquationSubSysDescriptor& eqSS = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
// unused //  const CFuint iEqSS = eqSS.getEqSS();
// unused //  const CFuint nbEqSS = eqSS.getTotalNbEqSS();
// unused //  const RealVector& refData =  _model->getReferencePhysicalData();
  
  // partial densities
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    result[ie] = sqRho*data[firstSpecies + ie];
  }
  
  result[nbSpecies] = sqRho*data[EulerTerm::VX];
  result[nbSpecies + 1] = sqRho*data[EulerTerm::VY];
  result[nbSpecies + 2] = sqRho*data[EulerTerm::H];
  
  // set the vibrational temperature
  const CFuint nbTv = _model->getNbScalarVars(1);
  const CFuint startTv = nbSpecies + 3;  
  const CFuint firstTv = _model->getFirstScalarVar(1);
  
// unused //  CFreal pdim = data[EulerTerm::P]*refData[EulerTerm::P];
// unused //  CFreal Tdim = data[EulerTerm::T]*refData[EulerTerm::T];
  for(CFuint ie = 0; ie < nbTv; ++ie) {
    // vibrational temperatures (in data there is y_i*E_i)
    result[startTv + ie] = sqRho*data[firstTv + ie];
  }
  
  if (_extraValues != CFNULL) {
    throw Common::NotImplementedException (FromHere(),"Euler2DNEQRhoivtTvToRoe::transformFromRef()");
  }
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
