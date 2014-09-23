#include "NEQ.hh"
#include "Euler2DNEQRoe.hh"
#include "Common/NotImplementedException.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQRoe, ConvectiveVarSet, NEQModule, 1>
euler2DNEQRoeProvider("Euler2DNEQRoe");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRoe::Euler2DNEQRoe(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<Euler2DVarSet>(term),
  _library(CFNULL),
  _Rgas(0.0),
  _dhe(3),
  _ys(),
  _rightEv(),
  _leftEv(),
  _alpha(),
  _RiGas(),
  _tempVib()
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint nbTv = getModel()->getNbScalarVars(1);

  const CFuint nbEqs = nbSpecies + 3 + nbTv;
  vector<std::string> names(nbEqs);
  for (CFuint ie = 0; ie < nbEqs; ++ie) {
    names[ie] = "z" + StringOps::to_str(ie);
  }

  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRoe::~Euler2DNEQRoe()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRoe::computeEigenValuesVectors(RealMatrix& rightEv,
                                           RealMatrix& leftEv,
                                           RealVector& eValues,
                                           const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQRoe::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DNEQRoe::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRoe::splitJacobian(RealMatrix& jacobPlus,
                                RealMatrix& jacobMin,
                                RealVector& eValues,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQRoe::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRoe::computePhysicalData(const State& state,
					RealVector& data)
{

  cout << "In Euler2DNEQRoe::computePhysicalData(const State& , RealVector& )" << endl;
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);

  // Set the mixture density (sum of the partial densities)
  CFreal sqRho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    sqRho += state[ie];
  }
  const CFreal rho = sqRho*sqRho;
  const CFreal ovSqrho = 1./sqRho;

  // Set the species
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ys[ie] = state[ie]*ovSqrho;
  }

  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!!
  _library->setSpeciesFractions(_ys);

  // set the species mass fractions
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    data[firstSpecies + ie] = _ys[ie];
  }

  const CFreal cvTr = _library->getCvTr();
  const CFuint uID = nbSpecies;
  const CFuint vID = nbSpecies+1;
  const CFuint hID = nbSpecies+2;
  const CFuint evID = nbSpecies+3;
  const CFreal u = state[uID]*ovSqrho;
  const CFreal v = state[vID]*ovSqrho;
  const CFreal V2 = u*u + v*v;
  data[EulerTerm::V] = sqrt(V2);
  data[EulerTerm::VX] = u;
  data[EulerTerm::VY] = v;
  data[EulerTerm::RHO] = rho;

  const CFreal mmass = _library->getMMass();
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFreal T = (state[hID] - state[evID] - 0.5*sqRho*V2)/(rho*sqRho*_Rgas/mmass + cvTr);
  CFreal rhodim = rho*refData[EulerTerm::RHO];
  CFreal Tdim = T*refData[EulerTerm::T];
  _tempVib[0] = state[evID];
  CFreal pdim = _library->pressure(rhodim, Tdim, &_tempVib[0]);
  CFreal p = pdim/refData[EulerTerm::P];
  
  // this is inconsistent
  data[EulerTerm::P] = p;
  data[EulerTerm::T] = T;
  data[EulerTerm::H] = state[nbSpecies +2]*ovSqrho;
  data[EulerTerm::E] = data[EulerTerm::H] - p/rho;


  ///@TODO AL:  here you should pass the vibrational temperatures
  _library->frozenGammaAndSoundSpeed(Tdim, pdim, rhodim,
				     data[EulerTerm::GAMMA],
				     data[EulerTerm::A],
				     CFNULL); // this sound speed is inconsistent

  const CFuint firstTv = getModel()->getFirstScalarVar(1);
  data[firstTv] = state[evID]/sqRho;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRoe::computeStateFromPhysicalData(const RealVector& data,
						 State& state)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQRoe::computeStateFromPhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler2DNEQRoe::getSpeed(const State& state) const
{
  throw Common::NotImplementedException
    (FromHere(), "Euler2DNEQRoe::getSpeed()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRoe::setDimensionalValues(const State& state,
                                          RealVector& result)
{
  throw Common::NotImplementedException
      (FromHere(), "Euler2DNEQRoe::setDimensionalValues()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRoe::setAdimensionalValues(const Framework::State& state,
                             RealVector& result)
{
  throw Common::NotImplementedException
      (FromHere(), "Euler2DNEQRoe::setAdimensionalValues() not implemented");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRoe::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  throw Common::NotImplementedException
    (FromHere(), "Euler2DNEQRoe::setDimensionalValuesPlusExtraValues()");
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> Euler2DNEQRoe::getExtraVarNames() const
{
  cf_assert(_library.isNotNull());

  vector<std::string> names;
throw Common::NotImplementedException (FromHere(),"Euler2DNEQRoe::getExtraVarNames()");
  return names;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRoe::setup()
{
  MultiScalarVarSet<Euler2DVarSet>::setup();

  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  // set the equation set data for each of the equation subsets
  // first equation subset
  MultiScalarVarSet<Euler2DVarSet>::getEqSetData().resize(2);
  MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[0].setup(0,0,nbSpecies);

  // second equation subset
  Euler2DVarSet::getEqSetData().resize(1);
  Euler2DVarSet::getEqSetData()[0].setup(1,nbSpecies,3);

  const CFuint nbTv = getModel()->getNbScalarVars(1);
  // third equation subset
  MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[1].setup(2,nbSpecies + 3,nbTv);

  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());

  _Rgas = _library->getRgas();

  _dhe.resize(3 + nbTv);
  _ys.resize(nbSpecies);
  const CFuint totNbEqs = 3 + nbTv + nbSpecies;

  _rightEv.resize(totNbEqs, totNbEqs);
  _rightEv = 0.0;

  _leftEv.resize(totNbEqs, totNbEqs);
  _leftEv = 0.0;

  _alpha.resize(nbSpecies);
  _RiGas.resize(nbSpecies);
  _library->setRiGas(_RiGas);
  
  _tempVib.resize(nbTv);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRoe::computePerturbedPhysicalData
(const State& state, const RealVector& pdataBkp, 
 RealVector& pdata, CFuint iVar) 
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQRoe::computePerturbedPhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
