#include "LTE.hh"
#include "Euler2DLTEDemixCons.hh"
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

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DLTEDemixCons, ConvectiveVarSet, LTEModule, 1>
euler2DLTEDemixConsProvider("Euler2DLTEDemixCons");

//////////////////////////////////////////////////////////////////////////////

Euler2DLTEDemixCons::Euler2DLTEDemixCons(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<Euler2DVarSet>(term),
  _library(CFNULL),
  _dhe(3),
  _x(),
  _ye()
{
  const CFuint nbElements = getModel()->getNbScalarVars(0);

  vector<std::string> names(4 + nbElements);
  names[0] = "rho";
  names[1] = "rhoU";
  names[2] = "rhoV";
  names[3] = "rhoE";

  // Names for the elemental fractions
  for (CFuint ie = 0; ie < nbElements; ++ie) {
    names[4 + ie] = "rhoYe" + StringOps::to_str(ie);
  }

  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DLTEDemixCons::~Euler2DLTEDemixCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixCons::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"Euler2DLTEDemixCons::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixCons::computeEigenValuesVectors(RealMatrix& rightEv,
                                           RealMatrix& leftEv,
                                           RealVector& eValues,
                                           const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DLTEDemixCons::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DLTEDemixCons::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixCons::splitJacobian(RealMatrix& jacobPlus,
                                RealMatrix& jacobMin,
                                RealVector& eValues,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DLTEDemixCons::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixCons::computePhysicalData
(const State& state, RealVector& data)
{

  throw Common::NotImplementedException (FromHere(),"Euler2DLTEDemixCons::computePhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixCons::computeStateFromPhysicalData
(const RealVector& data, State& state)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DLTEDemixCons::computeStateFromPhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler2DLTEDemixCons::getSpeed(const State& state) const
{
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];
  return sqrt(u*u + v*v);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixCons::setDimensionalValues(const State& state,
                                          RealVector& result)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DLTEDemixCons::setDimensionalValues()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixCons::setAdimensionalValues(const Framework::State& state,
                             RealVector& result)
{
  throw Common::NotImplementedException
      (FromHere(), "Euler2DLTEDemixCons::setAdimensionalValues() not implemented");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixCons::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DLTEDemixCons::setDimensionalValuesPlusExtraValues()");
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> Euler2DLTEDemixCons::getExtraVarNames() const
{
  cf_assert(_library.isNotNull());
  const CFuint nbSpecies = _library->getNbSpecies();

  vector<std::string> names(3 + nbSpecies);
  names[0] = "rho";
  names[1] = "H";
  names[2] = "M";

  for (CFuint is = 0; is < nbSpecies; ++is) {
    names[3 + is] = "xc" + StringOps::to_str(is);
  }

  return names;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixCons::setup()
{
  MultiScalarVarSet<Euler2DVarSet>::setup();
  
  // set the equation set data for each of the equation subsets
  // first equation subset
  Euler2DVarSet::getEqSetData().resize(1);
  Euler2DVarSet::getEqSetData()[0].setup(0,0,4);
  
  // second equation subset
  MultiScalarVarSet<Euler2DVarSet>::getEqSetData().resize(1);
  MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[0].setup
    (1,4,getModel()->getNbScalarVars(0));
  
  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixCons::computePerturbedPhysicalData
(const Framework::State& currState,
 const RealVector& bData,
 RealVector& data,
 CFuint iVar)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DLTEDemixCons::computePerturbedPhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DLTEDemixCons::setStateVelocityIDs (std::vector<CFuint>& velIDs) 
{
  velIDs.resize(2); velIDs[XX] = 1; velIDs[YY] = 2;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
