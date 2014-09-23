#include "LTE.hh"
#include "IncompEuler2DLTEDemixCons.hh"
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

Environment::ObjectProvider<IncompEuler2DLTEDemixCons, ConvectiveVarSet, LTEModule, 1>
incompEuler2DLTEDemixConsProvider("IncompEuler2DLTEDemixCons");

//////////////////////////////////////////////////////////////////////////////

IncompEuler2DLTEDemixCons::IncompEuler2DLTEDemixCons(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<IncompEuler2DVarSet>(term),
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

IncompEuler2DLTEDemixCons::~IncompEuler2DLTEDemixCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixCons::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"IncompEuler2DLTEDemixCons::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixCons::computeEigenValuesVectors(RealMatrix& rightEv,
                                           RealMatrix& leftEv,
                                           RealVector& eValues,
                                           const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"IncompEuler2DLTEDemixCons::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint IncompEuler2DLTEDemixCons::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixCons::splitJacobian(RealMatrix& jacobPlus,
                                RealMatrix& jacobMin,
                                RealVector& eValues,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"IncompEuler2DLTEDemixCons::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixCons::computePhysicalData
(const State& state, RealVector& data)
{

  throw Common::NotImplementedException (FromHere(),"IncompEuler2DLTEDemixCons::computePhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixCons::computeStateFromPhysicalData
(const RealVector& data, State& state)
{
  throw Common::NotImplementedException (FromHere(),"IncompEuler2DLTEDemixCons::computeStateFromPhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

CFreal IncompEuler2DLTEDemixCons::getSpeed(const State& state) const
{
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];
  return sqrt(u*u + v*v);
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixCons::setDimensionalValues(const State& state,
                                          RealVector& result)
{
  throw Common::NotImplementedException (FromHere(),"IncompEuler2DLTEDemixCons::setDimensionalValues()");
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixCons::setAdimensionalValues(const Framework::State& state,
                             RealVector& result)
{
  throw Common::NotImplementedException
      (FromHere(), "IncompEuler2DLTEDemixCons::setAdimensionalValues() not implemented");
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixCons::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  throw Common::NotImplementedException (FromHere(),"IncompEuler2DLTEDemixCons::setDimensionalValuesPlusExtraValues()");
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> IncompEuler2DLTEDemixCons::getExtraVarNames() const
{
  cf_assert (_library.isNotNull());
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

void IncompEuler2DLTEDemixCons::setup()
{
  MultiScalarVarSet<IncompEuler2DVarSet>::setup();
  
  // set the equation set data for each of the equation subsets
  // first equation subset
  IncompEuler2DVarSet::getEqSetData().resize(1);
  IncompEuler2DVarSet::getEqSetData()[0].setup(0,0,4);
  
  // second equation subset
  MultiScalarVarSet<IncompEuler2DVarSet>::getEqSetData().resize(1);
  MultiScalarVarSet<IncompEuler2DVarSet>::getEqSetData()[0].setup
    (1,4,getModel()->getNbScalarVars(0));
  
  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert (_library.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void IncompEuler2DLTEDemixCons::computePerturbedPhysicalData
(const Framework::State& state,
 const RealVector& pdataBkp,
 RealVector& pdata,
 CFuint iVar)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DChemNEQCons::computePerturbedStatesData()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
