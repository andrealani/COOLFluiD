#include "LTE.hh"
#include "Euler3DLTEDemixCons.hh"
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

Environment::ObjectProvider<Euler3DLTEDemixCons, ConvectiveVarSet, LTEModule, 1>
euler3DLTEDemixConsProvider("Euler3DLTEDemixCons");

//////////////////////////////////////////////////////////////////////////////

Euler3DLTEDemixCons::Euler3DLTEDemixCons(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<Euler3DVarSet>(term),
  _library(CFNULL),
  _dhe(3),
  _x(),
  _ye()
{
  const CFuint nbElements = getModel()->getNbScalarVars(0);

  vector<std::string> names(5 + nbElements);
  names[0] = "rho";
  names[1] = "rhoU";
  names[2] = "rhoV";
  names[3] = "rhoW";
  names[4] = "rhoE";

  // Names for the elemental fractions
  for (CFuint ie = 0; ie < nbElements; ++ie) {
    names[5 + ie] = "rhoYe" + StringOps::to_str(ie);
  }

  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler3DLTEDemixCons::~Euler3DLTEDemixCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixCons::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"Euler3DLTEDemixCons::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixCons::computeEigenValuesVectors(RealMatrix& rightEv,
					   RealMatrix& leftEv,
					   RealVector& eValues,
					   const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DLTEDemixCons::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler3DLTEDemixCons::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixCons::splitJacobian(RealMatrix& jacobPlus,
				RealMatrix& jacobMin,
				RealVector& eValues,
				const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DLTEDemixCons::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixCons::computePhysicalData
(const State& state, RealVector& data)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DLTEDemixCons::computePhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixCons::computeStateFromPhysicalData
(const RealVector& data, State& state)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DLTEDemixCons::computeStateFromPhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler3DLTEDemixCons::getSpeed(const State& state) const
{
  return sqrt((state[1]*state[1] + state[2]*state[2] + state[3]*state[3])/(state[0]*state[0]));
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixCons::setDimensionalValues(const State& state,
					  RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  const CFreal rhoRef = refData[EulerTerm::RHO];
  result[0] = state[0]*rhoRef;
  result[1] = state[1]*refData[EulerTerm::V]*rhoRef;
  result[2] = state[2]*refData[EulerTerm::V]*rhoRef;
  result[3] = state[3]*refData[EulerTerm::V]*rhoRef;
  result[4] = state[4]*refData[EulerTerm::H]*rhoRef;

  // Set the elemental fractions
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbElements; ++ie){
    result[5 + ie] = state[5 + ie]*rhoRef;
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixCons::setAdimensionalValues(const State& state,
					  RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();

  const CFreal rhoRef = refData[EulerTerm::RHO];
  result[0] = state[0]/rhoRef;
  result[1] = state[1]/(refData[EulerTerm::V]*rhoRef);
  result[2] = state[2]/(refData[EulerTerm::V]*rhoRef);
  result[3] = state[3]/(refData[EulerTerm::V]*rhoRef);
  result[4] = state[4]/(refData[EulerTerm::H]*rhoRef);

  // Set the elemental fractions
  const CFuint nbElements = getModel()->getNbScalarVars(0);
  for (CFuint ie = 0; ie < nbElements; ++ie){
    result[5 + ie] = state[5 + ie]/rhoRef;
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixCons::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  throw Common::NotImplementedException
    (FromHere(), "Euler3DLTEDemixCons::setDimensionalValuesPlusExtraValues()");
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> Euler3DLTEDemixCons::getExtraVarNames() const
{
  throw Common::NotImplementedException (FromHere(),"Euler3DLTEDemixCons::getExtraVarNames()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixCons::setup()
{
  MultiScalarVarSet<Euler3DVarSet>::setup();
  
  // set the equation set data for each of the equation subsets
  // first equation subset
  Euler3DVarSet::getEqSetData().resize(1);
  Euler3DVarSet::getEqSetData()[0].setup(0,0,5);
  
  // second equation subset
  MultiScalarVarSet<Euler3DVarSet>::getEqSetData().resize(1);
  MultiScalarVarSet<Euler3DVarSet>::getEqSetData()[0].setup
    (1,5,getModel()->getNbScalarVars(0));
  
  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixCons::computePerturbedPhysicalData
(const Framework::State& currState,
 const RealVector& bData,
 RealVector& data,
 CFuint iVar)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DLTEDemixCons::computePerturbedPhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DLTEDemixCons::setStateVelocityIDs (std::vector<CFuint>& velIDs) 
{
  velIDs.resize(3); velIDs[XX] = 1; velIDs[YY] = 2; velIDs[ZZ] = 3;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace LTE

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
