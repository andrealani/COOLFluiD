#include "Environment/ObjectProvider.hh"
#include "Maxwell/Maxwell.hh"
#include "Maxwell/Maxwell2DCons.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Maxwell2DCons, ConvectiveVarSet, MaxwellModule, 1>
maxwell2DConsProvider("Maxwell2DCons");

//////////////////////////////////////////////////////////////////////////////

Maxwell2DCons::Maxwell2DCons(Common::SafePtr<BaseTerm> term) :
  Maxwell2DVarSet(term),
  _rightEv(6,6),
  _leftEv(6,6)
{
  vector<std::string> names(6);
  names[0] = "Bx";
  names[1] = "By";
  names[2] = "Bz";
  names[3] = "Ex";
  names[4] = "Ey";
  names[5] = "Ez";
  
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Maxwell2DCons::~Maxwell2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DCons::setup()
{
  CFLog(NOTICE,"Maxwell2DCons::setup()\n");  
  
  Maxwell2DVarSet::setup();

  setConstJacob();
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DCons::setConstJacob()
{
  //CFLog(NOTICE, "setConstJacob()\n");  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DCons::computeProjectedJacobian(const RealVector& normal,
					   RealMatrix& jacob)
{
  //CFLog(NOTICE, "computeProjectedJacobian\n");
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DCons::computeJacobians()
{
  //CFLog(NOTICE, "computeJacobians()\n");
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DCons::computeEigenValuesVectors(RealMatrix& rightEv,
					    RealMatrix& leftEv,
					    RealVector& eValues,
					    const RealVector& normal)
{
  //CFLog(NOTICE, "computeEigenValuesVectors\n");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Maxwell2DCons::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DCons::splitJacobian(RealMatrix& jacobPlus,
           RealMatrix& jacobMin,
           RealVector& eValues,
           const RealVector& normal)
{
  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DCons::computePhysicalData(const State& state, RealVector& data)
{  
  //CFLog(NOTICE,"Maxwell2DCons::computePhysicalData()\n");
  data[ConvMaxwellTerm::BX] = state[0];
  data[ConvMaxwellTerm::BY] = state[1];
  data[ConvMaxwellTerm::BZ] = state[2];
  data[ConvMaxwellTerm::EX] = state[3];
  data[ConvMaxwellTerm::EY] = state[4];
  data[ConvMaxwellTerm::EZ] = state[5];
   
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DCons::computeStateFromPhysicalData(const RealVector& data,
					       State& state)
{
  state[0] = data[ConvMaxwellTerm::BX];
  state[1] = data[ConvMaxwellTerm::BY];
  state[2] = data[ConvMaxwellTerm::BZ];
  state[3] = data[ConvMaxwellTerm::EX];
  state[4] = data[ConvMaxwellTerm::EY];
  state[5] = data[ConvMaxwellTerm::EZ];
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DCons::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFreal refB = refData[ConvMaxwellTerm::BX];
  const CFreal refE = refData[ConvMaxwellTerm::EX];
  
  //CFLog(NOTICE, "setDimensionalValues\n");
  
  result[0] = state[0]*refB;
  result[1] = state[1]*refB;
  result[2] = state[2]*refB;
  result[3] = state[3]*refE;
  result[4] = state[4]*refE;
  result[5] = state[5]*refE;
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DCons::setAdimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFreal refB = refData[ConvMaxwellTerm::BX];
  const CFreal refE = refData[ConvMaxwellTerm::EX];

  //CFLog(NOTICE, "setAdimensionalValues\n");  
  
  result[0] = state[0]/refB;
  result[1] = state[1]/refB;
  result[2] = state[2]/refB;
  result[3] = state[3]/refE;
  result[4] = state[4]/refE;
  result[5] = state[5]/refE;
}
      
//////////////////////////////////////////////////////////////////////////////

void Maxwell2DCons::computePerturbedPhysicalData(const Framework::State& state,
						 const RealVector& pdataBkp,
						 RealVector& pdata,
						 CFuint iVar)
{
  throw Common::NotImplementedException
      (FromHere(), "Maxwell2DCons::computePerturbedPhysicalData() not implemented");
}

//////////////////////////////////////////////////////////////////////////////

bool Maxwell2DCons::isValid(const RealVector& data)
{
  return true;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
