#include "Environment/ObjectProvider.hh"
#include "Maxwell/Maxwell.hh"
#include "Maxwell/Maxwell3DProjectionCons.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Common/NotImplementedException.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Maxwell3DProjectionCons, ConvectiveVarSet, MaxwellModule, 1>
maxwell3DProjectionConsProvider("Maxwell3DProjectionCons");

//////////////////////////////////////////////////////////////////////////////

Maxwell3DProjectionCons::Maxwell3DProjectionCons(Common::SafePtr<BaseTerm> term) :
  Maxwell3DProjectionVarSet(term),
  _rightEv(8,8),
  _leftEv(8,8)
{
  vector<std::string> names(8);
  names[0] = "Bx";
  names[1] = "By";
  names[2] = "Bz";
  names[3] = "Ex";
  names[4] = "Ey";
  names[5] = "Ez";
  names[6] = "Psi";
  names[7] = "Phi";
  
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Maxwell3DProjectionCons::~Maxwell3DProjectionCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DProjectionCons::setup()
{
  CFLog(NOTICE,"Maxwell3DProjectionCons::setup()\n");  
  
  //Maxwell3DProjectionVarSet::setup();

  setConstJacob();
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> Maxwell3DProjectionCons::getExtraVarNames() const
{
  vector<std::string> names(8);
  names[0] = "Bx";
  names[1] = "By";
  names[2] = "Bz";
  names[3] = "Ex";
  names[4] = "Ey";
  names[5] = "Ez";
  names[6] = "Psi";
  names[7] = "Phi";
  return names;
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DProjectionCons::setConstJacob()
{
  //CFLog(NOTICE, "setConstJacob()\n");  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DProjectionCons::computeProjectedJacobian(const RealVector& normal,
					   RealMatrix& jacob)
{
  //CFLog(NOTICE, "computeProjectedJacobian\n");
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DProjectionCons::computeJacobians()
{
  //CFLog(NOTICE, "computeJacobians()\n");
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DProjectionCons::computeEigenValuesVectors(RealMatrix& rightEv,
					    RealMatrix& leftEv,
					    RealVector& eValues,
					    const RealVector& normal)
{
  //CFLog(NOTICE, "computeEigenValuesVectors\n");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Maxwell3DProjectionCons::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DProjectionCons::splitJacobian(RealMatrix& jacobPlus,
           RealMatrix& jacobMin,
           RealVector& eValues,
           const RealVector& normal)
{
  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DProjectionCons::computePhysicalData(const State& state, RealVector& data)
{  
//   CFLog(NOTICE, "computePhysicalData\n");  
  data[ConvMaxwellTerm::BX] = state[0];
  data[ConvMaxwellTerm::BY] = state[1];
  data[ConvMaxwellTerm::BZ] = state[2];
  data[ConvMaxwellTerm::EX] = state[3];
  data[ConvMaxwellTerm::EY] = state[4];
  data[ConvMaxwellTerm::EZ] = state[5];
  data[MaxwellProjectionTerm::PSI] = state[6];  
  data[MaxwellProjectionTerm::PHI] = state[7];
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DProjectionCons::computeStateFromPhysicalData(const RealVector& data,
					       State& state)
{
  
//   CFLog(NOTICE, "computeStateFromPhysicalData\n");
  
  state[0] = data[ConvMaxwellTerm::BX];
  state[1] = data[ConvMaxwellTerm::BY];
  state[2] = data[ConvMaxwellTerm::BZ];
  state[3] = data[ConvMaxwellTerm::EX];
  state[4] = data[ConvMaxwellTerm::EY];
  state[5] = data[ConvMaxwellTerm::EZ];
  state[6] = data[MaxwellProjectionTerm::PSI];
  state[7] = data[MaxwellProjectionTerm::PHI];
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DProjectionCons::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFreal refB = refData[ConvMaxwellTerm::BX];
  const CFreal refE = refData[ConvMaxwellTerm::EX];
  const CFreal refPsi = refData[MaxwellProjectionTerm::PSI];
  const CFreal refPhi = refData[MaxwellProjectionTerm::PHI];
  
//   CFLog(NOTICE, "setDimensionalValues\n");
  
  result[0] = state[0]*refB;
  result[1] = state[1]*refB;
  result[2] = state[2]*refB;
  result[3] = state[3]*refE;
  result[4] = state[4]*refE;
  result[5] = state[5]*refE;
  result[6] = state[6]*refPsi;
  result[7] = state[7]*refPhi;
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DProjectionCons::setAdimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFreal refB = refData[ConvMaxwellTerm::BX];
  const CFreal refE = refData[ConvMaxwellTerm::EX];
  const CFreal refPsi = refData[MaxwellProjectionTerm::PSI];
  const CFreal refPhi = refData[MaxwellProjectionTerm::PHI];

//   CFLog(NOTICE, "setAdimensionalValues\n");  
  
 result[0] = state[0]/refB;
 result[1] = state[1]/refB;
 result[2] = state[2]/refB;
 result[3] = state[3]/refE;
 result[4] = state[4]/refE;
 result[5] = state[5]/refE;
 result[6] = state[6]/refPsi;
 result[7] = state[7]/refPhi;
    
}
      
//////////////////////////////////////////////////////////////////////////////

void Maxwell3DProjectionCons::computePerturbedPhysicalData(const Framework::State& state,
						 const RealVector& pdataBkp,
						 RealVector& pdata,
						 CFuint iVar)
{
  throw Common::NotImplementedException
      (FromHere(), "Maxwell3DProjectionCons::computePerturbedPhysicalData() not implemented");
}

//////////////////////////////////////////////////////////////////////////////

bool Maxwell3DProjectionCons::isValid(const RealVector& data)
{
  return true;
}

/////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD
