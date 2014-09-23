#include "Environment/ObjectProvider.hh"
#include "Maxwell/Maxwell.hh"
#include "Maxwell/Maxwell3DAdimCons.hh"
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

Environment::ObjectProvider<Maxwell3DAdimCons, ConvectiveVarSet, MaxwellModule, 1>
maxwell3DAdimConsProvider("Maxwell3DAdimCons");

//////////////////////////////////////////////////////////////////////////////

Maxwell3DAdimCons::Maxwell3DAdimCons(Common::SafePtr<BaseTerm> term) :
  Maxwell3DVarSetAdim(term),
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

Maxwell3DAdimCons::~Maxwell3DAdimCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DAdimCons::setup()
{
  CFLog(NOTICE,"Maxwell3DAdimCons::setup()\n");  
  
  setConstJacob();
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> Maxwell3DAdimCons::getExtraVarNames() const
{
  vector<std::string> names(6);
  names[0] = "Bx";
  names[1] = "By";
  names[2] = "Bz";
  names[3] = "Ex";
  names[4] = "Ey";
  names[5] = "Ez";

  return names;
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DAdimCons::setConstJacob()
{
  //CFLog(NOTICE, "setConstJacob()\n");  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DAdimCons::computeProjectedJacobian(const RealVector& normal,
					   RealMatrix& jacob)
{
  //CFLog(NOTICE, "computeProjectedJacobian\n");
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DAdimCons::computeJacobians()
{
  //CFLog(NOTICE, "computeJacobians()\n");
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DAdimCons::computeEigenValuesVectors(RealMatrix& rightEv,
					    RealMatrix& leftEv,
					    RealVector& eValues,
					    const RealVector& normal)
{
  //CFLog(NOTICE, "computeEigenValuesVectors\n");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Maxwell3DAdimCons::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DAdimCons::splitJacobian(RealMatrix& jacobPlus,
           RealMatrix& jacobMin,
           RealVector& eValues,
           const RealVector& normal)
{
  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DAdimCons::computePhysicalData(const State& state, RealVector& data)
{  
  //CFLog(NOTICE, "computePhysicalData\n");  
  data[ConvMaxwellTerm::BX] = state[0];
  data[ConvMaxwellTerm::BY] = state[1];
  data[ConvMaxwellTerm::BZ] = state[2];
  data[ConvMaxwellTerm::EX] = state[3];
  data[ConvMaxwellTerm::EY] = state[4];
  data[ConvMaxwellTerm::EZ] = state[5];

}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DAdimCons::computeStateFromPhysicalData(const RealVector& data,
					       State& state)
{
  
  //CFLog(NOTICE, "computeStateFromPhysicalData\n");
  
  state[0] = data[ConvMaxwellTerm::BX];
  state[1] = data[ConvMaxwellTerm::BY];
  state[2] = data[ConvMaxwellTerm::BZ];
  state[3] = data[ConvMaxwellTerm::EX];
  state[4] = data[ConvMaxwellTerm::EY];
  state[5] = data[ConvMaxwellTerm::EZ];

}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DAdimCons::setDimensionalValues(const State& state,
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

void Maxwell3DAdimCons::setAdimensionalValues(const State& state,
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

void Maxwell3DAdimCons::computePerturbedPhysicalData(const Framework::State& state,
						 const RealVector& pdataBkp,
						 RealVector& pdata,
						 CFuint iVar)
{
  throw Common::NotImplementedException
      (FromHere(), "Maxwell3DAdimCons::computePerturbedPhysicalData() not implemented");
}

//////////////////////////////////////////////////////////////////////////////

bool Maxwell3DAdimCons::isValid(const RealVector& data)
{
  return true;
}

/////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD
