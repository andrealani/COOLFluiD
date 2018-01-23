#include "Environment/ObjectProvider.hh"
#include "Poisson/Poisson.hh"
#include "Poisson/PoissonConvCons.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Poisson {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<PoissonConvCons, ConvectiveVarSet, PoissonModule, 1>
poissonConv1DConsProvider("PoissonConv1DCons");

Environment::ObjectProvider<PoissonConvCons, ConvectiveVarSet, PoissonModule, 1>
poissonConv2DConsProvider("PoissonConv2DCons");

Environment::ObjectProvider<PoissonConvCons, ConvectiveVarSet, PoissonModule, 1>
poissonConv3DConsProvider("PoissonConv3DCons");

//////////////////////////////////////////////////////////////////////////////

PoissonConvCons::PoissonConvCons(Common::SafePtr<BaseTerm> term) :
  PoissonConvVarSet(term),
  _rightEv(1,1),
  _leftEv(1,1)
{
  vector<std::string> names(1);
  names[0] = "Phi";
  
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

PoissonConvCons::~PoissonConvCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void PoissonConvCons::setup()
{
  CFLog(NOTICE,"PoissonConvCons::setup()\n");
  
  PoissonConvVarSet::setup();
  
  setConstJacob();
}

//////////////////////////////////////////////////////////////////////////////

void PoissonConvCons::setConstJacob()
{
}
      
//////////////////////////////////////////////////////////////////////////////

void PoissonConvCons::computeProjectedJacobian(const RealVector& normal,
					   RealMatrix& jacob)
{
}
      
//////////////////////////////////////////////////////////////////////////////
      
void PoissonConvCons::computeJacobians()
{
}
      
//////////////////////////////////////////////////////////////////////////////

void PoissonConvCons::computeEigenValuesVectors(RealMatrix& rightEv,
					    RealMatrix& leftEv,
					    RealVector& eValues,
					    const RealVector& normal)
{
}

//////////////////////////////////////////////////////////////////////////////

CFuint PoissonConvCons::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void PoissonConvCons::splitJacobian(RealMatrix& jacobPlus,
				    RealMatrix& jacobMin,
				    RealVector& eValues,
				    const RealVector& normal)
{
}
      
//////////////////////////////////////////////////////////////////////////////

void PoissonConvCons::computePhysicalData(const State& state, RealVector& data)
{  
  data[PoissonConvTerm::PHI] = state[0];
}

//////////////////////////////////////////////////////////////////////////////

void PoissonConvCons::computeStateFromPhysicalData(const RealVector& data,
					       State& state)
{
  state[0] = data[PoissonConvTerm::PHI];
}

//////////////////////////////////////////////////////////////////////////////

void PoissonConvCons::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFreal refPhi = refData[PoissonConvTerm::PHI];
  result[0] = state[0]*refPhi;
}

//////////////////////////////////////////////////////////////////////////////

void PoissonConvCons::setAdimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFreal refPhi = refData[PoissonConvTerm::PHI];
  result[0] = state[0]/refPhi;
}
      
//////////////////////////////////////////////////////////////////////////////

void PoissonConvCons::computePerturbedPhysicalData(const Framework::State& state,
						 const RealVector& pdataBkp,
						 RealVector& pdata,
						 CFuint iVar)
{
  throw Common::NotImplementedException
      (FromHere(), "PoissonConvCons::computePerturbedPhysicalData() not implemented");
}

//////////////////////////////////////////////////////////////////////////////

bool PoissonConvCons::isValid(const RealVector& data)
{
  return true;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Poisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
