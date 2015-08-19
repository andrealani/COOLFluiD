#include "Environment/ObjectProvider.hh"
#include "Poisson/Poisson.hh"
#include "Poisson/PoissonConv3DCons.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Poisson {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<PoissonConv3DCons, ConvectiveVarSet, PoissonModule, 1>
PoissonConv3DConsProvider("PoissonConv3DCons");

//////////////////////////////////////////////////////////////////////////////

PoissonConv3DCons::PoissonConv3DCons(Common::SafePtr<BaseTerm> term) :
  PoissonConv3DVarSet(term),
  _rightEv(1,1),
  _leftEv(1,1)
{
  vector<std::string> names(1);
  names[0] = "Phi";
  
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

PoissonConv3DCons::~PoissonConv3DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void PoissonConv3DCons::setup()
{
  CFLog(NOTICE,"PoissonConv3DCons::setup()\n");
  
  PoissonConv3DVarSet::setup();

  setConstJacob();
}

//////////////////////////////////////////////////////////////////////////////

void PoissonConv3DCons::setConstJacob()
{
  //CFLog(NOTICE, "setConstJacob()\n");  
}

//////////////////////////////////////////////////////////////////////////////

void PoissonConv3DCons::computeProjectedJacobian(const RealVector& normal,
					   RealMatrix& jacob)
{
  //CFLog(NOTICE, "computeProjectedJacobian\n");
}

//////////////////////////////////////////////////////////////////////////////

void PoissonConv3DCons::computeJacobians()
{
  //CFLog(NOTICE, "computeJacobians()\n");
}

//////////////////////////////////////////////////////////////////////////////

void PoissonConv3DCons::computeEigenValuesVectors(RealMatrix& rightEv,
					    RealMatrix& leftEv,
					    RealVector& eValues,
					    const RealVector& normal)
{
  //CFLog(NOTICE, "computeEigenValuesVectors\n");
}

//////////////////////////////////////////////////////////////////////////////

CFuint PoissonConv3DCons::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void PoissonConv3DCons::splitJacobian(RealMatrix& jacobPlus,
           RealMatrix& jacobMin,
           RealVector& eValues,
           const RealVector& normal)
{
  
}

//////////////////////////////////////////////////////////////////////////////

void PoissonConv3DCons::computePhysicalData(const State& state, RealVector& data)
{  
  //CFLog(NOTICE,"PoissonConv3DCons::computePhysicalData()\n");
  data[PoissonConvTerm::PHI] = state[0];
   
}

//////////////////////////////////////////////////////////////////////////////

void PoissonConv3DCons::computeStateFromPhysicalData(const RealVector& data,
					       State& state)
{
  state[0] = data[PoissonConvTerm::PHI];

}

//////////////////////////////////////////////////////////////////////////////

void PoissonConv3DCons::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFreal refPhi = refData[PoissonConvTerm::PHI];
  
  //CFLog(NOTICE, "setDimensionalValues\n");
  
  result[0] = state[0]*refPhi;

}

//////////////////////////////////////////////////////////////////////////////

void PoissonConv3DCons::setAdimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData = getModel()->getReferencePhysicalData();
  const CFreal refPhi = refData[PoissonConvTerm::PHI];

  //CFLog(NOTICE, "setAdimensionalValues\n");  
  
  result[0] = state[0]/refPhi;

}
      
//////////////////////////////////////////////////////////////////////////////

void PoissonConv3DCons::computePerturbedPhysicalData(const Framework::State& state,
						 const RealVector& pdataBkp,
						 RealVector& pdata,
						 CFuint iVar)
{
  throw Common::NotImplementedException
      (FromHere(), "PoissonConv3DCons::computePerturbedPhysicalData() not implemented");
}

//////////////////////////////////////////////////////////////////////////////

bool PoissonConv3DCons::isValid(const RealVector& data)
{
  return true;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Poisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
