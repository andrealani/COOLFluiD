#include "HyperPoisson/HyperPoisson.hh"
#include "HyperPoisson3DCons.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace HyperPoisson {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<HyperPoisson3DCons, ConvectiveVarSet, HyperPoissonModule, 1>
HyperPoisson3DConsProvider("HyperPoisson3DCons");

//////////////////////////////////////////////////////////////////////////////

HyperPoisson3DCons::HyperPoisson3DCons(Common::SafePtr<BaseTerm> term) :
  HyperPoisson3DVarSet(term),
  _rightEv(4,4),
  _leftEv(4,4)
{
  vector<std::string> names(4);
  names[0] = "phi";
  names[1] = "Bx";
  names[2] = "By";
  names[3] = "Bz";
  setVarNames(names);
}


//////////////////////////////////////////////////////////////////////////////

HyperPoisson3DCons::~HyperPoisson3DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DCons::setup()
{
  HyperPoisson3DVarSet::setup();

  setConstJacob();
}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DCons::setConstJacob()
{
  vector<RealMatrix>* const jacobians = PhysicalModelStack::getActive()->getImplementor()->getJacobians();

  (*jacobians)[0](0,1) = 1.0;
  (*jacobians)[0](1,0) = 1.0;

  (*jacobians)[1](0,2) = 1.0;
  (*jacobians)[1](2,0) = 1.0;

  (*jacobians)[2](0,3) = 1.0;
  (*jacobians)[2](3,0) = 1.0;
}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DCons::computeJacobians()
{

}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DCons::computeEigenValuesVectors(RealMatrix& rightEv,
                                        RealMatrix& leftEv,
                                        RealVector& eValues,
                                        const RealVector& normal)
{

}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DCons::splitJacobian(RealMatrix& jacobPlus,
                             RealMatrix& jacobMin,
                             RealVector& eValues,
                             const RealVector& normal)
{

}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DCons::setEigenVect1(RealVector& r1,
                                State& state,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"HyperPoisson3DCons::setEigenVect1()");
}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DCons::setEigenVect2(RealVector& r2,
                                State& state,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"HyperPoisson3DCons::setEigenVect2()");
}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DCons::setEigenVect3(RealVector& r3,
                                State& state,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"HyperPoisson3DCons::setEigenVect3()");
}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DCons::setEigenVect4(RealVector& r4,
                                State& state,
                                const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"HyperPoisson3DCons::setEigenVect4()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint HyperPoisson3DCons::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DCons::computePhysicalData(const State& state, RealVector& data)
{
  const CFreal phi  = state[0];
  const CFreal bx = state[1];
  const CFreal by = state[2];
  const CFreal bz = state[3];

  data[HyperPTerm::PHI] = phi;
  data[HyperPTerm::BX] = bx;
  data[HyperPTerm::BY] = by;
  data[HyperPTerm::BZ] = bz;
}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DCons::computeStateFromPhysicalData(const RealVector& data,
					  State& state)
{
  const CFreal phi = data[HyperPTerm::PHI];
  state[0] = phi;
  state[1] = data[HyperPTerm::BX];
  state[2] = data[HyperPTerm::BY];
  state[3] = data[HyperPTerm::BZ];
}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DCons::setDimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData =
    getModel()->getReferencePhysicalData();

  result[0] = state[0];
  result[1] = state[1];
  result[2] = state[2];
  result[3] = state[3];
}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DCons::setAdimensionalValues(const State& state,
                                       RealVector& result)
{
  const RealVector& refData =
    getModel()->getReferencePhysicalData();

  result[0] = state[0];
  result[1] = state[1];
  result[2] = state[2];
  result[3] = state[3];
}

//////////////////////////////////////////////////////////////////////////////
      
void HyperPoisson3DCons::computeProjectedJacobian(const RealVector& normal,
					   RealMatrix& jacob)
{

}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace HyperPoisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
