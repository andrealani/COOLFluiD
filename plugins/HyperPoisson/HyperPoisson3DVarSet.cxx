#include "HyperPoisson/HyperPoisson.hh"
#include "HyperPoisson3DVarSet.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/CFPrintContainer.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace HyperPoisson {

//////////////////////////////////////////////////////////////////////////////

HyperPoisson3DVarSet::HyperPoisson3DVarSet(Common::SafePtr<BaseTerm> term) :
  HyperPoissonVarSet(term)
{
}

//////////////////////////////////////////////////////////////////////////////

HyperPoisson3DVarSet::~HyperPoisson3DVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DVarSet::setup()
{
  HyperPoissonVarSet::setup();

  // set EquationSetData
  HyperPoisson3DVarSet::getEqSetData().resize(1);
  HyperPoisson3DVarSet::getEqSetData()[0].setup(0,0,4);
}

//////////////////////////////////////////////////////////////////////////////

CFreal HyperPoisson3DVarSet::getMaxEigenValue(const RealVector& data,
				       const RealVector& normal)
{
  return 1;
}

//////////////////////////////////////////////////////////////////////////////

CFreal HyperPoisson3DVarSet::getMaxAbsEigenValue(const RealVector& data,
					  const RealVector& normal)
{
  return 1;
}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DVarSet::computeEigenValues(const RealVector& data,
				       const RealVector& normal,
				       RealVector& result)
{
  //const CFreal un = getNormalSpeed(data, normal);
  const vector<CFuint>& varIDs = HyperPoisson3DVarSet::getEqSetData()[0].getEqSetVarIDs();
  
    result[varIDs[0]] = 0.;
    result[varIDs[1]] = 1.;
    result[varIDs[2]] = 1.;
    result[varIDs[3]] = 1.;
  
}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DVarSet::computeFlux (const RealVector& data,
				 const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal nz = (normals.size() == 3) ? normals[ZZ] : 0.;
  const CFreal Bx = data[HyperPTerm::BX];
  const CFreal By = data[HyperPTerm::BY];
  const CFreal Bz = data[HyperPTerm::BZ];
  const CFreal Bn = (Bx*nx + By*ny + Bz*nz);
  const vector<CFuint>& varIDs = HyperPoisson3DVarSet::getEqSetData()[0].getEqSetVarIDs();
  
  //CFLog(DEBUG_MIN, "HyperPoisson3DVarSet::computeFlux() => " << CFPrintContainer<const vector<CFuint> >("varIDs = ", &varIDs));
  
    _fluxArray[varIDs[0]] = -Bn;
    _fluxArray[varIDs[1]] = -(data[HyperPTerm::PHI])*nx ;
    _fluxArray[varIDs[2]] = -(data[HyperPTerm::PHI])*ny ;
    _fluxArray[varIDs[3]] = -(data[HyperPTerm::PHI])*nz ;
}

//////////////////////////////////////////////////////////////////////////////

void HyperPoisson3DVarSet::computeStateFlux (const RealVector& data)
{

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace HyperPoisson

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
