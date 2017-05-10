#include "Maxwell/Maxwell.hh"
#include "Maxwell2DVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

Maxwell2DVarSet::Maxwell2DVarSet(Common::SafePtr<BaseTerm> term) :
  MaxwellVarSet(term)
{
}

//////////////////////////////////////////////////////////////////////////////

Maxwell2DVarSet::~Maxwell2DVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DVarSet::setup()
{  
  MaxwellVarSet::setup();
  
  // set EquationSetData
  Maxwell2DVarSet::getEqSetData().resize(1);
  Maxwell2DVarSet::getEqSetData()[0].setup(0,0,6);
}

//////////////////////////////////////////////////////////////////////////////

CFreal Maxwell2DVarSet::getMaxEigenValue(const RealVector& data,
				       const RealVector& normal)
{
  
  const CFreal ce = getModel()->getLightSpeed();
  return ce;  
}

//////////////////////////////////////////////////////////////////////////////

CFreal Maxwell2DVarSet::getMaxAbsEigenValue(const RealVector& data,
					  const RealVector& normal)
{
  const CFreal ce = getModel()->getLightSpeed();

  return ce;  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DVarSet::computeEigenValues(const RealVector& data,
				       const RealVector& normal,
				       RealVector& result)
{
  const CFreal ce = getModel()->getLightSpeed();
  const vector<CFuint>& varIDs = Maxwell2DVarSet::getEqSetData()[0].getEqSetVarIDs();
  
  result[varIDs[0]] = 0; 
  result[varIDs[1]] = ce;
  result[varIDs[2]] = - ce;
  result[varIDs[3]] = 0;
  result[varIDs[4]] = ce;
  result[varIDs[5]] = - ce;
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DVarSet::computeFlux (const RealVector& data,
				   const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal ce = getModel()->getLightSpeed();
  const vector<CFuint>& varIDs = Maxwell2DVarSet::getEqSetData()[0].getEqSetVarIDs();
  
  _fluxArray[varIDs[0]] = data[ConvMaxwellTerm::EZ]*ny;
  _fluxArray[varIDs[1]] = - data[ConvMaxwellTerm::EZ]*nx;
  _fluxArray[varIDs[2]] = data[ConvMaxwellTerm::EY]*nx-data[ConvMaxwellTerm::EX]*ny;
  _fluxArray[varIDs[3]] = - ce*ce*data[ConvMaxwellTerm::BZ]*ny;
  _fluxArray[varIDs[4]] = ce*ce*data[ConvMaxwellTerm::BZ]*nx;
  _fluxArray[varIDs[5]] = - ce*ce*data[ConvMaxwellTerm::BY]*nx + ce*ce*data[ConvMaxwellTerm::BX]*ny;
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DVarSet::computeStateFlux (const RealVector& data)
{
  
  const CFreal ce = getModel()->getLightSpeed();


  _physFlux(0,XX) = 0;
  _physFlux(0,YY) = data[ConvMaxwellTerm::EZ];

  _physFlux(1,XX) = - data[ConvMaxwellTerm::EZ];
  _physFlux(1,YY) = 0;

  _physFlux(2,XX) = data[ConvMaxwellTerm::EY];
  _physFlux(2,YY) = - data[ConvMaxwellTerm::EX];

  _physFlux(3,XX) = 0;
  _physFlux(3,YY) = - ce*ce*data[ConvMaxwellTerm::BZ];

  _physFlux(4,XX) = ce*ce*data[ConvMaxwellTerm::BZ];
  _physFlux(4,YY) = 0;

  _physFlux(5,XX) = - ce*ce*data[ConvMaxwellTerm::BY];
  _physFlux(5,YY) = ce*ce*data[ConvMaxwellTerm::BX];
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
