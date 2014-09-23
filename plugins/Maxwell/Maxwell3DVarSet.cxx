#include "Maxwell/Maxwell.hh"
#include "Maxwell3DVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

Maxwell3DVarSet::Maxwell3DVarSet(Common::SafePtr<BaseTerm> term) :
  MaxwellVarSet(term)
{
}

//////////////////////////////////////////////////////////////////////////////

Maxwell3DVarSet::~Maxwell3DVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DVarSet::setup()
{
  MaxwellVarSet::setup();

  // set EquationSetData
  Maxwell3DVarSet::getEqSetData().resize(1);
  Maxwell3DVarSet::getEqSetData()[0].setup(0,0,6);
}

//////////////////////////////////////////////////////////////////////////////

CFreal Maxwell3DVarSet::getMaxEigenValue(const RealVector& data,
				       const RealVector& normal)
{
  const CFreal ce = getModel()->getLightSpeed();

  return ce; 

}

//////////////////////////////////////////////////////////////////////////////

CFreal Maxwell3DVarSet::getMaxAbsEigenValue(const RealVector& data,
					  const RealVector& normal)
{
  const CFreal ce = getModel()->getLightSpeed();

  return ce;  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DVarSet::computeEigenValues(const RealVector& data,
				       const RealVector& normal,
				       RealVector& result)
{
  const CFreal ce = getModel()->getLightSpeed();

  result[0] = 0; 
  result[1] = ce;
  result[2] = - ce;
  result[3] = 0;
  result[4] = ce;
  result[5] = - ce;
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DVarSet::computeFlux (const RealVector& data,
				 const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal nz = normals[ZZ];
  const CFreal ce = getModel()->getLightSpeed();

  _fluxArray[0] = data[ConvMaxwellTerm::EZ]*ny - data[ConvMaxwellTerm::EY]*nz ;
  _fluxArray[1] = - data[ConvMaxwellTerm::EZ]*nx + data[ConvMaxwellTerm::EX]*nz;
  _fluxArray[2] = data[ConvMaxwellTerm::EY]*nx - data[ConvMaxwellTerm::EX]*ny;
  _fluxArray[3] = - ce*ce*data[ConvMaxwellTerm::BZ]*ny +  ce*ce*data[ConvMaxwellTerm::BY]*nz;
  _fluxArray[4] = ce*ce*data[ConvMaxwellTerm::BZ]*nx - ce*ce*data[ConvMaxwellTerm::BX]*nz ;
  _fluxArray[5] = - ce*ce*data[ConvMaxwellTerm::BY]*nx + ce*ce*data[ConvMaxwellTerm::BX]*ny;
 
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DVarSet::computeStateFlux (const RealVector& data)
{

  const CFreal ce = getModel()->getLightSpeed();
  
  _physFlux(0,XX) = 0;
  _physFlux(0,YY) = data[ConvMaxwellTerm::EZ];
  _physFlux(0,ZZ) = - data[ConvMaxwellTerm::EY]; 

  _physFlux(1,XX) = - data[ConvMaxwellTerm::EZ];
  _physFlux(1,YY) = 0;
  _physFlux(1,ZZ) = data[ConvMaxwellTerm::EX];  

  _physFlux(2,XX) = data[ConvMaxwellTerm::EY];
  _physFlux(2,YY) = - data[ConvMaxwellTerm::EX];
  _physFlux(2,ZZ) = 0;  

  _physFlux(3,XX) = 0;
  _physFlux(3,YY) = - ce*ce*data[ConvMaxwellTerm::BZ];
  _physFlux(3,ZZ) = ce*ce*data[ConvMaxwellTerm::BY];  

  _physFlux(4,XX) = ce*ce*data[ConvMaxwellTerm::BZ];
  _physFlux(4,YY) = 0;
  _physFlux(4,ZZ) = - ce*ce*data[ConvMaxwellTerm::BX];  

  _physFlux(5,XX) = - ce*ce*data[ConvMaxwellTerm::BY];
  _physFlux(5,YY) = ce*ce*data[ConvMaxwellTerm::BX];
  _physFlux(5,ZZ) = 0;  

  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
