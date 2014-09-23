#include "Maxwell/Maxwell.hh"
#include "Maxwell3DVarSetAdim.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

Maxwell3DVarSetAdim::Maxwell3DVarSetAdim(Common::SafePtr<BaseTerm> term) :
  MaxwellVarSet(term),
  _model(term.d_castTo<MaxwellAdimTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Maxwell3DVarSetAdim::~Maxwell3DVarSetAdim()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal Maxwell3DVarSetAdim::getMaxEigenValue(const RealVector& data,
				       const RealVector& normal)
{
  return 1;  
}

//////////////////////////////////////////////////////////////////////////////

CFreal Maxwell3DVarSetAdim::getMaxAbsEigenValue(const RealVector& data,
					  const RealVector& normal)
{
  return 1;  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DVarSetAdim::computeEigenValues(const RealVector& data,
				       const RealVector& normal,
				       RealVector& result)
{
  result[0] = 0; 
  result[1] = 1;
  result[2] = 1;
  result[3] = 0;
  result[4] = -1;
  result[5] = -1;
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DVarSetAdim::computeFlux (const RealVector& data,
				 const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal nz = normals[ZZ];
 

  _fluxArray[0] = data[ConvMaxwellTerm::EZ]*ny - data[ConvMaxwellTerm::EY]*nz;
  _fluxArray[1] = - data[ConvMaxwellTerm::EZ]*nx + data[ConvMaxwellTerm::EX]*nz;
  _fluxArray[2] = data[ConvMaxwellTerm::EY]*nx - data[ConvMaxwellTerm::EX]*ny;
  _fluxArray[3] = - data[ConvMaxwellTerm::BZ]*ny + data[ConvMaxwellTerm::BY]*nz;
  _fluxArray[4] = data[ConvMaxwellTerm::BZ]*nx - data[ConvMaxwellTerm::BX]*nz;
  _fluxArray[5] = - data[ConvMaxwellTerm::BY]*nx + data[ConvMaxwellTerm::BX]*ny;
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DVarSetAdim::computeStateFlux (const RealVector& data)
{
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
  _physFlux(3,YY) = - data[ConvMaxwellTerm::BZ];
  _physFlux(3,ZZ) = data[ConvMaxwellTerm::BY];  

  _physFlux(4,XX) = data[ConvMaxwellTerm::BZ];
  _physFlux(4,YY) = 0;
  _physFlux(4,ZZ) = - data[ConvMaxwellTerm::BX];  

  _physFlux(5,XX) = - data[ConvMaxwellTerm::BY];
  _physFlux(5,YY) = data[ConvMaxwellTerm::BX];
  _physFlux(5,ZZ) = 0;  
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
