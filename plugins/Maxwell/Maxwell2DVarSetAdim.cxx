#include "Maxwell/Maxwell.hh"
#include "Maxwell2DVarSetAdim.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

Maxwell2DVarSetAdim::Maxwell2DVarSetAdim(Common::SafePtr<BaseTerm> term) :
  MaxwellVarSet(term),
  _model(term.d_castTo<MaxwellAdimTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Maxwell2DVarSetAdim::~Maxwell2DVarSetAdim()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal Maxwell2DVarSetAdim::getMaxEigenValue(const RealVector& data,
				       const RealVector& normal)
{
  return 1;  
}

//////////////////////////////////////////////////////////////////////////////

CFreal Maxwell2DVarSetAdim::getMaxAbsEigenValue(const RealVector& data,
					  const RealVector& normal)
{
  return 1;  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DVarSetAdim::computeEigenValues(const RealVector& data,
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

void Maxwell2DVarSetAdim::computeFlux (const RealVector& data,
				 const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];

  _fluxArray[0] = data[ConvMaxwellTerm::EZ]*ny;
  _fluxArray[1] = - data[ConvMaxwellTerm::EZ]*nx;
  _fluxArray[2] = data[ConvMaxwellTerm::EY]*nx-data[ConvMaxwellTerm::EX]*ny;
  _fluxArray[3] = - data[ConvMaxwellTerm::BZ]*ny;
  _fluxArray[4] = data[ConvMaxwellTerm::BZ]*nx;
  _fluxArray[5] = - data[ConvMaxwellTerm::BY]*nx + data[ConvMaxwellTerm::BX]*ny;
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DVarSetAdim::computeStateFlux (const RealVector& data)
{
  _physFlux(0,XX) = 0;
  _physFlux(0,YY) = data[ConvMaxwellTerm::EZ];

  _physFlux(1,XX) = - data[ConvMaxwellTerm::EZ];
  _physFlux(1,YY) = 0;

  _physFlux(2,XX) = data[ConvMaxwellTerm::EY];
  _physFlux(2,YY) = - data[ConvMaxwellTerm::EX];

  _physFlux(3,XX) = 0;
  _physFlux(3,YY) = - data[ConvMaxwellTerm::BZ];

  _physFlux(4,XX) = data[ConvMaxwellTerm::BZ];
  _physFlux(4,YY) = 0;

  _physFlux(5,XX) = - data[ConvMaxwellTerm::BY];
  _physFlux(5,YY) = data[ConvMaxwellTerm::BX];
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
