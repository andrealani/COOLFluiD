#include "Maxwell2DProjectionVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

Maxwell2DProjectionVarSet::Maxwell2DProjectionVarSet(Common::SafePtr<BaseTerm> term) :
  MaxwellVarSet(term),
  _model(term.d_castTo<MaxwellProjectionTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Maxwell2DProjectionVarSet::~Maxwell2DProjectionVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal Maxwell2DProjectionVarSet::getMaxEigenValue(const RealVector& data,
					       const RealVector& normal)
{
  const CFreal ce = getModel()->getLightSpeed();
  const CFreal chi = getModel()->getDivECleaningConst();
  const CFreal gamma = getModel()->getDivBCleaningConst();
  
  if (gamma > 1){
    if (chi > gamma){
      return chi*ce;
    }
    else{
      return gamma*ce;
    }
  }
  else{
    if (chi > 1){
      return chi*ce;
    }
    else{
      return ce;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal Maxwell2DProjectionVarSet::getMaxAbsEigenValue(const RealVector& data,
					  const RealVector& normal)
{
  const CFreal ce = getModel()->getLightSpeed();
  const CFreal chi = getModel()->getDivECleaningConst();
  const CFreal gamma = getModel()->getDivBCleaningConst();
  
  if (gamma > 1.){
    if (chi > gamma){
      return chi*ce;
    }
    else{
      return gamma*ce;
    }
  }
  else{
    if (chi > 1.){
      return chi*ce;
    }
    else{
      return ce;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DProjectionVarSet::computeEigenValues(const RealVector& data,
				       const RealVector& normal,
				       RealVector& result)
{
  const CFreal ce = getModel()->getLightSpeed();
  const CFreal chi = getModel()->getDivECleaningConst();
  const CFreal gamma = getModel()->getDivBCleaningConst();

  result[0] = ce; 
  result[1] = ce;
  result[2] = -ce;
  result[3] = -ce;
  result[4] = gamma*ce;
  result[5] = -gamma*ce;
  result[6] = chi*ce;
  result[7] = -chi*ce;
  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DProjectionVarSet::computeFlux (const RealVector& data,
				 const RealVector& normals)
{
  CFLog(NOTICE, "computeFlux()\n");

  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal gamma = getModel()->getDivBCleaningConst();
  const CFreal ce = getModel()->getLightSpeed();
  const CFreal chi = getModel()->getDivECleaningConst();

  _fluxArray[0] = gamma*gamma*data[MaxwellProjectionTerm::PSI]*nx + data[ConvMaxwellTerm::EZ]*ny;
  _fluxArray[1] = - data[ConvMaxwellTerm::EZ]*nx + gamma*gamma*data[MaxwellProjectionTerm::PSI]*ny;
  _fluxArray[2] = data[ConvMaxwellTerm::EY]*nx - data[ConvMaxwellTerm::EX]*ny;
  _fluxArray[3] = (chi*chi*data[MaxwellProjectionTerm::PHI]*nx- data[ConvMaxwellTerm::BZ]*ny)*ce*ce;
  _fluxArray[4] = (data[ConvMaxwellTerm::BZ]*nx + chi*chi*data[MaxwellProjectionTerm::PHI]*ny)*ce*ce;
  _fluxArray[5] = (- data[ConvMaxwellTerm::BY]*nx + data[ConvMaxwellTerm::BX]*ny)*ce*ce;
  _fluxArray[6] = (data[ConvMaxwellTerm::BX]*nx + data[ConvMaxwellTerm::BY]*ny)*ce*ce;
  _fluxArray[7] = data[ConvMaxwellTerm::EX]*nx + data[ConvMaxwellTerm::EY]*ny;
  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DProjectionVarSet::computeStateFlux (const RealVector& data)
{
  
  const CFreal ce = getModel()->getLightSpeed();
  const CFreal gamma = getModel()->getDivBCleaningConst();
  const CFreal chi = getModel()->getDivECleaningConst();  


  _physFlux(0,XX) = gamma*gamma*data[MaxwellProjectionTerm::PSI];
  _physFlux(0,YY) = data[ConvMaxwellTerm::EZ];

  _physFlux(1,XX) = - data[ConvMaxwellTerm::EZ];
  _physFlux(1,YY) = gamma*gamma*data[MaxwellProjectionTerm::PSI];

  _physFlux(2,XX) = data[ConvMaxwellTerm::EY];
  _physFlux(2,YY) = - data[ConvMaxwellTerm::EX];

  _physFlux(3,XX) = chi*chi*data[MaxwellProjectionTerm::PHI]*ce*ce;
  _physFlux(3,YY) = - data[ConvMaxwellTerm::BZ]*ce*ce;

  _physFlux(4,XX) = data[ConvMaxwellTerm::BZ]*ce*ce;
  _physFlux(4,YY) = chi*chi*data[MaxwellProjectionTerm::PHI]*ce*ce;

  _physFlux(5,XX) = - data[ConvMaxwellTerm::BY]*ce*ce;
  _physFlux(5,YY) = data[ConvMaxwellTerm::BX]*ce*ce;
  
  _physFlux(6,XX) = data[ConvMaxwellTerm::BX]*ce*ce;
  _physFlux(6,YY) = data[ConvMaxwellTerm::BY]*ce*ce;

  _physFlux(7,XX) = data[ConvMaxwellTerm::EX];
  _physFlux(7,YY) = data[ConvMaxwellTerm::EY]; 
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
