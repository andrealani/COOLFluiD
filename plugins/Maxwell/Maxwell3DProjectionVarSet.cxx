#include "Maxwell3DProjectionVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

Maxwell3DProjectionVarSet::Maxwell3DProjectionVarSet(Common::SafePtr<BaseTerm> term) :
  MaxwellVarSet(term),
  _model(term.d_castTo<MaxwellProjectionTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Maxwell3DProjectionVarSet::~Maxwell3DProjectionVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal Maxwell3DProjectionVarSet::getMaxEigenValue(const RealVector& data,
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

CFreal Maxwell3DProjectionVarSet::getMaxAbsEigenValue(const RealVector& data,
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

void Maxwell3DProjectionVarSet::computeEigenValues(const RealVector& data,
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

void Maxwell3DProjectionVarSet::computeFlux (const RealVector& data,
				 const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal nz = normals[ZZ];

  const CFreal gamma = getModel()->getDivBCleaningConst();
  const CFreal ce = getModel()->getLightSpeed();
  const CFreal chi = getModel()->getDivECleaningConst();

  _fluxArray[0] = gamma*gamma*data[MaxwellProjectionTerm::PSI]*nx + data[ConvMaxwellTerm::EZ]*ny - data[ConvMaxwellTerm::EY]*nz ;
  _fluxArray[1] = - data[ConvMaxwellTerm::EZ]*nx + gamma*gamma*data[MaxwellProjectionTerm::PSI]*ny + data[ConvMaxwellTerm::EX]*nz;
  _fluxArray[2] = data[ConvMaxwellTerm::EY]*nx - data[ConvMaxwellTerm::EX]*ny + gamma*gamma*data[MaxwellProjectionTerm::PSI]*nz;
  _fluxArray[3] = chi*chi*ce*ce*data[MaxwellProjectionTerm::PHI]*nx - ce*ce*data[ConvMaxwellTerm::BZ]*ny +  ce*ce*data[ConvMaxwellTerm::BY]*nz;
  _fluxArray[4] = ce*ce*data[ConvMaxwellTerm::BZ]*nx + chi*chi*ce*ce*data[MaxwellProjectionTerm::PHI]*ny - ce*ce*data[ConvMaxwellTerm::BX]*nz ;
  _fluxArray[5] = - ce*ce*data[ConvMaxwellTerm::BY]*nx + ce*ce*data[ConvMaxwellTerm::BX]*ny + ce*ce*chi*chi*data[MaxwellProjectionTerm::PHI]*nz;
  _fluxArray[6] = (data[ConvMaxwellTerm::BX]*nx + data[ConvMaxwellTerm::BY]*ny + data[ConvMaxwellTerm::BZ]*nz)*ce*ce;
  _fluxArray[7] = data[ConvMaxwellTerm::EX]*nx + data[ConvMaxwellTerm::EY]*ny + data[ConvMaxwellTerm::EZ]*nz;
  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DProjectionVarSet::computeStateFlux (const RealVector& data)
{ 
  const CFreal ce = getModel()->getLightSpeed();
  const CFreal gamma = getModel()->getDivBCleaningConst();
  const CFreal chi = getModel()->getDivECleaningConst();  

  _physFlux(0,XX) = gamma*gamma*data[MaxwellProjectionTerm::PSI];
  _physFlux(0,YY) = data[ConvMaxwellTerm::EZ];
  _physFlux(0,ZZ) = - data[ConvMaxwellTerm::EY]; 

  _physFlux(1,XX) = - data[ConvMaxwellTerm::EZ];
  _physFlux(1,YY) = gamma*gamma*data[MaxwellProjectionTerm::PSI];
  _physFlux(1,ZZ) = data[ConvMaxwellTerm::EX];  

  _physFlux(2,XX) = data[ConvMaxwellTerm::EY];
  _physFlux(2,YY) = - data[ConvMaxwellTerm::EX];
  _physFlux(2,ZZ) = gamma*gamma*data[MaxwellProjectionTerm::PSI];  

  _physFlux(3,XX) = ce*ce*chi*chi*data[MaxwellProjectionTerm::PHI];
  _physFlux(3,YY) = - ce*ce*data[ConvMaxwellTerm::BZ];
  _physFlux(3,ZZ) = ce*ce*data[ConvMaxwellTerm::BY];  

  _physFlux(4,XX) = ce*ce*data[ConvMaxwellTerm::BZ];
  _physFlux(4,YY) = ce*ce*chi*chi*data[MaxwellProjectionTerm::PHI];
  _physFlux(4,ZZ) = - ce*ce*data[ConvMaxwellTerm::BX];  

  _physFlux(5,XX) = - ce*ce*data[ConvMaxwellTerm::BY];
  _physFlux(5,YY) = ce*ce*data[ConvMaxwellTerm::BX];
  _physFlux(5,ZZ) = ce*ce*chi*chi*data[MaxwellProjectionTerm::PHI];  
  
  _physFlux(6,XX) = data[ConvMaxwellTerm::BX]*ce*ce;
  _physFlux(6,YY) = data[ConvMaxwellTerm::BY]*ce*ce;
  _physFlux(6,ZZ) = data[ConvMaxwellTerm::BZ]*ce*ce;    

  _physFlux(7,XX) = data[ConvMaxwellTerm::EX];
  _physFlux(7,YY) = data[ConvMaxwellTerm::EY]; 
  _physFlux(8,ZZ) = data[ConvMaxwellTerm::EZ];  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
