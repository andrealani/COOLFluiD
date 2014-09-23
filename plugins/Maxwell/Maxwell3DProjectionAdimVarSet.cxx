#include "Maxwell3DProjectionAdimVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

Maxwell3DProjectionAdimVarSet::Maxwell3DProjectionAdimVarSet(Common::SafePtr<BaseTerm> term) :
  MaxwellVarSet(term),
  _model(term.d_castTo<MaxwellProjectionAdimTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Maxwell3DProjectionAdimVarSet::~Maxwell3DProjectionAdimVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal Maxwell3DProjectionAdimVarSet::getMaxEigenValue(const RealVector& data,
					       const RealVector& normal)
{
  const CFreal chi = getModel()->getDivECleaningConst();
  const CFreal gamma = getModel()->getDivBAdimCleaningConst();
  
  if (gamma > 1){
    if (chi > gamma){
      return chi;
    }
    else{
      return gamma;
    }
  }
  else{
    if (chi > 1){
      return chi;
    }
    else{
      return 1;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal Maxwell3DProjectionAdimVarSet::getMaxAbsEigenValue(const RealVector& data,
					  const RealVector& normal)
{
  const CFreal chi = getModel()->getDivECleaningConst();
  const CFreal gamma = getModel()->getDivBAdimCleaningConst();
  
  if (gamma > 1){
    if (chi > gamma){
      return chi;
    }
    else{
      return gamma;
    }
  }
  else{
    if (chi > 1){
      return chi;
    }
    else{
      return 1;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DProjectionAdimVarSet::computeEigenValues(const RealVector& data,
				       const RealVector& normal,
				       RealVector& result)
{
  const CFreal gamma = getModel()->getDivBAdimCleaningConst();
  const CFreal chi = getModel()->getDivECleaningConst();
 
  result[0] = 1; 
  result[1] = 1;
  result[2] = -1;
  result[3] = -1;
  result[4] = gamma;
  result[5] = -gamma;
  result[6] = chi;
  result[7] = -chi;
  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DProjectionAdimVarSet::computeFlux (const RealVector& data,
				 const RealVector& normals)
{
  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal nz = normals[ZZ];
  const CFreal gamma = getModel()->getDivBAdimCleaningConst();
  const CFreal chi = getModel()->getDivECleaningConst();

  _fluxArray[0] = gamma*gamma*data[MaxwellProjectionAdimTerm::PSI]*nx + data[ConvMaxwellTerm::EZ]*ny - data[ConvMaxwellTerm::EY]*nz ;
  _fluxArray[1] = - data[ConvMaxwellTerm::EZ]*nx + gamma*gamma*data[MaxwellProjectionAdimTerm::PSI]*ny + data[ConvMaxwellTerm::EX]*nz;
  _fluxArray[2] = data[ConvMaxwellTerm::EY]*nx - data[ConvMaxwellTerm::EX]*ny + gamma*gamma*data[MaxwellProjectionAdimTerm::PSI]*nz;
  _fluxArray[3] = - data[ConvMaxwellTerm::BZ]*ny +  data[ConvMaxwellTerm::BY]*nz + chi*chi*data[MaxwellProjectionAdimTerm::PHI]*nx;
  _fluxArray[4] = data[ConvMaxwellTerm::BZ]*nx - data[ConvMaxwellTerm::BX]*nz + chi*chi*data[MaxwellProjectionAdimTerm::PHI]*ny;
  _fluxArray[5] = - data[ConvMaxwellTerm::BY]*nx + data[ConvMaxwellTerm::BX]*ny + chi*chi*data[MaxwellProjectionAdimTerm::PHI]*nz;
  _fluxArray[6] = data[ConvMaxwellTerm::BX]*nx + data[ConvMaxwellTerm::BY]*ny + data[ConvMaxwellTerm::BZ]*nz;
  _fluxArray[7] = data[ConvMaxwellTerm::EX]*nx + data[ConvMaxwellTerm::EY]*ny + data[ConvMaxwellTerm::EZ]*nz;
  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell3DProjectionAdimVarSet::computeStateFlux (const RealVector& data)
{ 
  const CFreal gamma = getModel()->getDivBAdimCleaningConst();  
  const CFreal chi = getModel()->getDivECleaningConst();  


  _physFlux(0,XX) = gamma*gamma*data[MaxwellProjectionAdimTerm::PSI];
  _physFlux(0,YY) = data[ConvMaxwellTerm::EZ];
  _physFlux(0,ZZ) = - data[ConvMaxwellTerm::EY]; 

  _physFlux(1,XX) = - data[ConvMaxwellTerm::EZ];
  _physFlux(1,YY) = gamma*gamma*data[MaxwellProjectionAdimTerm::PSI];
  _physFlux(1,ZZ) = data[ConvMaxwellTerm::EX];  

  _physFlux(2,XX) = data[ConvMaxwellTerm::EY];
  _physFlux(2,YY) = - data[ConvMaxwellTerm::EX];
  _physFlux(2,ZZ) = gamma*gamma*data[MaxwellProjectionAdimTerm::PSI];  

  _physFlux(3,XX) = chi*chi*data[MaxwellProjectionAdimTerm::PHI];
  _physFlux(3,YY) = - data[ConvMaxwellTerm::BZ];
  _physFlux(3,ZZ) = data[ConvMaxwellTerm::BY];  

  _physFlux(4,XX) = data[ConvMaxwellTerm::BZ];
  _physFlux(4,YY) = chi*chi*data[MaxwellProjectionAdimTerm::PHI];
  _physFlux(4,ZZ) = - data[ConvMaxwellTerm::BX];  

  _physFlux(5,XX) = - data[ConvMaxwellTerm::BY];
  _physFlux(5,YY) = data[ConvMaxwellTerm::BX];
  _physFlux(5,ZZ) = chi*chi*data[MaxwellProjectionAdimTerm::PHI];  
  
  _physFlux(6,XX) = data[ConvMaxwellTerm::BX];
  _physFlux(6,YY) = data[ConvMaxwellTerm::BY];
  _physFlux(6,ZZ) = data[ConvMaxwellTerm::BZ];    

  _physFlux(7,XX) = data[ConvMaxwellTerm::EX];
  _physFlux(7,YY) = data[ConvMaxwellTerm::EY]; 
  _physFlux(8,ZZ) = data[ConvMaxwellTerm::EZ];  
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
