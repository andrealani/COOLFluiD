#include "Maxwell2DProjectionAdimVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {

//////////////////////////////////////////////////////////////////////////////

Maxwell2DProjectionAdimVarSet::Maxwell2DProjectionAdimVarSet(Common::SafePtr<BaseTerm> term) :
  MaxwellVarSet(term),
  _model(term.d_castTo<MaxwellProjectionAdimTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Maxwell2DProjectionAdimVarSet::~Maxwell2DProjectionAdimVarSet()
{
}

//////////////////////////////////////////////////////////////////////////////

CFreal Maxwell2DProjectionAdimVarSet::getMaxEigenValue(const RealVector& data,
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

CFreal Maxwell2DProjectionAdimVarSet::getMaxAbsEigenValue(const RealVector& data,
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

void Maxwell2DProjectionAdimVarSet::computeEigenValues(const RealVector& data,
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

void Maxwell2DProjectionAdimVarSet::computeFlux (const RealVector& data,
				 const RealVector& normals)
{
  CFLog(NOTICE, "computeFlux()\n");

  const CFreal nx = normals[XX];
  const CFreal ny = normals[YY];
  const CFreal gamma = getModel()->getDivBAdimCleaningConst();

  const CFreal chi = getModel()->getDivECleaningConst();

  _fluxArray[0] = gamma*gamma*data[MaxwellProjectionAdimTerm::PSI]*nx + data[ConvMaxwellTerm::EZ]*ny;
  _fluxArray[1] = - data[ConvMaxwellTerm::EZ]*nx + gamma*gamma*data[MaxwellProjectionAdimTerm::PSI]*ny;
  _fluxArray[2] = data[ConvMaxwellTerm::EY]*nx - data[ConvMaxwellTerm::EX]*ny;
  _fluxArray[3] = chi*chi*data[MaxwellProjectionAdimTerm::PHI]*nx- data[ConvMaxwellTerm::BZ]*ny;
  _fluxArray[4] = data[ConvMaxwellTerm::BZ]*nx + chi*chi*data[MaxwellProjectionAdimTerm::PHI]*nx;
  _fluxArray[5] = - data[ConvMaxwellTerm::BY]*nx + data[ConvMaxwellTerm::BX]*ny;
  _fluxArray[6] = data[ConvMaxwellTerm::BX]*nx + data[ConvMaxwellTerm::BY]*ny;
  _fluxArray[7] = data[ConvMaxwellTerm::EX]*nx + data[ConvMaxwellTerm::EY]*ny;
  
}

//////////////////////////////////////////////////////////////////////////////

void Maxwell2DProjectionAdimVarSet::computeStateFlux (const RealVector& data)
{
  
  const CFreal gamma = getModel()->getDivBAdimCleaningConst();
  const CFreal chi = getModel()->getDivECleaningConst();  


  _physFlux(0,XX) = gamma*gamma*data[MaxwellProjectionAdimTerm::PSI];
  _physFlux(0,YY) = data[ConvMaxwellTerm::EZ];

  _physFlux(1,XX) = - data[ConvMaxwellTerm::EZ];
  _physFlux(1,YY) = gamma*gamma*data[MaxwellProjectionAdimTerm::PSI];

  _physFlux(2,XX) = data[ConvMaxwellTerm::EY];
  _physFlux(2,YY) = - data[ConvMaxwellTerm::EX];

  _physFlux(3,XX) = chi*chi*data[MaxwellProjectionAdimTerm::PHI];
  _physFlux(3,YY) = - data[ConvMaxwellTerm::BZ];

  _physFlux(4,XX) = data[ConvMaxwellTerm::BZ];
  _physFlux(4,YY) = chi*chi*data[MaxwellProjectionAdimTerm::PHI];

  _physFlux(5,XX) = - data[ConvMaxwellTerm::BY];
  _physFlux(5,YY) = data[ConvMaxwellTerm::BX];
  
  _physFlux(6,XX) = data[ConvMaxwellTerm::BX];
  _physFlux(6,YY) = data[ConvMaxwellTerm::BY];  

  _physFlux(7,XX) = data[ConvMaxwellTerm::EX];
  _physFlux(7,YY) = data[ConvMaxwellTerm::EY]; 
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Maxwell

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
