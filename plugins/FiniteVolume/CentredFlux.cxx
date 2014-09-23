#include "CentredFlux.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "Framework/MethodStrategyProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<CentredFlux,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeModule>
centredFluxSplitterProvider("Centred");
      
//////////////////////////////////////////////////////////////////////////////
      
CentredFlux::CentredFlux(const std::string& name) :
  FVMCC_FluxSplitter(name),
  _sumFlux(),
  _tempUnitNormal(),
  _statesLR(2)
{
}

//////////////////////////////////////////////////////////////////////////////

CentredFlux::~CentredFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

void CentredFlux::setup()
{
  FVMCC_FluxSplitter::setup();
  
  _sumFlux.resize(PhysicalModelStack::getActive()->getNbEq());
  _tempUnitNormal.resize(PhysicalModelStack::getActive()->getDim());
}
      
//////////////////////////////////////////////////////////////////////////////

void CentredFlux::compute(RealVector& result)
{
  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVar();
  vector<RealVector>& pdata = getMethodData().getPolyReconstructor()->getExtrapolatedPhysicaData();
  CellCenterFVMData& data = this->getMethodData(); 

  // flux for the right and left state
  _sumFlux = updateVarSet->getFlux()(pdata[1], data.getUnitNormal());
  _sumFlux += updateVarSet->getFlux()(pdata[0], data.getUnitNormal());
  result = 0.5*_sumFlux;
  
  // compute update coefficient
  if (!getMethodData().isPerturb()) {
    GeometricEntity& face = *data.getCurrentFace();
    const CFreal faceArea = socket_faceAreas.getDataHandle()[face.getID()]/
      data.getPolyReconstructor()->nbQPoints();
    
    DataHandle<CFreal> updateCoeff = this->socket_updateCoeff.getDataHandle();
    
    // left contribution to update coefficient
    const CFreal maxEV =
      updateVarSet->getMaxEigenValue(pdata[0], data.getUnitNormal());
    const CFuint leftID = face.getState(0)->getLocalID();
    updateCoeff[leftID] += max(maxEV, (CFreal)0.0)*faceArea;
    
    if (!face.getState(1)->isGhost()) {
      // right contribution to update coefficient
      _tempUnitNormal = -1.0*data.getUnitNormal();
      const CFreal maxEV =
	updateVarSet->getMaxEigenValue(pdata[1], _tempUnitNormal);
      const CFuint rightID = face.getState(1)->getLocalID();
      updateCoeff[rightID] += max(maxEV, (CFreal)0.0)*faceArea;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
