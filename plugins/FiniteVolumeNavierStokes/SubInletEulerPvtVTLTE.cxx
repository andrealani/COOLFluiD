#include "MathTools/MathConsts.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubInletEulerPvtVTLTE.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "Common/BadValueException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SubInletEulerPvtVTLTE, CellCenterFVMData, FiniteVolumeNavierStokesModule>
subInletEulerPvtVTLTEFVMCCProvider("SubInletEulerPvtVTLTEFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubInletEulerPvtVTLTE::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("UseOld","Use the old (buggy) implementation.");
}
      
//////////////////////////////////////////////////////////////////////////////

SubInletEulerPvtVTLTE::SubInletEulerPvtVTLTE(const std::string& name) :
  SubInletEulerFunc(name),
  m_varSet(CFNULL),
  m_library(CFNULL)
{
  addConfigOptionsTo(this);
  
  m_useOld = true;
  setParameter("UseOld",&m_useOld);
}
      
//////////////////////////////////////////////////////////////////////////////

SubInletEulerPvtVTLTE::~SubInletEulerPvtVTLTE()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEulerPvtVTLTE::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
    
  if(_useFunction)
  {
    // coordinate of the boundary point
    _bCoord = (innerState->getCoordinates() +
               ghostState->getCoordinates());
    _bCoord *= 0.5;
    
    // (*ghostState) = 2*bcState - (*innerState)
    _vFunction.evaluate(_bCoord, _inletData);
    
    _inletData[1] /= m_varSet->getModel()->getTempRef();
  }
  else {
    _inletData[0] = _massFlow;
    _inletData[1] = _temperature;
  }
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint TID = dim+1;
  (*ghostState)[0] = (*innerState)[0];
  (*ghostState)[1] = (*innerState)[1];
  (*ghostState)[2] = (*innerState)[2];
  if (dim == DIM_3D) {
    (*ghostState)[3] = (*innerState)[3];
  }
  (*ghostState)[TID] = _inletData[1];
  
  m_varSet->computePhysicalData(*ghostState, _dataInnerState);
  
  const CFreal area = MathTools::MathConsts::CFrealPi()*(_inletRadii[1]*_inletRadii[1] - _inletRadii[0]*_inletRadii[0]);
  
  if (m_useOld) {
    const CFreal rho = _dataInnerState[EulerTerm::RHO];
    // 0.001 is conversion from cm^2 to m^2
    const CFreal uInf = .001*_inletData[0]/(area*rho);
    
    (*ghostState)[0] = (*innerState)[0];
    (*ghostState)[1] = uInf;
    (*ghostState)[2] = 0.;
    if (dim == DIM_3D) {(*ghostState)[3] = 0.;}
    (*ghostState)[TID] = _inletData[1];
  }
  else {
    const CFreal pInner = _dataInnerState[EulerTerm::P];
    const CFreal Tinlet = _inletData[1];
    const CFreal RovM   = m_library->getRgas()/m_library->getMMass();
    const CFreal rho    = m_varSet->getModel()->getPressureFromState(pInner)/(RovM*Tinlet);
    // 0.001 is conversion from cm^2 to m^2
    const CFreal uInf   = .001*_inletData[0]/(area*rho);
    
    (*ghostState)[0] = (*innerState)[0];
    (*ghostState)[1] = 2.*uInf - (*innerState)[1];
    (*ghostState)[2] = 0.;
    if (dim == DIM_3D) {(*ghostState)[3] = 0.;}
    (*ghostState)[TID] = _inletData[1]; 
    const CFreal Tin = 2.*Tinlet - (*innerState)[TID];
    if (Tin > 0.) {
      (*ghostState)[TID] = Tin;
    }
    else {
      CFLog(VERBOSE, "SubInletEulerPvtVTLTE::setGhostState() => ghost T <= 0 => T = " << Tin << "\n");
    }  
  }
}

//////////////////////////////////////////////////////////////////////////////
      
void SubInletEulerPvtVTLTE::setup()
{
  SubInletEulerFunc::setup();
  
  m_varSet = getMethodData().getUpdateVar().d_castTo<EulerVarSet>();

  m_library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(m_library.isNotNull());
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
