#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "SuperOutletMHD2DProjection.hh"
#include "MHD/MHD2DProjectionVarSet.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SuperOutletMHD2DProjection, CellCenterFVMData, FiniteVolumeMHDModule> 
superOutletMHD2DProjectionFVMCCProvider("SuperOutletMHD2DProjectionFVMCC");
    
//////////////////////////////////////////////////////////////////////////////
   
void SuperOutletMHD2DProjection::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("refPhi","Reference phi value imposed at the supersonic outlet.");
}

//////////////////////////////////////////////////////////////////////////////

SuperOutletMHD2DProjection::SuperOutletMHD2DProjection(const std::string& name) : 
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);
  _refPhi = 0.;
   setParameter("refPhi",&_refPhi);
}

//////////////////////////////////////////////////////////////////////////////

SuperOutletMHD2DProjection::~SuperOutletMHD2DProjection() 
{
}

//////////////////////////////////////////////////////////////////////

void SuperOutletMHD2DProjection::setup()
{
 FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD2DProjectionVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletMHD2DProjection::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void SuperOutletMHD2DProjection::setGhostState(GeometricEntity *const face)
 {
   State *const innerState = face->getState(0);
   State *const ghostState = face->getState(1);
   // watch out NOT to use the operator=, because in that case 
   // the overloaded version operator=(State) would be used =>
   // also the coordinates (Node) would be set equal!!! 
   
   // set the physical data starting from the inner state
   _varSet->computePhysicalData(*innerState, _dataInnerState);

   _dataGhostState = _dataInnerState;

   // there are two possible boundary conditions for phi

   // 1) a reference value should be imposed for phi
   _dataGhostState[MHDProjectionTerm::PHI] = _refPhi;

   // 2)
   //_dataGhostState[MHDProjectionTerm::PHI] = -_dataInnerState[MHDProjectionTerm::PHI];

   _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
 }

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
