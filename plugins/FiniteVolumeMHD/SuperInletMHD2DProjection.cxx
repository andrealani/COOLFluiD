  #include "FiniteVolume/FiniteVolume.hh"
#include "SuperInletMHD2DProjection.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SuperInletMHD2DProjection, CellCenterFVMData, FiniteVolumeModule> superInletMHD2DProjectionFVMCCProvider("SuperInletMHD2DProjectionFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SuperInletMHD2DProjection::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("refPhi","Reference phi value imposed at the supersonic inlet.");
}

//////////////////////////////////////////////////////////////////////////////

SuperInletMHD2DProjection::SuperInletMHD2DProjection(const std::string& name) :
  SuperInlet(name)
{
   addConfigOptionsTo(this);

  _refPhi = 0.;
   setParameter("refPhi",&_refPhi);
}

//////////////////////////////////////////////////////////////////////////////

SuperInletMHD2DProjection::~SuperInletMHD2DProjection()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletMHD2DProjection::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // coordinate of the boundary point
  _bCoord = (innerState->getCoordinates() +
             ghostState->getCoordinates());
  _bCoord *= 0.5;

  // ghostState = 2*bcState - innerState
  _vFunction.evaluate(_bCoord, *_dimState);

  if(_inputAdimensionalValues)
  {
    *ghostState = *_dimState;
    (*ghostState)[8] = (*innerState)[8];
  }
  else
  {
    _varSet->setAdimensionalValues(*_dimState, *ghostState);
    (*ghostState)[8] = (*innerState)[8];
  }

  *ghostState *= 2.;
  *ghostState -= *innerState;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
