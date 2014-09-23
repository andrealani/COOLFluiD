#include "FiniteVolume/FiniteVolume.hh"
#include "SuperInletProjection.hh"
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

MethodCommandProvider<SuperInletProjection, CellCenterFVMData, FiniteVolumeModule> 
superInletProjectionFVMCCProvider("SuperInletProjectionFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SuperInletProjection::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("refPhi","Reference phi value imposed at the supersonic inlet.");
  options.addConfigOption< vector<CFuint> >("ProjectionIDs","IDs corresponding to projection fields.");
}

//////////////////////////////////////////////////////////////////////////////

SuperInletProjection::SuperInletProjection(const std::string& name) :
  SuperInlet(name)
{
  addConfigOptionsTo(this);

  _refPhi = 0.;
  setParameter("refPhi",&_refPhi);
  
  _projectionIDs = vector<CFuint>();
  setParameter("ProjectionIDs",&_projectionIDs);
}

//////////////////////////////////////////////////////////////////////////////

SuperInletProjection::~SuperInletProjection()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletProjection::setGhostState(GeometricEntity *const face)
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
  }
  else
  {
    _varSet->setAdimensionalValues(*_dimState, *ghostState);
  }
  
  *ghostState *= 2.;
  *ghostState -= *innerState;

  cf_assert(_projectionIDs.size() > 0);
  for (CFuint i = 0; i < _projectionIDs.size(); ++i) {
    const CFuint varID = _projectionIDs[i];
    (*ghostState)[varID] = (*innerState)[varID];
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
