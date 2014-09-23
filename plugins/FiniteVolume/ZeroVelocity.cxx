#include "FiniteVolume/FiniteVolume.hh"
#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolume/ZeroVelocity.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ZeroVelocity, CellCenterFVMData,
		      FiniteVolumeModule>
zeroVelocity("ZeroVelocityFVMCC");

//////////////////////////////////////////////////////////////////////////////

ZeroVelocity::ZeroVelocity(const std::string& name) :
  FVMCC_BC(name),
  _isVelocityComp()
{
  addConfigOptionsTo(this);
  _velocityIDs = vector<CFuint>();
  setParameter("VelocityIDs",&_velocityIDs);
}

//////////////////////////////////////////////////////////////////////////////

void ZeroVelocity::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< vector<CFuint> >
    ("VelocityIDs", "Variable IDs of the velocity components");
}

//////////////////////////////////////////////////////////////////////////////

ZeroVelocity::~ZeroVelocity()
{
}

//////////////////////////////////////////////////////////////////////////////

void ZeroVelocity::setup()
{
  FVMCC_BC::setup();

  if(_velocityIDs.size() == 0) {
    CFLog(NOTICE, "ZeroVelocity::setup() => choosing default\n");
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    _velocityIDs.resize(dim);
    for (CFuint i = 0 ; i < dim; ++i) {
      _velocityIDs[i] = 1 + i;
    }
  }

  _isVelocityComp.resize(PhysicalModelStack::getActive()->getNbEq());
  _isVelocityComp = false;
  for (CFuint i = 0 ; i < _velocityIDs.size(); ++i) {
    _isVelocityComp[_velocityIDs[i]] = true;
  }
}

//////////////////////////////////////////////////////////////////////////////

void ZeroVelocity ::setGhostState
(Framework::GeometricEntity *const face)
{
  using namespace COOLFluiD::Framework;

  // unused // const CFuint dim = PhysicalModelStack::getActive()->getDim();
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  for (CFuint i = 0; i < _isVelocityComp.size(); ++i) {
    if (!_isVelocityComp[i]) {
      (*ghostState)[i] = (*innerState)[i];
    }
    else {
      (*ghostState)[i] = -(*innerState)[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
