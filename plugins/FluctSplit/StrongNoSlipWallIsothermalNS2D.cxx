#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "StrongNoSlipWallIsothermalNS2D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/CFL.hh"
#include "Framework/MeshData.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongNoSlipWallIsothermalNS2D, FluctuationSplitData, FluctSplitNavierStokesModule> strongNoSlipWallIsothermalNS2DProvider("StrongNoSlipWallIsothermalNS2D");

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallIsothermalNS2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >
    ("TWall","Wall temperature.");
}

//////////////////////////////////////////////////////////////////////////////

StrongNoSlipWallIsothermalNS2D::StrongNoSlipWallIsothermalNS2D(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  socket_isUpdated("isUpdated"),
  _useForInitialization(false)
{
   addConfigOptionsTo(this);

   _TWall = 0.0;
   setParameter("TWall",&_TWall);
}

//////////////////////////////////////////////////////////////////////////////

StrongNoSlipWallIsothermalNS2D::~StrongNoSlipWallIsothermalNS2D()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongNoSlipWallIsothermalNS2D::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_isUpdated);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallIsothermalNS2D::setup()
{

  FluctuationSplitCom::setup();

  cf_assert(_TWall > 0.);

  _TWall /= PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm().d_castTo<EulerTerm>()->getTempRef();
}

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallIsothermalNS2D::executeOnTrs()
{
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "StrongNoSlipWallAdiabaticNS2D::execute() called for TRS: "
        << trs->getName() << "\n");

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();

  SafePtr<EulerTerm> eulerModel = PhysicalModelStack::getActive()->
    getImplementor()->getConvectiveTerm().d_castTo<EulerTerm>();

  // (total) internal energy at the wall eWall = cv*Tw
  const CFreal gamma = eulerModel->getGamma();
  const CFreal eWall = eulerModel->getR()/(gamma - 1.)*_TWall;

  _useForInitialization = getMethodData().isInitializationPhase();

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint localStateID = (*statesIdx)[iState];
    State *const currState = states[localStateID];

    if (!_useForInitialization) {
      if (currState->isParUpdatable()) {
	if (!isUpdated[localStateID]) {
          // unused // const CFreal coeff = updateCoeff[localStateID]/cfl;
	  // velocity has to remain null
	  rhs(localStateID, 1, nbEqs) = 0.0;
	  rhs(localStateID, 2, nbEqs) = 0.0;
	  // d(RhoE) = Ewall*d(Rho)
	  // rhs(localStateID, 3, nbEqs) = coeff*rhs(localStateID, 0, nbEqs)*eWall;
	  rhs(localStateID, 3, nbEqs) = 0.0;

	  (*currState)[1] = 0.0;
	  (*currState)[2] = 0.0;
	  (*currState)[3] = (*currState)[0]*eWall;

	  isUpdated[localStateID] = true; // flagging is important!!!!!
	}
      }
    }
    else {
      (*currState)[1] = 0.0;
      (*currState)[2] = 0.0;
      (*currState)[3] = (*currState)[0]*eWall;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
