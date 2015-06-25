#include <numeric>

#include "FluctSplitLinEuler.hh"

#include "StrongSubOutletHedstrom2DCons.hh"
#include "FluctSplit/CreateBoundaryNodalNormals.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::LinearizedEuler;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongSubOutletHedstrom2DCons, FluctuationSplitData, FluctSplitLinEulerModule> StrongSubOutletHedstrom2DConsProvider("StrongSubOutletHedstrom2DCons");

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletHedstrom2DCons::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

StrongSubOutletHedstrom2DCons::StrongSubOutletHedstrom2DCons(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_updateCoeff("updateCoeff"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _varSet(),
  _r4(),
  _bcNormals()
{
   addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

StrongSubOutletHedstrom2DCons::~StrongSubOutletHedstrom2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSubOutletHedstrom2DCons::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletHedstrom2DCons::setup()
{
  FluctuationSplitCom::setup();

  _bcNormals.resize(getTrsList().size());
  _varSet->setup();
  _r4.resize(PhysicalModelStack::getActive()->getNbEq());

  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
  obj.setDataSockets(socket_normals,socket_faceNeighCell);
  obj.create(getTrsList(), _bcNormals);

}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletHedstrom2DCons::executeOnTrs()
{
  vector<RealVector>* bcNormalsInTrs = &_bcNormals[getCurrentTrsID()];

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();


  Common::SafePtr< vector<CFuint> > statesIdx = getCurrentTRS()->getStatesInTrs();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  const RealVector& linearData = _varSet->getModel()->getPhysicalData();
  const CFreal c     = linearData[LinEulerTerm::c];

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    State *const state = states[stateID];
    const RealVector* bcNormal = &(*bcNormalsInTrs)[stateID];
    const CFreal nx = (*bcNormal)[0];
    const CFreal ny = (*bcNormal)[1];

    if (!isUpdated[stateID]) {
      const CFreal updateCoeffValue = updateCoeff[stateID];
      if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps()) {

        _varSet->setEigenVect4(_r4, *state, *bcNormal);

        CFreal *const rhsStart = &rhs(stateID, 0, nbEqs);

        CFreal drho0u = rhsStart[1];
        CFreal drho0v = rhsStart[2];
        CFreal dp = rhsStart[3];

        const CFreal beta4 = -(nx*drho0u + ny*drho0v + dp/c);

        _r4 *= beta4;  // this is actually r3*beta3

        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          rhs(stateID, iEq, nbEqs) += _r4[iEq];
        }
      }
      isUpdated[stateID] = true; // flagging is important!!!!!
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletHedstrom2DCons::configure( Config::ConfigArgs& args)
{
  FluctuationSplitCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "LinEuler2DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::LinearizedEuler::LinEuler2DCons>());

  cf_assert(_varSet.isNotNull());

}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
