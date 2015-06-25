#include <numeric>

#include "FluctSplitLinEuler.hh"

#include "StrongSubOutletHedstrom3DCons.hh"
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

MethodCommandProvider<StrongSubOutletHedstrom3DCons, FluctuationSplitData, FluctSplitLinEulerModule> StrongSubOutletHedstrom3DConsProvider("StrongSubOutletHedstrom3DCons");

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletHedstrom3DCons::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

StrongSubOutletHedstrom3DCons::StrongSubOutletHedstrom3DCons(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_updateCoeff("updateCoeff"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _varSet(),
  _r5(),
  _bcNormals()
{
   addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

StrongSubOutletHedstrom3DCons::~StrongSubOutletHedstrom3DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSubOutletHedstrom3DCons::needsSockets()
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

void StrongSubOutletHedstrom3DCons::setup()
{
  FluctuationSplitCom::setup();

  _bcNormals.resize(getTrsList().size());
  _varSet->setup();
  _r5.resize(PhysicalModelStack::getActive()->getNbEq());

  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
  obj.setDataSockets(socket_normals,socket_faceNeighCell);
  obj.create(getTrsList(), _bcNormals);

}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletHedstrom3DCons::executeOnTrs()
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
    const CFreal nz = (*bcNormal)[2];

    if (!isUpdated[stateID]) {
      const CFreal updateCoeffValue = updateCoeff[stateID];
      if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps()) {

        _varSet->setEigenVect5(_r5, *state, *bcNormal);

        CFreal *const rhsStart = &rhs(stateID, 0, nbEqs);

        CFreal drho0u = rhsStart[1];
        CFreal drho0v = rhsStart[2];
	CFreal drho0w = rhsStart[3];
        CFreal dp = rhsStart[4];

        const CFreal beta5 = -(-nx*drho0u - ny*drho0v - nz*drho0w + dp/c);

        _r5 *= beta5;

        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          rhs(stateID, iEq, nbEqs) += _r5[iEq];
        }
      }
      isUpdated[stateID] = true; // flagging is important!!!!!
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletHedstrom3DCons::configure( Config::ConfigArgs& args)
{
  FluctuationSplitCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "LinEuler3DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::LinearizedEuler::LinEuler3DCons>());

  cf_assert(_varSet.isNotNull());

}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
