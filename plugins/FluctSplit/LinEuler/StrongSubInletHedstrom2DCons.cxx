#include "FluctSplitLinEuler.hh"
#include <numeric>

#include "StrongSubInletHedstrom2DCons.hh"
#include "FluctSplit/CreateBoundaryNodalNormals.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "MathTools/MathFunctions.hh"

#include "LinEuler/LinearizedEuler.hh"
#include "LinEuler/LinEuler2DVarSet.hh"
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::LinearizedEuler;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongSubInletHedstrom2DCons, FluctuationSplitData, FluctSplitLinEulerModule> StrongSubInletHedstrom2DConsProvider("StrongSubInletHedstrom2DCons");

//////////////////////////////////////////////////////////////////////////////

StrongSubInletHedstrom2DCons::StrongSubInletHedstrom2DCons(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _varSet(),
  _bcNormals(),
  _r1(),
  _r2(),
  _r3()
{
}

//////////////////////////////////////////////////////////////////////////////

StrongSubInletHedstrom2DCons::~StrongSubInletHedstrom2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSubInletHedstrom2DCons::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletHedstrom2DCons::setup()
{

  FluctuationSplitCom::setup();

  _bcNormals.resize(getTrsList().size());
  _varSet->setup();
  _r1.resize(PhysicalModelStack::getActive()->getNbEq());
  _r2.resize(PhysicalModelStack::getActive()->getNbEq());
  _r3.resize(PhysicalModelStack::getActive()->getNbEq());

  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
  obj.setDataSockets(socket_normals,socket_faceNeighCell);
  obj.create(getTrsList(), _bcNormals);

}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletHedstrom2DCons::executeOnTrs()
{
  vector<RealVector>* bcNormalsInTrs = &_bcNormals[getCurrentTrsID()];
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  Common::SafePtr< vector<CFuint> > const statesIdx = getCurrentTRS()->getStatesInTrs();

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
        _varSet->setEigenVect1(_r1, *state, *bcNormal);
        _varSet->setEigenVect2(_r2, *state, *bcNormal);
        _varSet->setEigenVect3(_r3, *state, *bcNormal);

        CFreal *const rhsStart = &rhs(stateID, 0, nbEqs);

        CFreal drho = rhsStart[0];
        CFreal drho0u = rhsStart[1];
        CFreal drho0v = rhsStart[2];
        CFreal dp = rhsStart[3];

        const CFreal beta1 = -(drho - dp/(c*c));
        const CFreal beta2 = -(ny*drho0u - nx*drho0v);
        const CFreal beta4 = -(-nx*drho0u - ny*drho0v + dp/c);

        _r1 *= beta1;
        _r2 *= beta2;
        _r3 *= beta4;
        _r3 += _r2 += _r1; // This is actually the whole boundary correction term

        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          rhs(stateID, iEq, nbEqs) += _r3[iEq];
        }
      }
      isUpdated[stateID] = true; // flagging is important!!!!!
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletHedstrom2DCons::configure( Config::ConfigArgs& args)
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
