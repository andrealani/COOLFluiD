#include <numeric>
#include <algorithm>

#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "StrongSubOutletEuler2DCons.hh"
#include "CreateBoundaryNodalNormals.hh"
#include "InwardNormalsData.hh"
#include "NavierStokes/Euler2DCons.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongSubOutletEuler2DCons, FluctuationSplitData, FluctSplitNavierStokesModule> strongSubOutletEuler2DConsProvider("StrongSubOutletEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletEuler2DCons::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("P","static pressure");
}

//////////////////////////////////////////////////////////////////////////////

StrongSubOutletEuler2DCons::StrongSubOutletEuler2DCons(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_updateCoeff("updateCoeff"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _varSet(),
  _r3(),
  _bcNormals(),
  _dPdU()
{
   addConfigOptionsTo(this);
  _pressure = 0.0;
   setParameter("P",&_pressure);
}

//////////////////////////////////////////////////////////////////////////////

StrongSubOutletEuler2DCons::~StrongSubOutletEuler2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSubOutletEuler2DCons::needsSockets()
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

void StrongSubOutletEuler2DCons::setup()
{
  FluctuationSplitCom::setup();

  _bcNormals.resize(getTrsList().size());
  _varSet->setup();
  _r3.resize(PhysicalModelStack::getActive()->getNbEq());
  _dPdU.resize(PhysicalModelStack::getActive()->getNbEq());

  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
  obj.setDataSockets(socket_normals,socket_faceNeighCell);
  obj.create(getTrsList(), _bcNormals);

  _pressure /= _varSet->getModel()->getPressRef();

}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletEuler2DCons::executeOnTrs()
{
  vector<RealVector>* bcNormalsInTrs = &_bcNormals[getCurrentTrsID()];

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();


  Common::SafePtr< vector<CFuint> > statesIdx = getCurrentTRS()->getStatesInTrs();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    State *const state = states[stateID];
    const RealVector* bcNormal = &(*bcNormalsInTrs)[stateID];
    if (!isUpdated[stateID]) {
      const CFreal updateCoeffValue = updateCoeff[stateID];
      if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps()) {
        _varSet->setEigenVect3(_r3, *state, *bcNormal);
        const CFreal invRho = 1./(*state)[0];
        const CFreal rhoK2 = 0.5*((*state)[1]*(*state)[1] +
                            (*state)[2]*(*state)[2])*invRho;

        _dPdU[0] = gammaMinus1*invRho*rhoK2;
        _dPdU[1] = -gammaMinus1*invRho*(*state)[1];
        _dPdU[2] = -gammaMinus1*invRho*(*state)[2];
        _dPdU[3] = gammaMinus1;

        const CFreal p = gammaMinus1*((*state)[3] - rhoK2);
        const CFreal a = sqrt(gamma*p*invRho);
        CFreal *const rhsStart = &rhs(stateID, 0, nbEqs);
        const CFreal dPdUrhs = std::inner_product(&_dPdU[0],
                                             &_dPdU[0] + nbEqs,
                                             rhsStart, 0.0);
        const CFreal beta = (_pressure - p - dPdUrhs)*invRho/a;
        _r3 *= beta;

        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          rhs(stateID, iEq, nbEqs) += _r3[iEq];
        }
      }
      isUpdated[stateID] = true; // flagging is important!!!!!
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubOutletEuler2DCons::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  std::string varSetName = "Euler2DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler2DCons>());

  cf_assert(_varSet.isNotNull());

}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
