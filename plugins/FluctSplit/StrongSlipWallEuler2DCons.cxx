#include "Framework/CFL.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "NavierStokes/Euler2DCons.hh"

#include "FluctSplit/StrongSlipWallEuler2DCons.hh"
#include "FluctSplit/CreateBoundaryNodalNormals.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/FluctSplitNavierStokes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StrongSlipWallEuler2DCons,
                       FluctuationSplitData,
                       FluctSplitNavierStokesModule>
aStrongSlipWallEuler2DConsProvider("StrongSlipWallEuler2DCons");

//////////////////////////////////////////////////////////////////////////////

StrongSlipWallEuler2DCons::StrongSlipWallEuler2DCons(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _varSet(),
  _r3(),
  _bcNormals(),
  _useForInitialization(false)
{
}

//////////////////////////////////////////////////////////////////////////////

StrongSlipWallEuler2DCons::~StrongSlipWallEuler2DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSlipWallEuler2DCons::needsSockets()
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

void StrongSlipWallEuler2DCons::setup()
{

  FluctuationSplitCom::setup();

  _bcNormals.resize(getTrsList().size());
  _varSet->setup();
  _r3.resize(PhysicalModelStack::getActive()->getNbEq());

  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
  obj.setDataSockets(socket_normals,socket_faceNeighCell);
  obj.create(getTrsList(), _bcNormals);


}

//////////////////////////////////////////////////////////////////////////////

void StrongSlipWallEuler2DCons::executeOnTrs()
{
  vector<RealVector>* bcNormalsInTrs = &_bcNormals[getCurrentTrsID()];

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  Common::SafePtr< vector<CFuint> > const statesIdx = getCurrentTRS()->getStatesInTrs();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _useForInitialization = getMethodData().isInitializationPhase();
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    State *const state = states[stateID];

    const RealVector* bcNormal = &(*bcNormalsInTrs)[stateID];

    if(!_useForInitialization)
    {
      if (!isUpdated[stateID])
      {
        const CFreal updateCoeffValue = updateCoeff[stateID];
        if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps())
        {
          _varSet->setEigenVect3(_r3, *state, *bcNormal);

          const CFreal nx = (*bcNormal)[0];
          const CFreal ny = (*bcNormal)[1];
          const CFreal k = cfl / updateCoeffValue;

          const CFreal beta = -(((*state)[1] + k*rhs(stateID, 1, nbEqs))*nx +
                                ((*state)[2] + k*rhs(stateID, 2, nbEqs))*ny)/
                                 (k*_r3[1]*nx + k*_r3[2]*ny);

          _r3 *= beta;

          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
            rhs(stateID, iEq, nbEqs) += _r3[iEq];
          }
        }
        isUpdated[stateID] = true; // flagging is important!!!!!
      }
    }
    else
    {
      const CFreal nx = (*bcNormal)[0];
      const CFreal ny = (*bcNormal)[1];
      const CFreal rhoU = (*state)[1];
      const CFreal rhoV = (*state)[2];
      const CFreal origRhoVel = sqrt(rhoU*rhoU + rhoV*rhoV);
      const CFreal rhoVelNormal = rhoU*nx + rhoV*ny;
      (*state)[1] -= rhoVelNormal*nx;
      (*state)[2] -= rhoVelNormal*ny;
      const CFreal newRhoVel = sqrt((*state)[1]*(*state)[1] + (*state)[2]*(*state)[2]);
      const CFreal coeff = origRhoVel/newRhoVel;
      (*state)[1] *= coeff;
      (*state)[2] *= coeff;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongSlipWallEuler2DCons::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace(name);
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
