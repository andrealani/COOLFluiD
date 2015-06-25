  #include "FluctSplit/FluctSplitNavierStokes.hh"
#include "StrongSlipWallEuler3DCons.hh"
#include "CreateBoundaryNodalNormals.hh"
#include "InwardNormalsData.hh"
#include "Framework/CFL.hh"
#include "NavierStokes/Euler3DCons.hh"
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

MethodCommandProvider<StrongSlipWallEuler3DCons, FluctuationSplitData, FluctSplitNavierStokesModule> strongSlipWallEuler3DConsProvider("StrongSlipWallEuler3DCons");

//////////////////////////////////////////////////////////////////////////////

StrongSlipWallEuler3DCons::StrongSlipWallEuler3DCons(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_updateCoeff("updateCoeff"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _varSet(),
  _r4(),
  _bcNormals(),
  _useForInitialization(false)
{
}

//////////////////////////////////////////////////////////////////////////////

StrongSlipWallEuler3DCons::~StrongSlipWallEuler3DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSlipWallEuler3DCons::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_states);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongSlipWallEuler3DCons::setup()
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

void StrongSlipWallEuler3DCons::executeOnTrs()
{
  vector<RealVector>* bcNormalsInTrs = &_bcNormals[getCurrentTrsID()];

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();

  Common::SafePtr< vector<CFuint> > const statesIdx = getCurrentTRS()->getStatesInTrs();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  _useForInitialization = getMethodData().isInitializationPhase();

  const CFreal CFL = getMethodData().getCFL()->getCFLValue();

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    State *const state = states[stateID];
    const RealVector* bcNormal = &(*bcNormalsInTrs)[stateID];

    if(!_useForInitialization)
    {
      if (!isUpdated[stateID]) {
        const CFreal updateCoeffValue = updateCoeff[stateID];
        if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps()) {
          _varSet->setEigenVect4(_r4, *state, *bcNormal);

          const CFreal nx = (*bcNormal)[XX];
          const CFreal ny = (*bcNormal)[YY];
          const CFreal nz = (*bcNormal)[ZZ];

          const CFreal k = CFL / updateCoeffValue;
          const CFreal beta = -(((*state)[1] + k*rhs(stateID, 1, nbEqs))*nx +
                                ((*state)[2] + k*rhs(stateID, 2, nbEqs))*ny +
                                ((*state)[3] + k*rhs(stateID, 3, nbEqs))*nz)/
            (k*_r4[1]*nx + k*_r4[2]*ny + k*_r4[3]*nz);

          _r4 *= beta;

          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
            rhs(stateID, iEq, nbEqs) += _r4[iEq];
          }
        }
        isUpdated[stateID] = true; // flagging is important!!!!!
      }
    }
    else
    {
      const CFreal nx = (*bcNormal)[0];
      const CFreal ny = (*bcNormal)[1];
      const CFreal nz = (*bcNormal)[2];

      const CFreal rhoU = (*state)[1];
      const CFreal rhoV = (*state)[2];
      const CFreal rhoW = (*state)[3];
      const CFreal origRhoVel = sqrt(rhoU*rhoU + rhoV*rhoV + rhoW*rhoW);
      const CFreal rhoVelNormal = rhoU*nx + rhoV*ny + rhoW*ny;
      (*state)[1] -= rhoVelNormal*nx;
      (*state)[2] -= rhoVelNormal*ny;
      (*state)[3] -= rhoVelNormal*nz;
      const CFreal newRhoVel = sqrt((*state)[1]*(*state)[1] + (*state)[2]*(*state)[2] + (*state)[3]*(*state)[3]);
      const CFreal coeff = origRhoVel/newRhoVel;
      (*state)[1] *= coeff;
      (*state)[2] *= coeff;
      (*state)[3] *= coeff;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongSlipWallEuler3DCons::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);
  
  const std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  const std::string varSetName = "Euler3DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
		 create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler3DCons>());
  
  cf_assert(_varSet.isNotNull());

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
