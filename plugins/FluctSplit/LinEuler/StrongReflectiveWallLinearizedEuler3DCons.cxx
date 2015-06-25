#include "FluctSplitLinEuler.hh"
#include <numeric>

#include "StrongReflectiveWallLinearizedEuler3DCons.hh"
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

MethodCommandProvider<StrongReflectiveWallLinearizedEuler3DCons, FluctuationSplitData, FluctSplitLinEulerModule> StrongReflectiveWallLinearizedEuler3DConsProvider("StrongReflectiveWallLinearizedEuler3DCons");

//////////////////////////////////////////////////////////////////////////////

StrongReflectiveWallLinearizedEuler3DCons::StrongReflectiveWallLinearizedEuler3DCons(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  socket_isUpdated("isUpdated"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  _varSet(),
  _bcNormals(),
  _r4()
{
}

//////////////////////////////////////////////////////////////////////////////

StrongReflectiveWallLinearizedEuler3DCons::~StrongReflectiveWallLinearizedEuler3DCons()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongReflectiveWallLinearizedEuler3DCons::needsSockets()
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

void StrongReflectiveWallLinearizedEuler3DCons::setup()
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

void StrongReflectiveWallLinearizedEuler3DCons::executeOnTrs()
{
  vector<RealVector>* bcNormalsInTrs = &_bcNormals[getCurrentTrsID()];
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  Common::SafePtr< vector<CFuint> > const statesIdx = getCurrentTRS()->getStatesInTrs();

  const CFuint m_nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
  RealVector Res(m_nbEqs);
  RealVector charRes(m_nbEqs);
  for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq) {
    Res[iEq] = 0.;
    charRes[iEq] = 0.;
  }

  const RealVector& linearData = _varSet->getModel()->getPhysicalData();
  const CFreal c     = linearData[LinEulerTerm::c];
  const CFreal oneoverc = 1./c;
  
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    State *const state = states[stateID];
    const RealVector* bcNormal = &(*bcNormalsInTrs)[stateID];
    const CFreal ncharx = -(*bcNormal)[0];
    const CFreal nchary = -(*bcNormal)[1];
    const CFreal ncharz = -(*bcNormal)[2];

    if (!isUpdated[stateID]) {
      const CFreal updateCoeffValue = updateCoeff[stateID];
      if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps()) {

//	transform the rhs vector to characteristic
      charRes[0] = -(rhs(stateID, 0, m_nbEqs) - rhs(stateID, 4, m_nbEqs)*oneoverc*oneoverc);
      
 if (nchary!=0.) { 
      charRes[1] = -(-nchary*rhs(stateID, 1, m_nbEqs) + ncharx*rhs(stateID, 2, m_nbEqs));
      charRes[2] = -(-ncharz*rhs(stateID, 2, m_nbEqs) + nchary*rhs(stateID, 3, m_nbEqs));
 }
 else if (ncharz!=0.) {
      charRes[1] = -(-ncharz*rhs(stateID, 1, m_nbEqs) + ncharx*rhs(stateID, 3, m_nbEqs));
      charRes[2] = -(-ncharz*rhs(stateID, 2, m_nbEqs) + nchary*rhs(stateID, 3, m_nbEqs));
 }
 else {
      charRes[1] = -(-ncharz*rhs(stateID, 1, m_nbEqs) + ncharx*rhs(stateID, 3, m_nbEqs));
      charRes[2] = -(-nchary*rhs(stateID, 1, m_nbEqs) + ncharx*rhs(stateID, 2, m_nbEqs));
 }
 
      charRes[3] = -(rhs(stateID, 1, m_nbEqs)*ncharx + rhs(stateID, 2, m_nbEqs)*nchary + rhs(stateID, 3, m_nbEqs)*ncharz + oneoverc*rhs(stateID, 4, m_nbEqs));

// apply the bc

      charRes[4] = charRes[3];


// transofrm back to conservative
      Res[0] = charRes[0]+0.5/c*(charRes[3]+charRes[4]);
      
           
 if (nchary!=0.) {      
      Res[1] = ((-(ncharz*ncharz+nchary*nchary)/nchary)*charRes[1] + (-(ncharx*ncharz)/nchary)*charRes[2] + 0.5*ncharx*(charRes[3]-charRes[4]));
      Res[2] = (ncharx*charRes[1] - ncharz*charRes[2] +0.5*nchary*(charRes[3]-charRes[4]));
      Res[3] = (ncharz*ncharx/nchary*charRes[1] + ((ncharx*ncharx+nchary*nchary)/nchary)*charRes[2] + 0.5*ncharz*(charRes[3]-charRes[4]));   
 }
 else if (ncharz!=0.) {
      Res[1] = ((-(ncharz*ncharz+nchary*nchary)/ncharz)*charRes[1] + ((ncharx*nchary)/ncharz)*charRes[2] + 0.5*ncharx*(charRes[3]-charRes[4]));
      Res[2] = (nchary*ncharx/ncharz*charRes[1] + (-(ncharx*ncharx+ncharz*ncharz)/ncharz)*charRes[2] + 0.5*nchary*(charRes[3]-charRes[4]));
      Res[3] = (ncharx*charRes[1] + nchary*charRes[2] + 0.5*ncharz*(charRes[3]-charRes[4]));
 }
 else {
      Res[1] = (-ncharz*charRes[1] -nchary*charRes[2] + 0.5*ncharx*(charRes[3]-charRes[4]));
      Res[2] = (-ncharz*nchary/ncharx*charRes[1] + ((ncharz*ncharz+ncharx*ncharx)/ncharx)*charRes[2] + 0.5*nchary*(charRes[3]-charRes[4])); 
      Res[3] =  (((ncharx*ncharx+nchary*nchary)/ncharx)*charRes[1] + (-(nchary*ncharz)/ncharx)*charRes[2] + 0.5*ncharz*(charRes[3]-charRes[4]));
 }
      
      Res[4] = (0.5*c*(charRes[3]+charRes[4]));
      
      
     for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq) {
       rhs(stateID, iEq, m_nbEqs) = -Res[iEq];
     }
      }
      isUpdated[stateID] = true; // flagging is important!!!!!
      
     
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StrongReflectiveWallLinearizedEuler3DCons::configure( Config::ConfigArgs& args)
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
