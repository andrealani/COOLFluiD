#include "FluctSplit/FluctSplit.hh"
#include "StrongNoSlipWallAdiabaticNSKOmega3D.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/CFL.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongNoSlipWallAdiabaticNSKOmega3D, FluctuationSplitData, FluctSplitModule> strongNoSlipWallAdiabaticNS3DProvider("StrongNoSlipWallAdiabaticNSKOmega3D");

//////////////////////////////////////////////////////////////////////////////

StrongNoSlipWallAdiabaticNSKOmega3D::StrongNoSlipWallAdiabaticNSKOmega3D(const std::string& name) :
  FluctuationSplitCom(name),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_updateCoeff("updateCoeff"),
  socket_isUpdated("isUpdated"),
  _useForInitialization(false)
{
  addConfigOptionsTo(this);

  m_u0 = 0.0;
  setParameter("u0",&m_u0);
  m_v0 = 0.0;
  setParameter("v0",&m_v0);
  m_w0 = 0.0;
  setParameter("w0", &m_w0);


}

//////////////////////////////////////////////////////////////////////////////

StrongNoSlipWallAdiabaticNSKOmega3D::~StrongNoSlipWallAdiabaticNSKOmega3D()
{
}

//////////////////////////////////////////////////////////////////////////////
void StrongNoSlipWallAdiabaticNSKOmega3D::defineConfigOptions(Config::OptionList& options)
{
     options.addConfigOption< CFreal >("u0","x-coponent of the velocity of the wall");
     options.addConfigOption< CFreal >("v0","y-coponent of the velocity of the wall");
     options.addConfigOption< CFreal >("w0","z-coponent of the velocity of the wall");

}


//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongNoSlipWallAdiabaticNSKOmega3D::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_isUpdated);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongNoSlipWallAdiabaticNSKOmega3D::executeOnTrs()
{

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "StrongNoSlipWallAdiabaticNSKOmega3DImpl::execute() called for TRS: "
        << trs->getName() << "\n");

  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal cfl = getMethodData().getCFL()->getCFLValue();


  _useForInitialization = getMethodData().isInitializationPhase();

  Common::SafePtr< vector<CFuint> > const statesIdx = trs->getStatesInTrs();
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint localStateID = (*statesIdx)[iState];
    State *const currState = states[localStateID];

    if (!_useForInitialization) {
      if (currState->isParUpdatable()) {
	if (!isUpdated[localStateID]) {
          const CFreal coeff = max(updateCoeff[localStateID]/cfl, 1./cfl);

          // velocity has to remain null
          rhs(localStateID, 1, nbEqs) = 0.0;
          rhs(localStateID, 2, nbEqs) = 0.0;
	  rhs(localStateID, 3, nbEqs) = 0.0;
	  
	  rhs(localStateID, 5, nbEqs) = 0.0;
	  rhs(localStateID, 6, nbEqs) = 0.0;

          (*currState)[1] = m_u0;
          (*currState)[2] = m_v0;
          (*currState)[3] = m_w0;
	  
          (*currState)[5] = 0.0;
          (*currState)[6] = 0.0;	  


          }
          isUpdated[localStateID] = true; // flagging is important!!!!!
        }
      }
    
    else {
      (*currState)[1] = m_u0;
      (*currState)[2] = m_v0;
      (*currState)[3] = m_w0;
      (*currState)[5] = 0.0;
      (*currState)[6] = 0.0;
      
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
