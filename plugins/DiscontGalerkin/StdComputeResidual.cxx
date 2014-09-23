#include "Framework/SubSystemStatus.hh"
#include "Common/COOLFluiD.hh"
#include "Framework/MethodCommandProvider.hh"

#include "DiscontGalerkin/StdComputeResidual.hh"
#include "DiscontGalerkin/DiscontGalerkin.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdComputeResidual,DiscontGalerkinSolverData,DiscontGalerkinModule >
  stdComputeResidualProvider("StdComputeResidual");

//////////////////////////////////////////////////////////////////////////////

StdComputeResidual::StdComputeResidual(const std::string& name)
  : DiscontGalerkinSolverCom(name),
    socket_rhs_norm("rhs_norm"),
    socket_states("states"),
    socket_old_states("old_states")
{
}

//////////////////////////////////////////////////////////////////////////////

StdComputeResidual::~StdComputeResidual()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdComputeResidual::setup()
{
  CFAUTOTRACE;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  DataHandle<CFreal> rhs_norm = socket_rhs_norm.getDataHandle();
  DataHandle<State*,GLOBAL> states = socket_states.getDataHandle();
  rhs_norm.resize(states.size()*nbEqs);
}

//////////////////////////////////////////////////////////////////////////////

void StdComputeResidual::unsetup()
{
  CFAUTOTRACE;
}
//////////////////////////////////////////////////////////////////////////////

void StdComputeResidual::execute()
{
  CFAUTOTRACE;
  DataHandle<State*> old_states = socket_old_states.getDataHandle();
  DataHandle<State*,GLOBAL> states = socket_states.getDataHandle();
  DataHandle<CFreal> rhs_norm = socket_rhs_norm.getDataHandle();
  const CFuint nbStates =states.size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  for(CFuint i=0; i<nbStates; ++i)
    for(CFuint j=0; j<nbEqs; ++j)
      rhs_norm[i*nbEqs + j] = (*states[i])[j] - (*old_states[i])[j];

  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  RealVector res(1);
  res[0]=getMethodData().getResidual();
  if (res[0]==0) res[0]=1;
  subSysStatus->setResidual(res);

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdComputeResidual::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_old_states);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
StdComputeResidual::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_rhs_norm);
  return result;
}

//////////////////////////////////////////////////////////////////////////////
  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

