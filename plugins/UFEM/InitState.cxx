#include "Common/CFLog.hh"

#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/NamespaceSwitcher.hh"

#include "UFEM/UFEM.hh"
#include "UFEM/InitState.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<InitState, UFEMSolverData, UFEMPlugin> initStateProvider("InitState");

//////////////////////////////////////////////////////////////////////////////

void InitState::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

InitState::InitState(const std::string& name) :
  UFEMSolverCom(name),
  socket_states("states"),
  socket_pastStates("pastStates"),
  socket_pastpastStates("pastpastStates")
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  m_func_str = std::vector<std::string>();
  setParameter("Def",&m_func_str);
}

//////////////////////////////////////////////////////////////////////////////

InitState::~InitState()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void InitState::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  UFEMSolverCom::configure(args);

  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance().getNamespace( getMethodData().getNamespace() );
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  const CFuint nbdim = physModel->getDim();

  std::vector<std::string> vars;
  if ( nbdim > DIM_0D ) vars.push_back ( "x" );
  if ( nbdim > DIM_1D ) vars.push_back ( "y" );
  if ( nbdim > DIM_2D ) vars.push_back ( "z" );
  cf_assert ( nbdim <= DIM_3D );
  vars.push_back ( "t" );

  m_init_func.setFunctions( m_func_str );
  m_init_func.setVariables( vars );
  m_init_func.parse();
}

//////////////////////////////////////////////////////////////////////////////

void InitState::executeOnTrs()
{
  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  CFLogDebugMin( "UFEM::InitState::executeOnTrs() called for TRS: " << trs->getName() << "\n");

cout << "UFEM::InitState::executeOnTrs" << endl << flush;

  DataHandle<State*,GLOBAL > states  = socket_states.getDataHandle();
  DataHandle<State*> pastStates = socket_pastStates.getDataHandle();
  DataHandle<State*> pastpastStates = socket_pastpastStates.getDataHandle();
  const CFreal dt=getMethodData().getDt();


  SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();
  std::vector<CFuint>::iterator itd;

  const CFuint nbdim = PhysicalModelStack::getActive()->getDim();
  const CFreal time = SubSystemStatusStack::getActive()->getCurrentTime();
  RealVector coord_t ( nbdim + 1 );
  for (itd = trsStates->begin(); itd != trsStates->end(); ++itd)
  {
    const RealVector& node = states[*itd]->getCoordinates();

    for ( CFuint i = 0; i < nbdim; ++i )
      coord_t[i] = node[i];

    coord_t[nbdim] = time;
    m_init_func.evaluate(coord_t,*states[*itd]);
    coord_t[nbdim] = time-dt;
    m_init_func.evaluate(coord_t,*pastStates[*itd]);
    coord_t[nbdim] = time-2.*dt;
    m_init_func.evaluate(coord_t,*pastpastStates[*itd]);
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
InitState::needsSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_states);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_pastpastStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace UFEM

} // namespace COOLFluiD
