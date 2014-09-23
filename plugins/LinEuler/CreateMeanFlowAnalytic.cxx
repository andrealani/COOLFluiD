#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"
#include "Framework/MethodCommandProvider.hh"
#include "LinEuler/LinearizedEuler.hh"
#include "LinEuler/CreateMeanFlowAnalytic.hh"
#include "LinEuler/LinEulerTerm.hh"
#include "Framework/SubSystemStatus.hh"

////////////////////////////////////////////////////////////////////////////
//

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

////////////////////////////////////////////////////////////////////////////
//

namespace COOLFluiD {

  namespace Physics {

    namespace LinearizedEuler {

////////////////////////////////////////////////////////////////////////////
//

MethodCommandProvider<CreateMeanFlowAnalytic,
                      DataProcessingData,
                      LinearizedEulerModule>
aCreateMeanFlowAnalyticProvider("CreateMeanFlowAnalytic");

////////////////////////////////////////////////////////////////////////////
//

void CreateMeanFlowAnalytic::defineConfigOptions(Config::OptionList&
options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Variable names.");
  options.addConfigOption< std::vector<std::string> >("MeanFlow","Function defining the mean flow.");
}

////////////////////////////////////////////////////////////////////////////
//

CreateMeanFlowAnalytic::CreateMeanFlowAnalytic(const std::string& name) :
  DataProcessingCom(name),
  socket_states("states"),
  socket_meanflow("meanflow")
{
  addConfigOptionsTo(this);

  m_vars_meanflow.resize(0);

  m_vars_meanflow = std::vector<std::string>();
  setParameter("Vars",&m_vars_meanflow);

  m_function_meanflow = vector<std::string>();
  setParameter("MeanFlow",&m_function_meanflow);

}

////////////////////////////////////////////////////////////////////////////
//

CreateMeanFlowAnalytic::~CreateMeanFlowAnalytic()
{
}

////////////////////////////////////////////////////////////////////////////
//

std::vector<Common::SafePtr<BaseDataSocketSink> >
CreateMeanFlowAnalytic::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  return result;
}

////////////////////////////////////////////////////////////////////////////
//
// This has to be deleted

std::vector<Common::SafePtr<BaseDataSocketSource> >
CreateMeanFlowAnalytic::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_meanflow);
  return result;
}
// This has to be deleted
////////////////////////////////////////////////////////////////////////////
//

void CreateMeanFlowAnalytic::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

  m_var_values.resize(m_vars_meanflow.size());

  if(m_function_meanflow.empty())
     throw BadValueException(FromHere(),"CreateMeanFlowAnalytic::setFuntion(): no mean flow functionprovided.");

  // configure the expression for the mean flow
  m_function_parser_meanflow.setFunctions(m_function_meanflow);
  m_function_parser_meanflow.setVariables(m_vars_meanflow);

  try
  {
    m_function_parser_meanflow.parse();
  }
  catch (Common::ParserException& e)
  {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

////////////////////////////////////////////////////////////////////////////
//

void CreateMeanFlowAnalytic::setup()
{
  DataHandle < Framework::State*, Framework::GLOBAL > states =
socket_states.getDataHandle();
// This has to be changed

  DataHandle<RealVector> meanflow = socket_meanflow.getDataHandle();
  meanflow.resize(states.size());
// This has to be changed


  const CFuint nb_funs = m_function_meanflow.size();

  for (CFuint i = 0; i < states.size(); ++i)
  {
    meanflow[i].resize(nb_funs);
  }

  executeOnTrs();

}

////////////////////////////////////////////////////////////////////////////
//

void CreateMeanFlowAnalytic::executeOnTrs()
{
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMin( "Meanflow::executeOnTrs() called for TRS: " <<
trs->getName() << "\n");

  DataHandle < Framework::State*, Framework::GLOBAL > states =
socket_states.getDataHandle();
// This has to be changed

  DataHandle<RealVector> meanflow = socket_meanflow.getDataHandle();
// This has to be changed


  const CFuint dim  = PhysicalModelStack::getActive()->getDim();
  const CFuint TT   = dim;
  const CFreal time =
SubSystemStatusStack::getActive()->getCurrentTimeDim();
  m_var_values[TT] = time;

  for (CFuint iState = 0; iState < states.size(); ++iState)
  {
    Node& coord = states[iState]->getCoordinates();
    for (CFuint iCoor = 0; iCoor < dim; ++iCoor)
    {
      m_var_values[iCoor] = coord[iCoor];
    }

    m_function_parser_meanflow.evaluate(m_var_values, meanflow[iState]);
  }

  SafePtr<LinEulerTerm> lterm = PhysicalModelStack::getActive()->getImplementor()-> getConvectiveTerm().d_castTo<LinEulerTerm>();

  lterm->setMeanFlowArray(meanflow);

}

////////////////////////////////////////////////////////////////////////////
//

void CreateMeanFlowAnalytic::unsetup()
{
  DataHandle<RealVector> meanflow = socket_meanflow.getDataHandle();

  meanflow.resize(0);
}

////////////////////////////////////////////////////////////////////////////
//

    } // namespace LinearizedEuler

  } // namespace Physics

} // namespace COOLFluiD
