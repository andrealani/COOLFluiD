#include "Environment/DirPaths.hh"

#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PathAppender.hh"

#include "AnalyticalEE/AnalyticalEE.hh"
#include "AnalyticalEE/ComputeDiscreteError.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AnalyticalEE {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeDiscreteError, AnalyticEEData, AnalyticalEEModule> computeDiscreteErrorProvider("ComputeDiscreteError");

//////////////////////////////////////////////////////////////////////////////

void ComputeDiscreteError::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Functions","Definition of the Functions.");
  options.addConfigOption< std::string >("OutputFile","Name of Output File to write the Error Norm.");
}

//////////////////////////////////////////////////////////////////////////////

ComputeDiscreteError::ComputeDiscreteError(std::string name) : AnalyticEECom(name),
  socket_nodes("nodes"),
  socket_states("states"),
  m_L2(0),
  m_result(0)
{
  m_file = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  addConfigOptionsTo(this);

  m_vars = std::vector<std::string>();
  setParameter("Vars",&m_vars);

//   m_vars[0] = "x";
//   m_vars[1] = "y";

  m_functions = std::vector<std::string>();
  setParameter("Functions",&m_functions);

  m_nameOutputFile = "DiscreteErrorNorm.plt";
  setParameter("OutputFile",&m_nameOutputFile);

}

//////////////////////////////////////////////////////////////////////////////

ComputeDiscreteError::~ComputeDiscreteError()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiscreteError::setup()
{
  CFAUTOTRACE;

  // first call parent method
  AnalyticEECom::setup();

  const CFuint nbeqs = PhysicalModelStack::getActive()->getNbEq();

  m_exact.resize(nbeqs);
  m_error.resize(nbeqs);

  m_L2.resize(nbeqs);
  m_result.resize(nbeqs);
  prepareOutputFile();

}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiscreteError::unsetup()
{
  CFAUTOTRACE;

  m_exact.resize(0);
  m_error.resize(0);
  m_L2.resize(0);
  m_result.resize(0);
  m_file->close();

  // last call parent method
  AnalyticEECom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiscreteError::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  AnalyticEECom::configure(args);

  m_vFunction.setFunctions(m_functions);


//   const CFuint dim = PhysicalModelStack::getActive()->getDim();
//   const CFuint dim = 2;
//
//   m_vars = std::vector<std::string>();
//   m_vars.resize(dim);
//   m_vars[0] = "x";
//   m_vars[1] = "y";
//
//   if(dim == DIM_3D) m_vars[2] = "z";
//
//   CFout << "These are the variables:\n";
//   CFout << m_vars[0] << "\n";
//   CFout << m_vars[1] << "\n";
//   CFout << m_vars[2] << "\n";


  m_vFunction.setVariables(m_vars);
  try {
    m_vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }


//   CFout << "Compute discrete error:\n";
//   CFout << "Functions: " << m_functions[0] << "\n";
//   CFout << "Variables:\n";
//   CFout << "\t\t" << m_vars[0] << "\n";
//   CFout << "\t\t" << m_vars[1] << "\n";
//   CFout << "Size of variables' field: " << m_vars.size() << "\n";
//   CFout << "Output file: " << m_nameOutputFile << "\n";
//   CF_DEBUG_EXIT;

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ComputeDiscreteError::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ComputeDiscreteError::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiscreteError::execute()
{
  CFAUTOTRACE;

  DataHandle<Node*,GLOBAL>  nodes  = socket_nodes.getDataHandle();
  DataHandle<State*,GLOBAL> states = socket_states.getDataHandle();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  RealVector variables(dim);

  bool isUnsteady(false);
  if(SubSystemStatusStack::getActive()->getDT() > 0.) isUnsteady = true;
  if(isUnsteady) variables.resize(dim+1);

  m_L2 = 0.;
  for ( CFuint iState = 0; iState < states.size(); ++iState )
  {

  for (CFuint iDim = 0; iDim < dim ;++iDim){
        variables[iDim] = states[iState]->getCoordinates()[iDim];
      }
  if(isUnsteady) variables[states[iState]->getCoordinates().size()] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

  m_vFunction.evaluate(variables,m_exact);

  m_error = *states[iState] - m_exact;

  for(CFuint iEq =0; iEq < m_L2.size(); iEq ++)
    m_L2[iEq] += m_error[iEq]*m_error[iEq];

  }

  m_L2 /= states.size();

  for(CFuint iEq =0; iEq < m_L2.size(); iEq ++){
    m_L2[iEq] = std::sqrt(m_L2[iEq]);
    m_result[iEq] = std::log(m_L2[iEq]);
  }

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  openFile(true);

  cf_assert(m_file->isopen());
  ofstream& out = m_file->get();

  out << iter << " " ;
  for(CFuint iEq =0; iEq < m_L2.size(); iEq ++){
    out << m_result[iEq] << " "<< m_L2[iEq] << " ";
  }
  out <<"\n";
/*

for(CFuint iEq =0; iEq < m_L2.size(); iEq ++){
     cout <<  m_result[iEq] << "\t" << m_L2[iEq] << "\t";
} */
  out <<"\n";

  m_file->close();
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiscreteError::openFile(bool append)
{
  CFAUTOTRACE;

  cf_assert(!m_file->isopen());

  boost::filesystem::path fpath = Environment::DirPaths::getInstance().getResultsDir() /
boost::filesystem::path(m_nameOutputFile);

#ifdef CF_HAVE_BOOST_1_85
  fpath.replace_extension(".plt");
#else
  boost::filesystem::change_extension(fpath,".plt");
#endif
  fpath = Framework::PathAppender::getInstance().appendParallel( fpath );

  if (append)
  {
    m_file->open(fpath,ios::app);
  }
  else
  {
    m_file->open(fpath);
  }

  m_file->get().precision(16);
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiscreteError::prepareOutputFile()
{
  CFAUTOTRACE;

  openFile(false);

  cf_assert(m_file->isopen());
  ofstream& out = m_file->get();

  out << "TITLE  =  Analytical Error Norm" << "\n";
  out << "VARIABLES = Iter ";
  for (CFuint iEq = 0; iEq < m_L2.size(); ++iEq)
    out << "LogL2_"<<iEq<< " L2_" <<iEq<< " ";

  out <<"\n";

  m_file->close();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace AnalyticalEE

} // namespace COOLFluiD
