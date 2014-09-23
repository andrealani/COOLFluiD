#include "Environment/DirPaths.hh"

#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/PathAppender.hh"

#include "AnalyticalEE/AnalyticalEE.hh"
#include "AnalyticalEE/ComputeIntegralError.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AnalyticalEE {

//////////////////////////////////////////////////////////////////////////////

void ComputeIntegralError::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Functions","Definition of the Functions.");
  options.addConfigOption< std::string >("OutputFile","Name of Output File to write the Error Norm.");
}

//////////////////////////////////////////////////////////////////////////////

// CFreal ComputeIntegralError::qdWeightsTriag[4] = { -27./48., 25./48., 25./48., 25./48. };
// CFreal ComputeIntegralError::xQdPtsTriag[4] = { 1./3., 0.2, 0.6, 0.2 };
// CFreal ComputeIntegralError::yQdPtsTriag[4] = { 1./3., 0.2, 0.2, 0.6 };
// CFuint ComputeIntegralError::nbQdPts = 4;

CFreal ComputeIntegralError::a1 = 0.797426985353;
CFreal ComputeIntegralError::b1 = 0.101286507323;
CFreal ComputeIntegralError::a2 = 0.059715871789;
CFreal ComputeIntegralError::b2 = 0.470142064105;
CFreal ComputeIntegralError::w1 = 0.225;
CFreal ComputeIntegralError::w2 = 0.125939180544;
CFreal ComputeIntegralError::w3 = 0.132394152788;



CFreal ComputeIntegralError::qdWeightsTriag[7] = { w1, w2, w2, w2, w3, w3, w3 };
CFreal ComputeIntegralError::xQdPtsTriag[7] = { 1./3., b1, a1, b1, b2, a2, b2 };
CFreal ComputeIntegralError::yQdPtsTriag[7] = { 1./3., b1, b1, a1, b2, b2, a2 };
CFuint ComputeIntegralError::nbQdPts = 7;




//////////////////////////////////////////////////////////////////////////////

ComputeIntegralError::ComputeIntegralError(std::string name) : AnalyticEECom(name),
  socket_nodes("nodes"),
  socket_states("states"),
  m_dim(DIM_2D),
  m_nbEqs(1),
  m_mappedCoord(DIM_2D),
  m_physicalCoord(DIM_2D),
  m_numerical(0),
  m_exact(0),
  m_L2(0),
  m_result(0),
  m_qdState(0),
  m_stdTrsGeoBuilder()
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeIntegralError::~ComputeIntegralError()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeIntegralError::setup()
{
  CFAUTOTRACE;
  // first call parent method
  AnalyticEECom::setup();

  m_dim = PhysicalModelStack::getActive()->getDim();
  m_nbEqs = PhysicalModelStack::getActive()->getNbEq();

  m_mappedCoord.resize(m_dim);
  m_physicalCoord.resize(m_dim);
  m_numerical.resize(m_nbEqs);
  m_exact.resize(m_nbEqs);
  m_error.resize(m_nbEqs);

  m_L2.resize(m_nbEqs);
  m_result.resize(m_nbEqs);

  prepareOutputFile();

  m_stdTrsGeoBuilder.setup();

}

//////////////////////////////////////////////////////////////////////////////

void ComputeIntegralError::unsetup()
{
  CFAUTOTRACE;

  m_numerical.resize(0);
  m_exact.resize(0);
  m_error.resize(0);
  m_L2.resize(0);
  m_result.resize(0);
  m_file->close();

  // last call parent method
  AnalyticEECom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void ComputeIntegralError::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  AnalyticEECom::configure(args);

  m_vFunction.setFunctions(m_functions);

  m_vFunction.setVariables(m_vars);
  try {
    m_vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
ComputeIntegralError::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = AnalyticEECom::providesSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ComputeIntegralError::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = AnalyticEECom::needsSockets();

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeIntegralError::execute()
{ }

//////////////////////////////////////////////////////////////////////////////

void ComputeIntegralError::openFile(bool append)
{
  CFAUTOTRACE;

  cf_assert(!m_file->isopen());

  boost::filesystem::path fpath = Environment::DirPaths::getInstance().getResultsDir() / 
boost::filesystem::path(m_nameOutputFile);
  boost::filesystem::change_extension(fpath,".plt");
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

void ComputeIntegralError::prepareOutputFile()
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
