#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "FluxReconstructionMethod/AnalyticalSourceTerm.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<AnalyticalSourceTerm, FluxReconstructionSolverData, FluxReconstructionModule>
AnalyticalSourceTermProvider("AnalyticalSourceTerm");

//////////////////////////////////////////////////////////////////////////////

void AnalyticalSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

AnalyticalSourceTerm::AnalyticalSourceTerm(const std::string& name) :
    StdSourceTerm(name),
    m_vFunctionInputVars(),
    m_srcTerm(),
    m_dim()
{
  addConfigOptionsTo(this);
  m_functions = vector<std::string>();
  setParameter("Def",&m_functions);

  m_vars = vector<std::string>();
  setParameter("Vars",&m_vars);
}

//////////////////////////////////////////////////////////////////////////////

AnalyticalSourceTerm::~AnalyticalSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void AnalyticalSourceTerm::getSourceTermData()
{
  StdSourceTerm::getSourceTermData();

  // get current time
  m_vFunctionInputVars[m_dim] = SubSystemStatusStack::getActive()->getCurrentTimeDim();
}

//////////////////////////////////////////////////////////////////////////////

void AnalyticalSourceTerm::addSourceTerm(RealVector& resUpdates)
{
//   // get the datahandle of the rhs
//   DataHandle< CFreal > rhs = socket_rhs.getDataHandle();
// 
//   // get residual factor
//   const CFreal resFactor = getMethodData().getResFactor();
// 
//   // loop over solution points in this cell to add the source term
//   CFuint resID = m_nbrEqs*( (*m_cellStates)[0]->getLocalID() );
  const CFuint nbrSol = m_cellStates->size();
  for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
  {
    // get solution point coordinates and reference length
    const Node& node = (*m_cellStates)[iSol]->getCoordinates();
    const CFreal refLength = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

    // copy coordinates
    for (CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      m_vFunctionInputVars[iDim] = node[iDim]*refLength;
    }

    // current time is put in getSourceTermData

    // evaluate the function at the state coordinate
    m_vFunction.evaluate(m_vFunctionInputVars,m_srcTerm);

    // loop over physical variables
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
//       rhs[resID] += resFactor*m_solPntJacobDets[iSol]*m_srcTerm[iEq];
      resUpdates[m_nbrEqs*iSol + iEq] = m_srcTerm[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void AnalyticalSourceTerm::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  StdSourceTerm::configure(args);

  // parsing the functions that the user inputed
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

void AnalyticalSourceTerm::setup()
{
  CFAUTOTRACE;
  
  // setup the parent class
  StdSourceTerm::setup();

  // get dimensionality
  m_dim = PhysicalModelStack::getActive()->getDim ();

  // resize m_vFunctionInputVars
  m_vFunctionInputVars.resize(m_dim+1);

  // resize m_srcTerm
  m_srcTerm.resize(m_nbrEqs);
}

//////////////////////////////////////////////////////////////////////////////

void AnalyticalSourceTerm::unsetup()
{
  CFAUTOTRACE;
  
  // unsetup the parent class
  StdSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    AnalyticalSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
