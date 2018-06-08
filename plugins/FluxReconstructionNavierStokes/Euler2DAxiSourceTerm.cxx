#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/Euler2DVarSet.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/Euler2DAxiSourceTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Euler2DAxiSourceTerm, FluxReconstructionSolverData, FluxReconstructionNavierStokesModule>
Euler2DAxiSourceTermProvider("Euler2DAxiSourceTerm");

//////////////////////////////////////////////////////////////////////////////

void Euler2DAxiSourceTerm::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DAxiSourceTerm::Euler2DAxiSourceTerm(const std::string& name) :
    StdSourceTerm(name),
    m_srcTerm(),
    m_dim(),
    m_eulerVarSet(CFNULL),
    m_solPhysData()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DAxiSourceTerm::~Euler2DAxiSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DAxiSourceTerm::getSourceTermData()
{
  StdSourceTerm::getSourceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DAxiSourceTerm::addSourceTerm(RealVector& resUpdates)
{
//   // get the datahandle of the rhs
//   DataHandle< CFreal > rhs = socket_rhs.getDataHandle();
// 
//   // get residual factor
//   const CFreal resFactor = getMethodData().getResFactor();

  // loop over solution points in this cell to add the source term
  CFuint resID = m_nbrEqs*( (*m_cellStates)[0]->getLocalID() );
  const CFuint nbrSol = m_cellStates->size();
  for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
  { 
    m_eulerVarSet->computePhysicalData(*((*m_cellStates)[iSol]), m_solPhysData);
    
    resUpdates[m_nbrEqs*iSol+2] = m_solPhysData[EulerTerm::P];

//     rhs[resID+2] += resFactor*m_solPntJacobDets[iSol]*m_solPhysData[EulerTerm::P];
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DAxiSourceTerm::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  StdSourceTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DAxiSourceTerm::setup()
{
  CFAUTOTRACE;
  StdSourceTerm::setup();

  // get dimensionality
  m_dim = PhysicalModelStack::getActive()->getDim ();

  // resize m_srcTerm
  m_srcTerm.resize(m_nbrEqs);
  
  // get Euler 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in Euler2DAxiSourceTerm!");
  }
  cf_assert(m_nbrEqs == 4);
  
  m_eulerVarSet->getModel()->resizePhysicalData(m_solPhysData);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DAxiSourceTerm::unsetup()
{
  CFAUTOTRACE;
  StdSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    Euler2DAxiSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
