#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"

#include "FluxReconstructionNavierStokes/FluxReconstructionNavierStokes.hh"
#include "FluxReconstructionNavierStokes/NS2DAxiSourceTerm.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NS2DAxiSourceTerm, FluxReconstructionSolverData, FluxReconstructionNavierStokesModule>
NS2DAxiSourceTermProvider("NS2DAxiSourceTerm");

//////////////////////////////////////////////////////////////////////////////

void NS2DAxiSourceTerm::defineConfigOptions(Config::OptionList& options)
{
    options.template addConfigOption< CFreal >("CutoffR", "Value of R for which 1/R is cut-off to a maximum.");
}

//////////////////////////////////////////////////////////////////////////////

NS2DAxiSourceTerm::NS2DAxiSourceTerm(const std::string& name) :
    StdSourceTerm(name),
    socket_gradients("gradients"),
    m_srcTerm(),
    m_dim(),
    m_eulerVarSet(CFNULL),
    m_diffVarSet(CFNULL),
    m_solPhysData(),
    m_dummyGradients(),
    m_cellGrads(),
    m_nbrSolPnts()
{
  addConfigOptionsTo(this);
  
  m_cutoffR = 0.0;
  setParameter("CutoffR",&m_cutoffR);
}

//////////////////////////////////////////////////////////////////////////////

NS2DAxiSourceTerm::~NS2DAxiSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void NS2DAxiSourceTerm::getSourceTermData()
{
  StdSourceTerm::getSourceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void NS2DAxiSourceTerm::addSourceTerm(RealVector& resUpdates)
{ 
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // set gradients
  for (CFuint iState = 0; iState < m_nbrSolPnts; ++iState)
  {
    const CFuint stateID = (*m_cellStates)[iState]->getLocalID();
    
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      *((*(m_cellGrads[iState]))[iVar]) = gradients[stateID][iVar];
    }
  }
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  { 
    m_eulerVarSet->computePhysicalData(*((*m_cellStates)[iSol]), m_solPhysData);
    
//    CFreal dUdX = 0.0;
//    CFreal dVdR = 0.0;
//    
//    const CFuint uID = 1;
//    const CFuint vID = 2;
//    
////    if (true) //Puvt
////    {
//      dUdX = (*((*(m_cellGrads[iSol]))[uID]))[XX];
//      dVdR = (*((*(m_cellGrads[iSol]))[vID]))[YY];
////    }
////    else
////    {
////      const CFreal invRho = 1.0/(*((*m_cellStates)[iSol]))[0];
////      const CFreal u = invRho*(*((*m_cellStates)[iSol]))[uID];
////      const CFreal v = invRho*(*((*m_cellStates)[iSol]))[vID];
////      
////      // apply chain rule
////      dUdX = invRho*((*(m_cellGrads[iSol]))[uID][XX] - u*(*(m_cellGrads[iSol]))[0][XX]);
////      dVdR = invRho*((*(m_cellGrads[iSol]))[vID][YY] - v*(*(m_cellGrads[iSol]))[0][YY]);
////    }
//    
//    const CFreal avV = m_solPhysData[EulerTerm::VY];
//    
//    // @todo this will not work if gradients are needed (Menter SST turb model)
//    const CFreal mu = m_diffVarSet->getDynViscosity(*((*m_cellStates)[iSol]), *(m_cellGrads[iSol]));//m_dummyGradients
//    const CFreal coeffMu = m_diffVarSet->getModel().getCoeffTau()*2.0/3.0*mu;
//    const CFreal invR = max(1.0/m_cutoffR,1.0/(m_cell->computeCoordFromMappedCoord((*m_solPntsLocalCoords)[iSol]))[YY]);
//    const CFreal tauThetaTheta = -coeffMu*(dUdX + dVdR - 2.0*avV*invR);
//    
//    // AL: check this in Hontzatko's report (dp!)
//    //m_srcTerm[vID] = m_solPhysData[EulerTerm::P] - tauThetaTheta;
//    resUpdates[m_nbrEqs*iSol + vID] = m_solPhysData[EulerTerm::P] - tauThetaTheta;
 
    const CFreal r = max(m_cutoffR,(m_cell->computeCoordFromMappedCoord((*m_solPntsLocalCoords)[iSol]))[YY]);
        
    m_diffVarSet->getAxiSourceTerm(m_solPhysData,*((*m_cellStates)[iSol]),*(m_cellGrads[iSol]),r,m_srcTerm);
      
    m_srcTerm /= r;
      
    for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq)
    {
      resUpdates[m_nbrEqs*iSol + iEq] = m_srcTerm[iEq];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NS2DAxiSourceTerm::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  StdSourceTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void NS2DAxiSourceTerm::setup()
{
  CFAUTOTRACE;
  StdSourceTerm::setup();

  // get dimensionality
  m_dim = PhysicalModelStack::getActive()->getDim ();
  
  // get NS 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in Euler2DAxiSourceTerm!");
  }
  cf_assert(m_nbrEqs == 4);
  
  m_eulerVarSet->getModel()->resizePhysicalData(m_solPhysData);
  
  m_diffVarSet = getMethodData().getDiffusiveVar().d_castTo<NavierStokes2DVarSet>();
  
  m_srcTerm.resize(m_nbrEqs);
  
  // get the local spectral FD data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  // number of sol points
  m_nbrSolPnts = frLocalData[0]->getNbrOfSolPnts();
  
  m_cellGrads.resize(m_nbrSolPnts);
  
  for (CFuint iSol = 0; iSol < m_nbrSolPnts; ++iSol)
  {
    m_cellGrads[iSol] = new std::vector< RealVector* >(m_nbrEqs);
    
    for (CFuint iVar = 0; iVar < m_nbrEqs; ++iVar)
    {
      (*(m_cellGrads[iSol]))[iVar] = new RealVector(m_dim);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NS2DAxiSourceTerm::unsetup()
{
  CFAUTOTRACE;
  StdSourceTerm::unsetup();
  
  m_cellGrads.clear();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    NS2DAxiSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();

  result.push_back(&socket_gradients);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
