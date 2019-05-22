#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/NavierStokes2DVarSet.hh"

#include "FluxReconstructionTurb/FluxReconstructionKOmega.hh"
#include "FluxReconstructionTurb/KOmega2DSourceTerm.hh"
#include "KOmega/NavierStokesKOmegaVarSetTypes.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::KOmega;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<KOmega2DSourceTerm, FluxReconstructionSolverData, FluxReconstructionKOmegaModule>
KOmega2DSourceTermProvider("KOmega2DSourceTerm");

//////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

KOmega2DSourceTerm::KOmega2DSourceTerm(const std::string& name) :
    StdSourceTerm(name),
    socket_gradients("gradients"),
    m_dim(),
    m_eulerVarSet(CFNULL),
    m_diffVarSet(CFNULL),
    m_solPhysData(),
    m_dummyGradients(),
    m_cellGrads()
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

KOmega2DSourceTerm::~KOmega2DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::getSourceTermData()
{
  StdSourceTerm::getSourceTermData();
}

//////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::addSourceTerm(RealVector& resUpdates)
{
  m_diffVarSet = getMethodData().getDiffusiveVar();
  
  SafePtr< NavierStokes2DKOmega > navierStokesVarSet = m_diffVarSet.d_castTo< NavierStokes2DKOmega >();

  const CFuint nbrSol = m_cellStates->size();
  const bool Puvt = getMethodData().getUpdateVarStr() == "Puvt";
  
  // get the gradients datahandle
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();

  // set gradients
  const CFuint nbrStates = m_cellStates->size();
  m_cellGrads.resize(nbrStates);
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const CFuint stateID = (*m_cellStates)[iState]->getLocalID();
    m_cellGrads[iState] = &gradients[stateID];
  }
  
  const EquationSubSysDescriptor& eqData = PhysicalModelStack::getActive()->getEquationSubSysDescriptor();
  const CFuint nbEqs = eqData.getNbEqsSS();
  const CFuint totalNbEqs = PhysicalModelStack::getActive()->getNbEq();  
  
  for (CFuint iSol = 0; iSol < nbrSol; ++iSol)
  {
    m_eulerVarSet->computePhysicalData(*((*m_cellStates)[iSol]), m_solPhysData);
    
    const CFreal R = m_eulerVarSet->getModel()->getR();
    
    
    
    
    
    CFreal dUdX = 0.0;
    CFreal dVdR = 0.0;
    
    const CFuint uID = 1;
    const CFuint vID = 2;
    
    if (true) //Puvt
    {
      dUdX = (*(m_cellGrads[iSol]))[uID][XX];
      dVdR = (*(m_cellGrads[iSol]))[vID][YY];
    }
    else
    {
      const CFreal invRho = 1.0/(*((*m_cellStates)[iSol]))[0];
      const CFreal u = invRho*(*((*m_cellStates)[iSol]))[uID];
      const CFreal v = invRho*(*((*m_cellStates)[iSol]))[vID];
      
      // apply chain rule
      dUdX = invRho*((*(m_cellGrads[iSol]))[uID][XX] - u*(*(m_cellGrads[iSol]))[0][XX]);
      dVdR = invRho*((*(m_cellGrads[iSol]))[vID][YY] - v*(*(m_cellGrads[iSol]))[0][YY]);
    }
    
    const CFreal avV = m_solPhysData[EulerTerm::VY];
    
    // @todo this will not work if gradients are needed (Menter SST turb model)
    const CFreal mu = navierStokesVarSet->getDynViscosity(*((*m_cellStates)[iSol]), m_dummyGradients);
    const CFreal coeffMu = navierStokesVarSet->getModel().getCoeffTau()*2.0/3.0*mu;
    const CFreal invR = 1.0/(m_cell->computeCoordFromMappedCoord((*m_solPntsLocalCoords)[iSol]))[YY];
    const CFreal tauThetaTheta = -coeffMu*(dUdX + dVdR - 2.0*avV*invR);
    
    // AL: check this in Hontzatko's report (dp!)
    //m_srcTerm[vID] = m_solPhysData[EulerTerm::P] - tauThetaTheta;
    resUpdates[m_nbrEqs*iSol + vID] = m_solPhysData[EulerTerm::P] - tauThetaTheta;

//     for (CFuint iEq = 0; iEq < m_nbrEqs; ++iEq, ++resID)
//     {
//       rhs[resID] += resFactor*m_solPntJacobDets[iSol]*m_srcTerm[iEq];
//     }
  }
}

//////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  StdSourceTerm::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::setup()
{
  CFAUTOTRACE;
  StdSourceTerm::setup();

  // get dimensionality
  m_dim = PhysicalModelStack::getActive()->getDim ();
  
  // get NS 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in KOmega2DSourceTerm!");
  }
  cf_assert(m_nbrEqs == 4);
  
  m_eulerVarSet->getModel()->resizePhysicalData(m_solPhysData);
}

//////////////////////////////////////////////////////////////////////////////

void KOmega2DSourceTerm::unsetup()
{
  CFAUTOTRACE;
  StdSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    KOmega2DSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();

  result.push_back(&socket_gradients);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD
